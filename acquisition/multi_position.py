import numpy as np
import random
import time
from pycromanager import Bridge
from skimage.measure import label, regionprops
from skimage import morphology

from shared.analysis import central_pixel_without_cells, bleach_location
from shared.find_blobs import select

# variables
from shared.find_organelles import find_organelle

nr = 40
nr_between_projector_checks = 2
cal_exposure = 200
cal_offset = 5
n_curve = 500

# build up pycromanager bridge
bridge = Bridge()
mmc = bridge.get_core()
mm = bridge.get_studio()
projector = bridge.construct_java_object("org.micromanager.projector.ProjectorAPI")
projector_device = projector.get_projection_device()


def snap_and_get_bleach_location(exposure, cutoff):
    """
    Takes an image with the current settings.  Finds a location close to the center where there are
    no objects (as defined in function central_picel_without_cells).  If no such location is found, returns -1.
    Targets bleacher to this location, exposes and takes an image of that exposure.  Finds the center of the
    actual bleach spot.  When the square of the distance between the intended target and the actual bleach spot is
    greater than provided offset, will execute a full calibration.
    :param exposure: exposure time to use for bleaching
    :param cutoff: square of distance.  When offset is higher than this code should execute a calibration
    :return: tuple with first Boolean indicating if a calibration took place, second variable the square of the offset distance
    """
    p_exposure = projector_device.get_exposure()
    c_exposure = mmc.get_exposure()
    test_img = mm.live().snap(True).get(0)
    test_np_img = np.reshape(test_img.get_raw_pixels(), newshape=[test_img.get_height(), test_img.get_width()])
    location = central_pixel_without_cells(test_np_img)
    if location:
        auto_shutter = mm.shutter().get_auto_shutter()
        mm.shutter().set_auto_shutter(False)
        projector.set_exposure(projector_device, exposure)
        mmc.set_exposure(exposure)
        projector.enable_point_and_shoot_mode(True)
        pre_img = mm.live().snap(True).get(0)
        pre_np_img = np.reshape(pre_img.get_raw_pixels(), newshape=[pre_img.get_height(), pre_img.get_width()])
        projector.add_point_to_point_and_shoot_queue(int(location[1]), int(location[0]))
        post_img = mm.live().snap(True).get(0)
        post_np_img = np.reshape(post_img.get_raw_pixels(), newshape=[post_img.get_height(), post_img.get_width()])
        measured_location = bleach_location(pre_np_img, post_np_img, location, [100, 100])
        offset = (measured_location[0] - location[0], measured_location[1] - location[1])
        print(offset)
        cal = False
        if offset[0] * offset[0] + offset[1] * offset[1] > cutoff:
            projector.calibrate(True)
            cal = True
            print("Calibrated")
        projector.set_exposure(projector_device, p_exposure)
        mmc.set_exposure(c_exposure)
        mm.shutter().set_auto_shutter(auto_shutter)
        return cal, offset[0] * offset[0] + offset[1] * offset[1]
    return False, -1


# TODO We may want to configure the acquisition settings to ensure they are what we want

pm = mm.positions()
pos_list = pm.get_position_list()
well = pos_list.get_position(0).get_label().split('-')[0]
well_count = 0
ds = mm.data().create_ram_datastore()
count = 0

for idx in range(pos_list.get_number_of_positions()):
    pos = pos_list.get_position(idx)

    well_temp = pos.get_label().split('-')[0]
    if well_temp == well:
        if well_count >= n_curve:
            continue
    else:
        well_count = 0
        well = well_temp

    # Close DataViewer opened during previous run
    dv = mm.displays().close_displays_for(ds)
    pos.go_to_position(pos, mmc)

    time.sleep(0.1)
    if count >= nr_between_projector_checks:
        calibrated, error = snap_and_get_bleach_location(cal_exposure, cal_offset)
        if error < 0:
            count -= 1
        else:
            count = 0
        if calibrated:
            continue
    count += 1
    img = mm.live().snap(False).get(0)
    pixels = np.reshape(img.get_raw_pixels(), newshape=[img.get_height(), img.get_width()])
    # find organelles using a combination of thresholding and watershed
    _, segmented = find_organelle(pixels, 'local-nucleoli', 500, 200, 10, 1000)
    label_img = label(segmented)
    label_img = morphology.remove_small_objects(label_img, 5)
    blobs = regionprops(label_img)
    centered = select(blobs, 'centroid', img.get_width() / 10, 0.9 * img.get_width())

    if len(centered) > (nr // 2):
        projector.enable_point_and_shoot_mode(True)
        ssb = mm.acquisitions().get_acquisition_settings().copy_builder()
        mm.acquisitions().set_acquisition_settings(ssb.prefix(pos.get_label()).build())
        ds = mm.acquisitions().run_acquisition_nonblocking()
        # Trick to get timing right.  Wait for Core to report that Sequence is running
        while not mmc.is_sequence_running(mmc.get_camera_device()):
            time.sleep(0.1)
        time.sleep(1.5)

        for region_list in [centered]:
            nr_shots = nr if len(region_list) >= (2 * nr) else int(len(region_list) / 2)
            well_count += nr_shots
            shots = random.sample(region_list, nr_shots)
            # shots = region_list[0:10]
            for shot in shots:
                # Note that MM has x-y coordinates, and Python uses row-column (equivalent to y-x)
                projector.add_point_to_point_and_shoot_queue(shot['centroid'][1], shot['centroid'][0])
                # print(shot['centroid'][1], " ", shot['centroid'][0])
                time.sleep(0.07)
            print(pos.get_label(), ": Shots ", len(shots))
            while mmc.is_sequence_running(mmc.get_camera_device()):
                time.sleep(0.5)
            time.sleep(1)

print("Done!")