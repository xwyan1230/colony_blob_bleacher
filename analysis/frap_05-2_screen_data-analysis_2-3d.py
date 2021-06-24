import os
import pandas as pd
import shared.display as dis
import numpy as np

multi_data_source = "D:/Xiaowei/data/20210607_screen/dataFiles/"
WT_source = "D:/Xiaowei/data/20210607_screen/WTFiles/NT/"
save_source = "D:/Xiaowei/data/20210607_screen/dataPlots/NT/"

analyze_organelle = 'nucleoli'  # only accepts 'sg' or 'nucleoli'
analysis_mode = 'single_exp'
export_figure_mode = 'on'

bounds_mob = np.array([0, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1, 1.2])
bounds_t_half = np.array([0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0])
bounds_slope = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 1.0])

multi_dirs = [x for x in os.listdir(multi_data_source)]
if '.DS_Store' in multi_dirs:
    multi_dirs.remove('.DS_Store')

data_name = []
data_frap_n_curve = []
mob = []
curve_mob = []
t_half = []
curve_t_half = []
slope = []
curve_slope = []
data_organelle_n = []
size = []
raw_int = []
data_organelle_n_circ = []
circ = []
ecce = []

for r in range(len(multi_dirs)):
    name = multi_dirs[r]
    data_name.append(name)
    print('# Analyzing %s ... (%d/%d)' % (name, r+1, len(multi_dirs)))
    data_source = ("%s%s/" % (multi_data_source, name))
    save_path = ("%s%s/" % (save_source, name))

    data_sample = pd.read_csv(("%s%s_data_full.txt" % (data_source, name)), na_values=['.'], sep='\t')
    data_sample_ft = data_sample[data_sample['frap_filter_%s' % analysis_mode] == 1]
    data_frap_n_curve.append(len(data_sample_ft))

    data_WT = pd.read_csv(("%sWT_data_full.txt" % WT_source), na_values=['.'], sep='\t')
    data_WT_ft = data_WT[data_WT['frap_filter_%s' % analysis_mode] == 1]

    # save average data
    mob.append(np.mean(data_sample_ft['%s_mobile_fraction' % analysis_mode]))
    curve_mob.append(np.mean(data_sample_ft['mobile_fraction']))
    t_half.append(np.mean(data_sample_ft['%s_t_half' % analysis_mode]))
    curve_t_half.append(np.mean(data_sample_ft['t_half']))
    slope.append(np.mean(data_sample_ft['%s_slope' % analysis_mode]))
    curve_slope.append(np.mean(data_sample_ft['linear_slope']))

    data_sample_organelle = pd.read_csv(("%s%s_data_%s.txt" % (data_source, name, analyze_organelle)), na_values=['.'],
                                        sep='\t')
    data_WT_organelle = pd.read_csv(("%sWT_data_%s.txt" % (WT_source, analyze_organelle)), na_values=['.'], sep='\t')

    data_organelle_n.append(len(data_sample_organelle))
    data_organelle_n_circ.append(len(data_sample_organelle[data_sample_organelle['size'] > 50]))

    # save average data
    size.append(np.mean(data_sample_organelle['size']))
    raw_int.append(np.mean(data_sample_organelle['raw_int']))
    circ.append(np.mean(data_sample_organelle[data_sample_organelle['size'] > 50]['circ']))
    ecce.append(np.mean(data_sample_organelle['eccentricity']))

    if export_figure_mode == 'on':
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        # export FRAP 3d images
        print("# Export mobile_fraction 3d plot ...")
        dis.plot_frap_3d(data_WT_ft, data_sample_ft, 'mobile_fraction', bounds_mob, name, 'curve_mob', save_path)
        dis.plot_frap_3d(data_WT_ft, data_sample_ft, '%s_mobile_fraction' % analysis_mode, bounds_mob, name, 'mob',
                         save_path)
        print("# Export t_half 3d plot ...")
        dis.plot_frap_3d(data_WT_ft, data_sample_ft, 't_half', bounds_t_half, name, 'curve_t_half', save_path)
        dis.plot_frap_3d(data_WT_ft, data_sample_ft, '%s_t_half' % analysis_mode, bounds_t_half, name, 't_half',
                         save_path)
        print("# Export slope 3d plot ...")
        dis.plot_frap_3d(data_WT_ft, data_sample_ft, 'linear_slope', bounds_slope, name, 'curve_slope', save_path)
        dis.plot_frap_3d(data_WT_ft, data_sample_ft, '%s_slope' % analysis_mode, bounds_slope, name, 'slope',
                         save_path)

        # export organelle 2d images
        print("# Export organelle 2d images ...")
        dis.plot_organelle_2d(data_WT_organelle, data_sample_organelle, 'raw_int', 'log', 'ln(raw intensity)',
                              'size', 'linear', 'size (pixel)', name, save_path)
        dis.plot_organelle_2d(data_WT_organelle, data_sample_organelle, 'raw_int', 'log', 'ln(raw intensity)',
                              'circ', 'linear', 'circularity', name, save_path)
        dis.plot_organelle_2d(data_WT_organelle, data_sample_organelle, 'raw_int', 'log', 'ln(raw intensity)',
                              'eccentricity', 'linear', 'eccentricity', name, save_path)
        dis.plot_organelle_2d(data_WT_organelle, data_sample_organelle, 'size', 'linear', 'size (pixel)',
                              'eccentricity', 'linear', 'eccentricity', name, save_path)
        dis.plot_organelle_2d(data_WT_organelle, data_sample_organelle, 'circ', 'linear', 'circularity',
                              'eccentricity', 'linear', 'eccentricity', name, save_path)
        dis.plot_organelle_2d(data_WT_organelle, data_sample_organelle, 'size', 'linear', 'size (pixel)',
                              'circ', 'linear', 'circularity', name, save_path)

data_frap = pd.DataFrame({'sample': data_name, 'n_curve': data_frap_n_curve, 'mob': mob, 'curve_mob': curve_mob,
                          't_half': t_half, 'curve_t_half': curve_t_half, 'slope': slope, 'curve_slope': curve_slope,
                          'n_organelle': data_organelle_n, 'size_organelle': size, 'raw_int_organelle': raw_int,
                          'n_organelle_circ': data_organelle_n_circ, 'circ_organelle': circ, 'ecce_organelle': ecce})

save_path = ("%ssummary/" % save_source)
if not os.path.exists(save_path):
    os.makedirs(save_path)

data_frap.to_csv('%ssummary_value.txt' % save_path, index=False, sep='\t')

print("DONE!")
