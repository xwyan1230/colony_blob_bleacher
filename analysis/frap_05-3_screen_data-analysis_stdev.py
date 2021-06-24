import os
import pandas as pd
import numpy as np
from datetime import datetime

master_folder = "D:/Xiaowei/data/20210607_screen/"
multi_data_source = "D:/Xiaowei/data/20210607_screen/dataFiles/"
save_source = "D:/Xiaowei/data/20210607_screen/dataStdev/NT/"

analyze_organelle = 'nucleoli'  # only accepts 'sg' or 'nucleoli'
analysis_mode = 'single_exp'

# log all the running info
if not os.path.exists("%sscreen_processing_log.txt" % master_folder):
    f = open("%sscreen_processing_log.txt" % master_folder, "w+")
else:
    f = open("%sscreen_processing_log.txt" % master_folder, "a+")

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
f.write("datetime: %s\n" % dt_string)
f.write("code_running: %s\n" % __file__)
f.write("master_folder: %s\n" % master_folder)
f.write("multi_data_source: %s\n" % multi_data_source)
f.write("save_source: %s\n" % save_source)
f.write("analyze_organelle: %s\n" % analyze_organelle)
f.write("analysis_mode: %s\n\n" % analysis_mode)

f.close()

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

    # save stdev data
    mob.append(np.std(data_sample_ft['%s_mobile_fraction' % analysis_mode]))
    curve_mob.append(np.std(data_sample_ft['mobile_fraction']))
    t_half.append(np.std(data_sample_ft['%s_t_half' % analysis_mode]))
    curve_t_half.append(np.std(data_sample_ft['t_half']))
    slope.append(np.std(data_sample_ft['%s_slope' % analysis_mode]))
    curve_slope.append(np.std(data_sample_ft['linear_slope']))

    data_sample_organelle = pd.read_csv(("%s%s_data_%s.txt" % (data_source, name, analyze_organelle)), na_values=['.'],
                                        sep='\t')
    data_organelle_n.append(len(data_sample_organelle))
    data_organelle_n_circ.append(len(data_sample_organelle[data_sample_organelle['size'] > 50]))

    # save stdev data
    size.append(np.std(data_sample_organelle['size']))
    raw_int.append(np.std(data_sample_organelle['raw_int']))
    circ.append(np.std(data_sample_organelle[data_sample_organelle['size'] > 50]['circ']))
    ecce.append(np.std(data_sample_organelle['eccentricity']))

data_frap = pd.DataFrame({'sample': data_name, 'n_curve': data_frap_n_curve, 'mob': mob, 'curve_mob': curve_mob,
                          't_half': t_half, 'curve_t_half': curve_t_half, 'slope': slope, 'curve_slope': curve_slope,
                          'n_organelle': data_organelle_n, 'size_organelle': size, 'raw_int_organelle': raw_int,
                          'n_organelle_circ': data_organelle_n_circ, 'circ_organelle': circ, 'ecce_organelle': ecce})

save_path = ("%ssummary/" % save_source)
if not os.path.exists(save_path):
    os.makedirs(save_path)

data_frap.to_csv('%ssummary_stdev.txt' % save_path, index=False, sep='\t')

print("DONE!")
