import pandas as pd
import shared.display as dis
import numpy as np
import os

data_source = "D:/Xiaowei/data/20210607_screen/summary/NT/"
save_source = "D:/Xiaowei/data/20210607_screen/volcanoPlots/NT/"
WT_source = "D:/Xiaowei/data/20210607_screen/WTFiles/NT/"

analyze_organelle = 'nucleoli'
analysis_mode = 'single_exp'
features = ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope', 'size_organelle', 'raw_int_organelle',
            'circ_organelle', 'ecce_organelle']

data_WT = pd.read_csv(("%sWT_data_full.txt" % WT_source), na_values=['.'], sep='\t')
data_WT_ft = data_WT[data_WT['frap_filter_%s' % analysis_mode] == 1]
data_WT_organelle = pd.read_csv(("%sWT_data_%s.txt" % (WT_source, analyze_organelle)), na_values=['.'], sep='\t')
centers = [np.mean(data_WT_ft['%s_mobile_fraction' % analysis_mode]), np.mean(data_WT_ft['mobile_fraction']),
           np.mean(data_WT_ft['%s_t_half' % analysis_mode]), np.mean(data_WT_ft['t_half']),
           np.mean(data_WT_ft['%s_slope' % analysis_mode]), np.mean(data_WT_ft['linear_slope']),
           np.mean(data_WT_organelle['size']), np.mean(data_WT_organelle['raw_int']),
           np.mean(data_WT_organelle[data_WT_organelle['size'] > 50]['circ']),
           np.mean(data_WT_organelle['eccentricity'])]

data_p = pd.read_csv(("%ssummary_p.txt" % data_source), na_values=['.'], sep='\t')
data_value = pd.read_csv(("%ssummary_value.txt" % data_source), na_values=['.'], sep='\t')
data_gamma = pd.read_csv(("%ssummary_gamma.txt" % data_source), na_values=['.'], sep='\t')

if not os.path.exists(save_source):
    os.makedirs(save_source)

# export volcano hit plots
for i in range(len(features)):
    dis.plot_volcano_hit(data_p, data_value, data_gamma, features[i], centers[i], 60, 'Y', save_source)

print("DONE!")