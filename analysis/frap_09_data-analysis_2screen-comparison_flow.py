import pandas as pd
import shared.display as dis
import shared.dataframe as dat
import numpy as np
import os

data1_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210610_CBB_nucleoliFRAPscreenComparison/summaryFiles/4/summary/NT/"
data2_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210610_CBB_nucleoliFRAPscreenComparison/summaryFiles/6-1/summary/NT/"
save_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210610_CBB_nucleoliFRAPscreenComparison/comparison/4_vs_6-1/"
WT1_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210610_CBB_nucleoliFRAPscreenComparison/summaryFiles/4/WTFiles/NT/"
WT2_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210610_CBB_nucleoliFRAPscreenComparison/summaryFiles/6-1/WTFiles/NT/"

copy_mode = 'table'  # accepts 'sample' or 'table'
features = ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope', 'size_organelle', 'raw_int_organelle',
            'circ_organelle', 'ecce_organelle']
analysis_mode = 'single_exp'
analyze_organelle = 'nucleoli'

if not os.path.exists(save_source):
    os.makedirs(save_source)

# get centers
data1_WT = pd.read_csv(("%sWT_data_full.txt" % WT1_source), na_values=['.'], sep='\t')
data1_WT_ft = data1_WT[data1_WT['frap_filter_%s' % analysis_mode] == 1]
data1_WT_organelle = pd.read_csv(("%sWT_data_%s.txt" % (WT1_source, analyze_organelle)), na_values=['.'], sep='\t')
centers1 = [np.mean(data1_WT_ft['%s_mobile_fraction' % analysis_mode]), np.mean(data1_WT_ft['mobile_fraction']),
            np.mean(data1_WT_ft['%s_t_half' % analysis_mode]), np.mean(data1_WT_ft['t_half']),
            np.mean(data1_WT_ft['%s_slope' % analysis_mode]), np.mean(data1_WT_ft['linear_slope']),
            np.mean(data1_WT_organelle['size']), np.mean(data1_WT_organelle['raw_int']),
            np.mean(data1_WT_organelle[data1_WT_organelle['size'] > 50]['circ']),
            np.mean(data1_WT_organelle['eccentricity'])]

data2_WT = pd.read_csv(("%sWT_data_full.txt" % WT2_source), na_values=['.'], sep='\t')
data2_WT_ft = data2_WT[data2_WT['frap_filter_%s' % analysis_mode] == 1]
data2_WT_organelle = pd.read_csv(("%sWT_data_%s.txt" % (WT2_source, analyze_organelle)), na_values=['.'], sep='\t')
centers2 = [np.mean(data2_WT_ft['%s_mobile_fraction' % analysis_mode]), np.mean(data2_WT_ft['mobile_fraction']),
            np.mean(data2_WT_ft['%s_t_half' % analysis_mode]), np.mean(data2_WT_ft['t_half']),
            np.mean(data2_WT_ft['%s_slope' % analysis_mode]), np.mean(data2_WT_ft['linear_slope']),
            np.mean(data2_WT_organelle['size']), np.mean(data2_WT_organelle['raw_int']),
            np.mean(data2_WT_organelle[data2_WT_organelle['size'] > 50]['circ']),
            np.mean(data2_WT_organelle['eccentricity'])]

# combine value data
data1_value = pd.read_csv(("%ssummary_value.txt" % data1_source), na_values=['.'], sep='\t')
data2_value = pd.read_csv(("%ssummary_value.txt" % data2_source), na_values=['.'], sep='\t')

if copy_mode == 'sample':
    data2_value['sample_matched'] = data2_value['sample']
elif copy_mode == 'table':
    table = pd.read_csv(("%s2screenSampleTable.txt" % save_source), na_values=['.'], sep='\t')
    data2_value = dat.copy_based_on_index(data2_value, table, 'sample', 'screen2', ['sample_matched'], ['screen1'])

data_value = data1_value.copy()
data_value = dat.copy_based_on_index(data_value, data2_value, 'sample', 'sample_matched',
                                     ['mob1', 'curve_mob1', 't_half1', 'curve_t_half1', 'slope1', 'curve_slope1',
                                      'size_organelle1', 'raw_int_organelle1', 'circ_organelle1', 'ecce_organelle1'],
                                     ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope',
                                      'size_organelle', 'raw_int_organelle', 'circ_organelle', 'ecce_organelle'])
data_value = data_value.dropna()

for i in range(len(features)):
    data_value['%s_centered' % features[i]] = data_value[features[i]] - centers1[i]
    data_value['%s1_centered' % features[i]] = data_value['%s1' % features[i]] - centers2[i]

# combine p data
data1_p = pd.read_csv(("%ssummary_p.txt" % data1_source), na_values=['.'], sep='\t')
data2_p = pd.read_csv(("%ssummary_p.txt" % data2_source), na_values=['.'], sep='\t')

if copy_mode == 'sample':
    data2_p['sample_matched'] = data2_p['sample']
elif copy_mode == 'table':
    data2_p = dat.copy_based_on_index(data2_p, table, 'sample', 'screen2', ['sample_matched'], ['screen1'])

data_p = data1_p.copy()
data_p = dat.copy_based_on_index(data_p, data2_p, 'sample', 'sample_matched',
                                 ['mob1', 'curve_mob1', 't_half1', 'curve_t_half1', 'slope1', 'curve_slope1',
                                  'size_organelle1', 'raw_int_organelle1', 'circ_organelle1', 'ecce_organelle1'],
                                 ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope',
                                  'size_organelle', 'raw_int_organelle', 'circ_organelle', 'ecce_organelle'])
data_p = data_p.dropna()

# export flow plots
for i in range(len(features)):
    dis.plot_comparison_flow(data_p, data_value, features[i], save_source)

print("DONE!")
