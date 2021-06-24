import pandas as pd
import shared.display as dis
import shared.dataframe as dat
import os

data1_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210610_CBB_nucleoliFRAPscreenComparison/summaryFiles/4/summary/NT/"
data2_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210610_CBB_nucleoliFRAPscreenComparison/summaryFiles/6-1/summary/NT/"
save_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210610_CBB_nucleoliFRAPscreenComparison/comparison/4_vs_6-1/"


copy_mode = 'table'  # accepts 'sample' or 'table'
features = ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope', 'size_organelle', 'raw_int_organelle',
            'circ_organelle', 'ecce_organelle']

if not os.path.exists(save_source):
    os.makedirs(save_source)

data1 = pd.read_csv(("%ssummary_value.txt" % data1_source), na_values=['.'], sep='\t')
data2 = pd.read_csv(("%ssummary_value.txt" % data2_source), na_values=['.'], sep='\t')

if copy_mode == 'sample':
    data2['sample_matched'] = data2['sample']
elif copy_mode == 'table':
    table = pd.read_csv(("%s2screenSampleTable.txt" % save_source), na_values=['.'], sep='\t')
    data2 = dat.copy_based_on_index(data2, table, 'sample', 'screen2', ['sample_matched'], ['screen1'])

data = data1.copy()
data = dat.copy_based_on_index(data, data2, 'sample', 'sample_matched',
                               ['mob1', 'curve_mob1', 't_half1', 'curve_t_half1', 'slope1', 'curve_slope1',
                                'size_organelle1', 'raw_int_organelle1', 'circ_organelle1', 'ecce_organelle1'],
                               ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope',
                                'size_organelle', 'raw_int_organelle', 'circ_organelle', 'ecce_organelle'])
data = data.dropna()

# export dot plots
for i in range(len(features)):
    dis.plot_comparison_dot(data, features[i], save_source)

print("DONE!")
