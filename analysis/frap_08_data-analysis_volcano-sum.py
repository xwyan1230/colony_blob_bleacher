import pandas as pd
import shared.display as dis
import os

data_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210408_CBB_nucleoliFRAPscreen1/plate1/summary/WT/"
save_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210408_CBB_nucleoliFRAPscreen1/plate1/volcanoPlots/WT/"

features = ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope', 'size_organelle', 'raw_int_organelle',
            'circ_organelle', 'ecce_organelle']

data_p = pd.read_csv(("%ssummary_p.txt" % data_source), na_values=['.'], sep='\t')
data_value = pd.read_csv(("%ssummary_value.txt" % data_source), na_values=['.'], sep='\t')
data_gamma = pd.read_csv(("%ssummary_gamma.txt" % data_source), na_values=['.'], sep='\t')

if not os.path.exists(save_source):
    os.makedirs(save_source)

# export volcano plots
for i in range(len(features)):
    dis.plot_volcano(data_p, data_value, features[i], save_source)

print("DONE!")
