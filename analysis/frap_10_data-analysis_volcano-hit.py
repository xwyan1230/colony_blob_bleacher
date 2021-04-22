import pandas as pd
import shared.display as dis
import os

data_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210408_CBB_nucleoliFRAPscreen1/plate1/summary/WT/"
save_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210408_CBB_nucleoliFRAPscreen1/plate1/volcanoPlots/WT/"

features = ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope', 'size_organelle', 'raw_int_organelle',
            'circ_organelle', 'ecce_organelle']
centers = [0.67, 0.68, 1.325, 1.17, 0.35, 0.475, 130, 2250, 0.8, 0.65]

data_p = pd.read_csv(("%ssummary_p.txt" % data_source), na_values=['.'], sep='\t')
data_value = pd.read_csv(("%ssummary_value.txt" % data_source), na_values=['.'], sep='\t')
data_gamma = pd.read_csv(("%ssummary_gamma.txt" % data_source), na_values=['.'], sep='\t')

if not os.path.exists(save_source):
    os.makedirs(save_source)

# export volcano hit plots
for i in range(len(features)):
    dis.plot_volcano_hit(data_p, data_value, data_gamma, features[i], centers[i], 60, 'Y', save_source)

print("DONE!")
