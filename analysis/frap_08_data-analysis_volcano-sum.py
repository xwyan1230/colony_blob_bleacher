import pandas as pd
import shared.display as dis
import os

data_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210415_CBB_nucleoliFRAPscreen2/plate2/summary/NT/"
save_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210415_CBB_nucleoliFRAPscreen2/plate2/volcanoPlots/NT/"

data_p = pd.read_csv(("%ssummary_p.txt" % data_source), na_values=['.'], sep='\t')
data_value = pd.read_csv(("%ssummary_value.txt" % data_source), na_values=['.'], sep='\t')

if not os.path.exists(save_source):
    os.makedirs(save_source)

# export volcano plots
dis.plot_volcano(data_p, data_value, 'mob', save_source)
dis.plot_volcano(data_p, data_value, 'curve_mob', save_source)
dis.plot_volcano(data_p, data_value, 't_half', save_source)
dis.plot_volcano(data_p, data_value, 'curve_t_half', save_source)
dis.plot_volcano(data_p, data_value, 'slope', save_source)
dis.plot_volcano(data_p, data_value, 'curve_slope', save_source)
dis.plot_volcano(data_p, data_value, 'size_organelle', save_source)
dis.plot_volcano(data_p, data_value, 'raw_int_organelle', save_source)
dis.plot_volcano(data_p, data_value, 'circ_organelle', save_source)
dis.plot_volcano(data_p, data_value, 'ecce_organelle', save_source)

print("DONE!")
