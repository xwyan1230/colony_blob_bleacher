import pandas as pd
import os
import shared.dataframe as dat

data_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/20210415_CBB_nucleoliFRAPscreen2/"\
    "plate2/"
WT = 'WT'

data_p_source = ("%sdataSummary/%s/summary/" % (data_source, WT))
data_value_source = ("%sdataPlots/%s/summary/" % (data_source, WT))
data_p = pd.read_csv(("%ssummary_p.txt" % data_p_source), na_values=['.'], sep='\t')
data_value = pd.read_csv(("%ssummary_value.txt" % data_value_source), na_values=['.'], sep='\t')
genetable = pd.read_csv(("%sgenetable.txt" % data_source), na_values=['.'], sep='\t')
save_path = ("%ssummary/%s/" % (data_source, WT))
if not os.path.exists(save_path):
    os.makedirs(save_path)

data_p = dat.correlate_genetable(data_p, genetable)
data_value = dat.correlate_genetable(data_value, genetable)

data_p.to_csv('%ssummary_p.txt' % save_path, index=False, sep='\t')
data_value.to_csv('%ssummary_value.txt' % save_path, index=False, sep='\t')

print("DONE!")