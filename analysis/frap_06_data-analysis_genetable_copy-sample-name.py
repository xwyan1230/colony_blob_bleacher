import pandas as pd
import os

data_source = "D:/Xiaowei/data/20210510_ctrlComparison/"
WT = 'WT'

data_p_source = ("%sdataSummary/%s/summary/" % (data_source, WT))
data_value_source = ("%sdataPlots/%s/summary/" % (data_source, WT))
data_stdev_source = ("%sdataStdev/%s/summary" % (data_source, WT))
data_p = pd.read_csv(("%ssummary_p.txt" % data_p_source), na_values=['.'], sep='\t')
data_value = pd.read_csv(("%ssummary_value.txt" % data_value_source), na_values=['.'], sep='\t')
data_stdev = pd.read_csv(("%ssummary_stdev.txt" % data_stdev_source), na_values=['.'], sep='\t')
save_path = ("%ssummary/%s/" % (data_source, WT))
if not os.path.exists(save_path):
    os.makedirs(save_path)

data_p['gene'] = data_p['sample']
data_value['gene'] = data_value['sample']
data_stdev['gene'] = data_stdev['sample']

data_p.to_csv('%ssummary_p.txt' % save_path, index=False, sep='\t')
data_value.to_csv('%ssummary_value.txt' % save_path, index=False, sep='\t')
data_stdev.to_csv('%ssummary_stdev.txt' % save_path, index=False, sep='\t')

print("DONE!")