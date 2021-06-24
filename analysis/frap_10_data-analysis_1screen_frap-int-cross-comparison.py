import pandas as pd
from matplotlib import pyplot as plt
import os

data_source = "D:/Xiaowei/data/20210607_screen/"
WT = 'NT'
save_source = "D:/Xiaowei/data/20210607_screen/crossComparison/"

features = ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope']

if not os.path.exists(save_source):
    os.makedirs(save_source)

data_value_source = ("%sdataPlots/%s/summary/" % (data_source, WT))
data_stdev_source = ("%sdataStdev/%s/summary/" % (data_source, WT))
data_value = pd.read_csv(("%ssummary_value.txt" % data_value_source), na_values=['.'], sep='\t')
data_stdev = pd.read_csv(("%ssummary_stdev.txt" % data_stdev_source), na_values=['.'], sep='\t')

for i in features:
    plt.figure(figsize=(9, 6))
    plt.errorbar(data_value['raw_int_organelle'], data_value[i], data_stdev[i], linestyle='None', marker='^')
    plt.xlabel('raw_int')
    plt.ylabel(i)

    plt.savefig('%s%s_comparedWithRawInt.pdf' % (save_source, i))
    plt.close()

print("DONE!")