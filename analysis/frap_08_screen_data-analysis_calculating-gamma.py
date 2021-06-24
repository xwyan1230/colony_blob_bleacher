import pandas as pd
import numpy as np
from datetime import datetime
import os

master_folder = "D:/Xiaowei/data/20210607_screen/"
data_source = "D:/Xiaowei/data/20210607_screen/summary/NT/"
WT_source = "D:/Xiaowei/data/20210607_screen/WTFiles/NT/"

analyze_organelle = 'nucleoli'
analysis_mode = 'single_exp'
features = ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope', 'size_organelle', 'raw_int_organelle',
            'circ_organelle', 'ecce_organelle']

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
f.write("data_source: %s\n" % data_source)
f.write("WT_source: %s\n" % WT_source)
f.write("analyze_organelle: %s\n" % analyze_organelle)
f.write("analysis_mode: %s\n" % analysis_mode)
f.write("features: %s\n\n" % features)

f.close()


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


def get_phenotypic_score(df_p: pd.DataFrame, df_val: pd.DataFrame, features: list, centers: list):
    df = pd.DataFrame()
    df['sample'] = df_p['sample']
    df['n_curve'] = df_p['n_curve']
    df['limit'] = df_p['limit']
    df['n_organelle'] = df_p['n_organelle']
    df['limit_organelle'] = df_p['limit_organelle']
    df['limit_organelle_circ'] = df_p['limit_organelle_circ']
    df['gene'] = df_p['gene']
    for i in range(len(features)):
        df[features[i]] = np.abs(df_val[features[i]] - centers[i]) * df_p[features[i]]
    return df


data_gamma = get_phenotypic_score(data_p, data_value, features, centers)
data_gamma.to_csv('%ssummary_gamma.txt' % data_source, index=False, sep='\t')

print("DONE!")




