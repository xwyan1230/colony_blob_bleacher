import pandas as pd
import shared.display as dis
import numpy as np
import os

data_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210408_CBB_nucleoliFRAPscreen1/plate1/summary/WT/"

features = ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope', 'size_organelle', 'raw_int_organelle',
            'circ_organelle', 'ecce_organelle']
centers = [0.67, 0.68, 1.325, 1.17, 0.35, 0.475, 130, 2250, 0.8, 0.65]

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




