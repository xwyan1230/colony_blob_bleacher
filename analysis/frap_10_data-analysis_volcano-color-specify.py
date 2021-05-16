import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os

data_source = "D:/Xiaowei/data/20210510_ctrlComparison/summary/WT/"
save_source = "D:/Xiaowei/data/20210510_ctrlComparison/volcanoPlots/vs_threeONrun/"
WT_source = "D:/Xiaowei/data/20210510/WTFiles/WT_B2-B3-B4/"

analyze_organelle = 'nucleoli'
analysis_mode = 'single_exp'
features = ['mob', 'curve_mob', 't_half', 'curve_t_half', 'slope', 'curve_slope', 'size_organelle', 'raw_int_organelle',
            'circ_organelle', 'ecce_organelle']

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
data_gamma = pd.read_csv(("%ssummary_gamma.txt" % data_source), na_values=['.'], sep='\t')

if not os.path.exists(save_source):
    os.makedirs(save_source)


def volcano(pd_p: pd.DataFrame, pd_value: pd.DataFrame, pd_gamma: pd.DataFrame, feature: str, center: float,
            threshold: float, show_gene: str, save_path: str):
    """
    Plot volcano plot for FRAP screen (with hit gene labeled)

    :param pd_p: pd.DataFrame, -ln(p) value dataframe
    :param pd_value: pd.DataFrame, average value dataframe
    :param pd_gamma: pd.DataFrame, score gamma (value-center)*[-ln(p)] dataframe
    :param feature: str, plotting feature, column name
    :param center: float, center used for calculating gamma, get from volcano plot-sum
    :param threshold: float, threshold used to call hit
                      fold of value difference from center/ std(control values) * 3 (roughly -ln(0.05))
    :param show_gene: str, only accepts 'N' or 'Y', whether display gene name on volcano plot
    :param save_path: str, saving path
    :return:
    """
    plt.figure(figsize=(9, 6))

    ctrl_std = np.std(pd_value[pd_value['gene'].str.contains('P')][feature])
    ishit = pd_gamma[feature] / ctrl_std >= threshold

    plt.scatter(pd_value[pd_value['gene'].str.contains('P')][feature], pd_p[pd_p['gene'].str.contains('P')][feature],
                s=4, alpha=0.8, c='#808080', label='20210510')
    #plt.scatter(pd_value[pd_value['gene'].str.contains('I')][feature], pd_p[pd_p['gene'].str.contains('I')][feature],
    #            s=4, alpha=0.8, c='#FFA500', label='Trial1 Ctrl')
    #plt.scatter(pd_value[pd_value['gene'].str.contains('J')][feature], pd_p[pd_p['gene'].str.contains('J')][feature],
    #            s=4, alpha=0.8, c='#DC143C', label='Trial1 sample')
    #plt.scatter(pd_value[pd_value['gene'].str.contains('K|Q')][feature],
    #            pd_p[pd_p['gene'].str.contains('K|Q')][feature], s=4, alpha=0.8, c='#9ACD32', label='Trial2 Ctrl')
    #plt.scatter(pd_value[pd_value['gene'].str.contains('L|R')][feature],
    #            pd_p[pd_p['gene'].str.contains('L|R')][feature], s=4, alpha=0.8, c='#008000', label='Trial2 sample')
    plt.scatter(pd_value[pd_value['gene'].str.contains('M')][feature], pd_p[pd_p['gene'].str.contains('M')][feature],
                s=4, alpha=0.8, c='#87CEEB', label='20210503')
    plt.scatter(pd_value[pd_value['gene'].str.contains('N')][feature], pd_p[pd_p['gene'].str.contains('N')][feature],
                s=4, alpha=0.8, c='#4169E1', label='20210504')
    plt.scatter(pd_value[pd_value['gene'].str.contains('O')][feature], pd_p[pd_p['gene'].str.contains('O')][feature],
                s=4, alpha=0.8, c='#00008B', label='20210505')


    ymax = np.ceil(max(pd_p[feature])) * 1.02
    xmin = min(pd_value[feature]) * 0.95
    xmax = max(pd_value[feature]) * 1.1

    plt.plot(np.linspace(xmin, xmax, 1000),
             np.abs(threshold / np.linspace((xmin-center) / ctrl_std, (xmax-center) / ctrl_std, 1000)), 'k--', lw=.5)

    if show_gene == 'Y':
        x_lst = pd_value[ishit][feature].tolist()
        y_lst = pd_p[ishit][feature].tolist()
        gene_lst = pd_value[ishit]['gene'].tolist()
        for j in range(len(x_lst)):
            plt.text(x_lst[j], y_lst[j], gene_lst[j], fontsize=6)

    plt.xlim((xmin, xmax))
    plt.ylim((0, ymax))
    plt.xlabel(feature)
    plt.ylabel('-ln(p)')
    plt.legend(loc=4, bbox_to_anchor=(0.7, 0, 0.3, 0.3))

    plt.savefig('%s%s_hit.pdf' % (save_path, feature))
    plt.close()


# export volcano hit plots
for i in range(len(features)):
    volcano(data_p, data_value, data_gamma, features[i], centers[i], 60, 'N', save_source)

print("DONE!")
