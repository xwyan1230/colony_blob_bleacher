import pandas as pd
import os
from datetime import datetime

master_folder = "D:/Xiaowei/data/20210607_screen/"
data_source = "D:/Xiaowei/data/20210607_screen/dataFiles/"
save_source = "D:/Xiaowei/data/20210607_screen/WTFiles/NT/"
WT_lst = ['E12', 'F12', 'G12']
analyze_organelle = 'nucleoli'  # only accepts 'nucleoli' or 'sg'

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
f.write("save_source: %s\n" % save_source)
f.write("WT_list: %s\n" % WT_lst)
f.write("analyze_organelle: %s\n\n" % analyze_organelle)

f.close()

# script starts
data_WT = pd.DataFrame()
data_ctrl_WT = pd.DataFrame()
data_log_WT = pd.DataFrame()
data_nuclear_WT = pd.DataFrame()
data_nucleoli_WT = pd.DataFrame()
data_sg_WT = pd.DataFrame()

for i in range(len(WT_lst)):
    data_WT_temp = pd.read_csv(("%s%s/%s_data_full.txt" % (data_source, WT_lst[i], WT_lst[i])),
                               na_values=['.'], sep='\t')
    data_ctrl_WT_temp = pd.read_csv(("%s%s/%s_data_ctrl.txt" % (data_source, WT_lst[i], WT_lst[i])),
                                    na_values=['.'], sep='\t')
    data_log_WT_temp = pd.read_csv(("%s%s/%s_data_log.txt" % (data_source, WT_lst[i], WT_lst[i])),
                                   na_values=['.'], sep='\t')
    data_WT = pd.concat([data_WT, data_WT_temp])
    data_ctrl_WT = pd.concat([data_ctrl_WT, data_ctrl_WT_temp])
    data_log_WT = pd.concat([data_log_WT, data_log_WT_temp])

    if analyze_organelle == 'nucleoli':
        data_nuclear_WT_temp = pd.read_csv(("%s%s/%s_data_nuclear.txt" % (data_source, WT_lst[i], WT_lst[i])),
                                           na_values=['.'], sep='\t')
        data_nucleoli_WT_temp = pd.read_csv(("%s%s/%s_data_nucleoli.txt" % (data_source, WT_lst[i], WT_lst[i])),
                                            na_values=['.'], sep='\t')
        data_nuclear_WT = pd.concat([data_nuclear_WT, data_nuclear_WT_temp])
        data_nucleoli_WT = pd.concat([data_nucleoli_WT, data_nucleoli_WT_temp])
    elif analyze_organelle == 'sg':
        data_sg_WT_temp = pd.read_csv(("%s%s/%s_data_sg.txt" % (data_source, WT_lst[i], WT_lst[i])),
                                      na_values=['.'], sep='\t')
        data_sg_WT = pd.concat([data_sg_WT, data_sg_WT_temp])

if not os.path.exists(save_source):
    os.makedirs(save_source)
data_WT.to_csv('%sWT_data_full.txt' % save_source, index=False, sep='\t')
data_ctrl_WT.to_csv('%sWT_data_ctrl.txt' % save_source, index=False, sep='\t')
data_log_WT.to_csv('%sWT_data_log.txt' % save_source, index=False, sep='\t')
if analyze_organelle == 'nucleoli':
    data_nuclear_WT.to_csv('%sWT_data_nuclear.txt' % save_source, index=False, sep='\t')
    data_nucleoli_WT.to_csv('%sWT_data_nucleoli.txt' % save_source, index=False, sep='\t')
elif analyze_organelle == 'sg':
    data_sg_WT.to_csv('%sWT_data_sg.txt' % save_source, index=False, sep='\t')

print("DONE!")