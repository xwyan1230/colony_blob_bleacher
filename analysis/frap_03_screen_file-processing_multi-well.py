import os
import shared.dataframe as dat
from datetime import datetime

master_folder = "D:/Xiaowei/data/20210607_screen/"
multi_data_source = "D:/Xiaowei/data/20210607_screen/dataAnalysis/"
save_source = "D:/Xiaowei/data/20210607_screen/dataFiles/"

analyze_organelle = 'nucleoli'  # only accepts 'sg' or 'nucleoli'

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
f.write("multi_data_source: %s\n" % multi_data_source)
f.write("save_source: %s\n" % save_source)
f.write("analyze_organelle: %s\n\n" % analyze_organelle)

f.close()

# script starts
multi_dirs = [x for x in os.listdir(multi_data_source)]
if '.DS_Store' in multi_dirs:
    multi_dirs.remove('.DS_Store')

for r in range(len(multi_dirs)):
    data_source = ("%s%s/" % (multi_data_source, multi_dirs[r]))
    dirs = [x[0] for x in os.walk(data_source)]
    dirs.pop(0)
    num_dir = len(dirs)
    name = multi_dirs[r]

    save_path = ("%s%s/" % (save_source, name))
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    f_full = open("%s/%s_data_full.txt" % (save_path, name), 'w+')
    f_log = open("%s/%s_data_log.txt" % (save_path, name), 'w+')
    f_ctrl = open("%s/%s_data_ctrl.txt" % (save_path, name), 'w+')
    if analyze_organelle == 'nucleoli':
        f_nuclear = open("%s/%s_data_nuclear.txt" % (save_path, name), 'w+')
        f_nucleoli = open("%s/%s_data_nucleoli.txt" % (save_path, name), 'w+')
    elif analyze_organelle == 'sg':
        f_sg = open("%s/%s_data_sg.txt" % (save_path, name), 'w+')

    count = 0
    for s in range(len(dirs)):
        data_path = dirs[s]
        if os.path.exists("%s/data_full.txt" % data_path):
            f1_full = open("%s/data_full.txt" % data_path, 'r')
            dat.append_data(f_full, f1_full, count)
            f1_full.close()
            f1_log = open("%s/data_log.txt" % data_path, 'r')
            dat.append_data(f_log, f1_log, count)
            f1_log.close()
            f1_ctrl = open("%s/data_ctrl.txt" % data_path, 'r')
            dat.append_data(f_ctrl, f1_ctrl, count)
            f1_ctrl.close()
            if analyze_organelle == 'nucleoli':
                f1_nuclear = open("%s/data_nuclear.txt" % data_path, 'r')
                dat.append_data(f_nuclear, f1_nuclear, count)
                f1_nuclear.close()
                f1_nucleoli = open("%s/data_nucleoli.txt" % data_path, 'r')
                dat.append_data(f_nucleoli, f1_nucleoli, count)
                f1_nucleoli.close()
            elif analyze_organelle == 'sg':
                f1_sg = open("%s/data_sg.txt" % data_path, 'r')
                dat.append_data(f_sg, f1_sg, count)
                f1_sg.close()
            count = count + 1

    f_full.close()
    f_log.close()
    f_ctrl.close()
    if analyze_organelle == 'nucleoli':
        f_nuclear.close()
        f_nucleoli.close()
    elif analyze_organelle == 'sg':
        f_sg.close()

print('Done')
