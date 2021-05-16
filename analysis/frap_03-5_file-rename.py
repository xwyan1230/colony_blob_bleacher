import os
import pandas as pd

multi_data_source = "D:/Xiaowei/data/J_1_treatment/"

transform_name = 'J'

analyze_organelle = 'nucleoli'  # only accepts 'sg' or 'nucleoli'

transform_table = pd.DataFrame({'row': ['B', 'C', 'D', 'E', 'F', 'G'], 'number': [0, 10, 20, 30, 40, 50]})

multi_dirs = [x for x in os.listdir(multi_data_source)]
if '.DS_Store' in multi_dirs:
    multi_dirs.remove('.DS_Store')

for r in range(len(multi_dirs)):
    target_number = int(transform_table[transform_table['row'] == multi_dirs[r][0]]['number']+int(multi_dirs[r][1:])-1)
    target = ("%s%d" % (transform_name, target_number))
    data_source = ("%s%s/" % (multi_data_source, multi_dirs[r]))
    os.rename("%s%s_data_ctrl.txt" % (data_source, multi_dirs[r]), "%s%s_data_ctrl.txt" % (data_source, target))
    os.rename("%s%s_data_full.txt" % (data_source, multi_dirs[r]), "%s%s_data_full.txt" % (data_source, target))
    os.rename("%s%s_data_log.txt" % (data_source, multi_dirs[r]), "%s%s_data_log.txt" % (data_source, target))
    os.rename("%s%s_data_%s.txt" % (data_source, multi_dirs[r], analyze_organelle),
              "%s%s_data_%s.txt" % (data_source, target, analyze_organelle))
    if analyze_organelle == 'nucleoli':
        os.rename("%s%s_data_nuclear.txt" % (data_source, multi_dirs[r]), "%s%s_data_nuclear.txt" %
                  (data_source, target))
    os.rename(data_source, "%s%s/" % (multi_data_source, target))

print("DONE!")