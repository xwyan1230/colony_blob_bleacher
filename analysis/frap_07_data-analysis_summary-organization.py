import pandas as pd
import os
import shared.dataframe as dat

multi_data_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/20210408_CBB_nucleoliFRAPscreen1/"
save_source = "/Users/xiaoweiyan/Dropbox/LAB/ValeLab/Projects/Blob_bleacher/Exp/"\
    "20210408_CBB_nucleoliFRAPscreen1_summary/summary/"
WT = 'WT'

multi_dirs = [x for x in os.listdir(multi_data_source)]
if '.DS_Store' in multi_dirs:
    multi_dirs.remove('.DS_Store')

data_p = pd.DataFrame()
data_value = pd.DataFrame()

for r in range(len(multi_dirs)):
    name = multi_dirs[r]
    number = name[5:]
    print('# Organizing %s ... (%d/%d)' % (name, r+1, len(multi_dirs)))
    data_p_source = ("%s%s/dataSummary/%s/summary/" % (multi_data_source, name, WT))
    data_value_source = ("%s%s/dataPlots/%s/summary/" % (multi_data_source, name, WT))
    genetable_source = ("%s%s/" % (multi_data_source, name))
    save_path = ("%s%s/" % (save_source, WT))
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    data_p_temp = pd.read_csv(("%ssummary_p.txt" % data_p_source), na_values=['.'], sep='\t')
    data_p_temp['plate'] = [number] * len(data_p_temp)
    data_value_temp = pd.read_csv(("%ssummary_value.txt" % data_value_source), na_values=['.'], sep='\t')
    data_value_temp['plate'] = [number] * len(data_value_temp)
    genetable_temp = pd.read_csv(("%sgenetable.txt" % genetable_source), na_values=['.'], sep='\t')
    data_p_temp = dat.correlate_genetable(data_p_temp, genetable_temp)
    data_value_temp = dat.correlate_genetable(data_value_temp, genetable_temp)

    data_p = pd.concat([data_p, data_p_temp], ignore_index=True)
    data_value = pd.concat([data_value, data_value_temp], ignore_index=True)

data_p.to_csv('%ssummary_p.txt' % save_path, index=False, sep='\t')
data_value.to_csv('%ssummary_value.txt' % save_path, index=False, sep='\t')

print("DONE!")