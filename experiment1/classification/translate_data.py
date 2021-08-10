import pandas as pd
import numpy as np

################################################################################
#
# This script takes the output generated from the 
# ../raw_current_analysis/third_file_parse.R script and translates it into 
# a format thta can be used for classification in python.
#
################################################################################

data = pd.read_csv("expanded_current_row_variable_dataset.csv")
current_row = 9
current_base = 'I'
current_id = '06ec157a-337d-4114-baff-4442bf860a97'
current_dwell = 122
row_values = []
row_count = 0
base_count = 0
#new_df_w_mean = pd.DataFrame(columns=["mean","row_nr","base"]) 
#new_df_w_all_signals = pd.DataFrame(columns=["values","row_nr","base"])
new_df_mean_sd_dwell = pd.DataFrame(columns=["mean","sd","dwell","read_id","row_nr","base"])
total = 13535

print("Starting data processing")

for idx,row in data.iterrows():
	if base_count > 2707 or current_base != row[4]:
		if current_base == row[4]:
			continue
		else:
			print("new base", row[4])
			current_base = row[4]
			#current_row = int(row[1])
			current_id = row[3]
			current_dwell = int(row[2])
			row_values = [float(row[0])]
			base_count = 1
	if current_row == int(row[1]):
		if(current_base == 'T'):
			print("on the same row")
		row_values.append(float(row[0]))
	else:
		if current_base == 'T':
			print("changing row from",current_row,"to",row[1])
		if current_base == row[4]:
			base_count+=1
		else:
			base_count==1
		#new_df_w_mean.loc[len(new_df_w_mean)] = [np.mean(row_values)]+[current_row]+[current_base]
		#new_df_w_all_signals.loc[len(new_df_w_all_signals)] = [row_values]+[current_row]+[current_base]
		new_df_mean_sd_dwell.loc[len(new_df_mean_sd_dwell)] =[np.mean(row_values)]+[np.std(row_values)]+[current_dwell]+[current_id]+[current_row]+[current_base]
		row_values = [float(row[0])]
		current_row = int(row[1])
		current_id = row[3]
		current_base = row[4]
		current_dwell = int(row[2])
		row_count+=1
		if row_count%100 == 0:
			print("Done",str(row_count)+"/"+str(total))

print("All done, writing to file.")
#new_df_w_mean.to_csv('balanced_dataframe_w_mean.csv')
#new_df_w_all_signals.to_csv('balanced_dataframe_w_row.csv')
new_df_mean_sd_dwell.to_csv('balanced_dataframe_w_mean_sd_dwell.csv')
