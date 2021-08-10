import csv
import pandas as pd
import os

################################################################################
#
# This script simply prepares the data gathered by the fast5_processing.py 
# script for further analysis. The formatting here is conducive to further 
# position-by-position analysis of each of the 9-base segments. 
#
################################################################################

path = "/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/construct4/perfs/third_pos_data_out.csv"

cm = "wc -l"
cmd = " ".join([cm,path])

stream = os.popen(cmd)
output = stream.read()
lines = int(output.split(" ")[0])

file = open(path)
read_csv = csv.reader(file, delimiter=",")

curr_read = ""
perfect = []
#first = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
#['0', '001e5f4f-afcb-4c20-a96d-3c11325b0d68', 'G', '56.895743991794255', '3.2241089110144263', '30.0']
first = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
second = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
third = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
fourth = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
fifth = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
sixth = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
seventh = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
eight = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
ninth = pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]})
read = [] #pd.DataFrame({"id": [], "base": [], "current_avg": [], "current_stdev":[], "dwell":[]}) 
print(lines,"to read")
for index, row in enumerate(read_csv):
	if index == 0:
		continue
	#if index == 29:
	#	print(first)
	#	print(fifth)
	#	print(ninth)
	#	exit()
	#if index == 50:
		#print(len(perfect))
		#print(first)
		#print(ninth)
		#print(fifth)
	if index%1000==0:
		print(index,"out of",lines)
	if curr_read != row[1]: #new read
		if read != [] and count == 9:
			#print("here with", curr_read)
			perfect.append(read)
			newFirst = {"id":read[0][1],"base":read[0][2],"current_avg":read[0][3],"current_stdev":read[0][4],"dwell":read[0][5]}
			first = first.append(newFirst,ignore_index=True)
			second = second.append({"id":read[1][1],"base":read[1][2],"current_avg":read[1][3],"current_stdev":read[1][4],"dwell":read[1][5]},ignore_index=True)
			third = third.append({"id":read[2][1],"base":read[2][2],"current_avg":read[2][3],"current_stdev":read[2][4],"dwell":read[2][5]},ignore_index=True)
			fourth = fourth.append({"id":read[3][1],"base":read[3][2],"current_avg":read[3][3],"current_stdev":read[3][4],"dwell":read[3][5]},ignore_index=True)
			fifth = fifth.append({"id":read[4][1],"base":read[4][2],"current_avg":read[4][3],"current_stdev":read[4][4],"dwell":read[4][5]},ignore_index=True)
			sixth = sixth.append({"id":read[5][1],"base":read[5][2],"current_avg":read[5][3],"current_stdev":read[5][4],"dwell":read[5][5]},ignore_index=True)
			seventh = seventh.append({"id":read[6][1],"base":read[6][2],"current_avg":read[6][3],"current_stdev":read[6][4],"dwell":read[6][5]},ignore_index=True)
			eight = eight.append({"id":read[7][1],"base":read[7][2],"current_avg":read[7][3],"current_stdev":read[7][4],"dwell":read[7][5]},ignore_index=True)
			ninth = ninth.append({"id":read[8][1],"base":read[8][2],"current_avg":read[8][3],"current_stdev":read[8][4],"dwell":read[8][5]},ignore_index=True)
		count = 1
		#print("here with", row[1],"count is",count)
		curr_read = row[1]
		read = [row]
	elif curr_read == row[1]:
		count+=1
		#print("here with", row[1],"count is",count)
		read.append(row)
	#if count == 9:
               	#print("Successful read")

print("\nwriting all to file.\n")
first.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/thirdPos/firstData.csv")
second.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/thirdPos/secondData.csv")
third.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/thirdPos/thirdData.csv")
fourth.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/thirdPos/fourthData.csv")
fifth.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/thirdPos/fifthData.csv")
sixth.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/thirdPos/sixthData.csv")
seventh.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/thirdPos/seventhData.csv")
eight.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/thirdPos/eightData.csv")
ninth.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/thirdPos/ninthData.csv")
print("\ndone.\n")
