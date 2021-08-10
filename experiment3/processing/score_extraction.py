import csv
import os
import sys 
import pandas as pd

################################################################################
#
# This script extract the relevant basecall quality scores from the reads which
# have reached this stage of analysis. Thus, we are able to generate the score 
# plots.
#
################################################################################

con1_first = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/firstData.csv")
print(con1_first.head())
del con1_first["current_stdev"]
del con1_first["current_avg"]
del con1_first["dwell"]
con1_first.drop(con1_first.columns[0],axis=1,inplace=True)
print(con1_first.head())
dicti = con1_first.set_index("id").T.to_dict("list")
print(len(list(dicti.items())))
print(list(dicti.items())[0:2])

path = "/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/2_tsv_conversion/construct1_for_map_sort.tsv"

cm = "wc -l"
cmd = " ".join([cm,path])

stream = os.popen(cmd)
output = stream.read()
lines = int(output.split(" ")[0])

csv.field_size_limit(sys.maxsize)

print("Opening in and out files.\n")
with open(path,"rt") as src:
	csrc = csv.reader(src, delimiter='\t')
	currentRead = ""
	first = pd.DataFrame({"id": [], "score": []})
	second = pd.DataFrame({"id": [], "score": []})
	third = pd.DataFrame({"id": [], "score": []})
	fourth = pd.DataFrame({"id": [], "score": []})
	fifth = pd.DataFrame({"id": [], "score": []})
	sixth = pd.DataFrame({"id": [], "score": []})
	seventh = pd.DataFrame({"id": [], "score": []})
	eight = pd.DataFrame({"id": [], "score": []})
	ninth = pd.DataFrame({"id": [], "score": []})
	inFlag = False
	inDictCount = 0
	for num, row in enumerate(csrc):
		if num == 0:
			continue
		if num%1000000==0:
			print(num,"of",lines)
		if row[0] != currentRead:
			currentRead = row[0]
			inFlag = False
			if currentRead in dicti:
				#print(currentRead,"is in dictionary")
				#print(dicti[currentRead])
				inDictCount+=1
				inFlag = True
			else:
				continue
		elif inFlag == True:
			try:
				int(row[7])
			except ValueError:
				continue
			if int(row[7]) >= 166 and int(row[7]) <= 174:
				if int(row[7]) == 166:
					first = first.append({"id":currentRead,"score":row[6]},ignore_index=True)
				elif int(row[7]) == 167:
					second = second.append({"id":currentRead,"score":row[6]},ignore_index=True)
				elif int(row[7]) == 168:
					third = third.append({"id":currentRead,"score":row[6]},ignore_index=True)
				elif int(row[7]) == 169:
					fourth = fourth.append({"id":currentRead,"score":row[6]},ignore_index=True)
				elif int(row[7]) == 170:
					fifth = fifth.append({"id":currentRead,"score":row[6]},ignore_index=True)
				elif int(row[7]) == 171:
					sixth = sixth.append({"id":currentRead,"score":row[6]},ignore_index=True)
				elif int(row[7]) == 172:
					seventh = seventh.append({"id":currentRead,"score":row[6]},ignore_index=True)
				elif int(row[7]) == 173:
					eight = eight.append({"id":currentRead,"score":row[6]},ignore_index=True)
				elif int(row[7]) == 174:
					ninth = ninth.append({"id":currentRead,"score":row[6]},ignore_index=True)
				dicti.pop(currentRead,None)
		else:
			continue
		#if num == 2:
			#exit()

first.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/test_firstScore/firstData.csv")
second.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/test_firstScore/secondData.csv")
third.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/test_firstScore/thirdData.csv")
fourth.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/test_firstScore/fourthData.csv")
fifth.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/test_firstScore/fifthData.csv")
sixth.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/test_firstScore/sixthData.csv")
seventh.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/test_firstScore/seventhData.csv")
eight.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/test_firstScore/eightData.csv")
ninth.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/test_firstScore/ninthData.csv")

print("\nCompleted with",inDictCount,"reads detected in dictionary")
