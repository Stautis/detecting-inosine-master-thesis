import csv
import sys
import pandas as pd

################################################################################
#
# This script is written to align the basecall score quality with the other
# features extracted from the data prior to plotting the comparisons. 
#
################################################################################

def reader(con, seg, fold, pos):
        path = "/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct"+con+"/"+seg+fold+"/"+pos+"Data.csv"
        return pd.read_csv(path)

list_of = ["first","second","third","fourth","fifth","sixth","seventh","eight","ninth"]

for number in list_of:
	
	print("on position",number)	

	df_con1_pos1 = reader("5","third","Pos",number)
	df_con1_score1 = reader("5","third","Score",number)

	print(df_con1_pos1.head())
	print(df_con1_score1.head())
	print((df_con1_score1["score"].map(ord)-33).head())

	sorted_df_con1_score1 = df_con1_score1.sort_values("id",inplace=False)
	#sorted_df_con1_pos1 = df_con1_pos1.sort_values("id",inplace=False)
	print(sorted_df_con1_score1.head())
	#print(sorted_df_con1_pos1.head())


	score_df = pd.DataFrame(df_con1_pos1[["id","dwell"]])
	score_df = score_df.rename(columns={"id":"id","dwell":"score"})
	print(score_df.head())


	for idx, i in df_con1_score1.iterrows():
		if idx%1000==0:#(len(df_con1_score1)/10)==0:
			print(idx ,"of",len(df_con1_score1))
		if i["id"] in list(df_con1_pos1["id"]):
			#print("Present with score",i["score"])
			#print(list(df_con1_pos1["id"]).index(i["id"]))
			#print(df_con1_pos1.loc[list(df_con1_pos1["id"]).index(i["id"])])
			#print(score_df.loc[list(df_con1_pos1["id"]).index(i["id"])])
			score_df.at[list(df_con1_pos1["id"]).index(i["id"]),'score']=(ord(i["score"])-33) #score_df["score"].loc[list(df_con1_pos1["id"]).index(i["id"])]=i["score"]
			#print(score_df.loc[list(df_con1_pos1["id"]).index(i["id"])])
		else:
			print(i ,"is Not")
			exit()

	sorted_score_df = score_df.sort_values("id",inplace=False)
	print(sorted_score_df.head())	
	df_con1_pos1["score"] = score_df["score"]
	df_con1_pos1.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis//5_classification/score_merge/construct5/third/"+number+"PosScoreMerge.csv")
	#print(score_df.head())
print("done.")
