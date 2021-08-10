import csv
import sys
import pandas as pd

################################################################################
#
# This script normalizes the data used for classification. It takes 10% of the
# data, averages it, and uses that to normalize the remaining 90%.
#
################################################################################

def reader(con, seg, pos):
	path = "/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/5_classification/score_merge/construct"+con+"/"+seg+"/"+pos+"PosScoreMerge.csv"
	return pd.read_csv(path)

def writer(con, seg, pos):
	df.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/5_classification/normalized/construct"+con+"/"+seg+"/"+pos+"PosNorm.csv")

constructs = ["2","3","4","5"]
segments = ["first","second","third"]
positions = ["first","second","third","fourth","fifth","sixth","seventh","eight","ninth"]

for construct in constructs:
	for segment in segments:
		for position in positions:

			con, seg, pos = construct,segment,position

			df = reader(con,seg,pos)
			print(df.head())
			print("sampling:",round(len(df.index)/10))
			sample = df.sample(frac=0.1)
			print(sample.head())
			df = df.drop(sample.index)

			normalization_current = sample["current_avg"].mean(axis=0)
			normalization_dwell = sample["dwell"].mean(axis=0)
			print(normalization_current, normalization_dwell)
			df["current_avg"] = df["current_avg"]-normalization_current
			df["dwell"] = df["dwell"]-normalization_dwell
			df = df.drop(df.columns[0], axis=1)
			df = df.drop(df.columns[0], axis=1)
			print(df.head())


			writer(con,seg,pos)

