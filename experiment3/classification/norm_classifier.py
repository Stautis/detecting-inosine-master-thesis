import csv
import sys
import pandas as pd
import h2o
from h2o.estimators import H2ORandomForestEstimator
h2o.init()

################################################################################
#
# This program written for classification of segments based on data present at
# all 9 positions in the segment. The response label can be specified to change
# from a 2-way to a 5-way classification task. Here we train on the first 
# segment, then we evaluate on the second segment. The response label here is 
# the construct label, and so we are doing a 5-way classification. This data has
# been normalized.
#
################################################################################

def reader(con, seg, pos):
        path = "/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/5_classification/normalized/construct"+con+"/"+seg+"/"+pos+"PosNorm.csv"
        return pd.read_csv(path)

file_list = ["first","second","third","fourth","fifth","sixth","seventh","eight","ninth"]

df_con1_pos1 = reader("1","first","first")
df_con1_pos2 = reader("1","first","second")
df_con1_pos3 = reader("1","first","third")
df_con1_pos4 = reader("1","first","fourth")
df_con1_pos5 = reader("1","first","fifth")
df_con1_pos6 = reader("1","first","sixth")
df_con1_pos7 = reader("1","first","seventh")
df_con1_pos8 = reader("1","first","eight")
df_con1_pos9 = reader("1","first","ninth")

first_cons = [df_con1_pos1,df_con1_pos2,df_con1_pos3,df_con1_pos4,df_con1_pos5,df_con1_pos6,df_con1_pos7,df_con1_pos8,df_con1_pos9]

#print(first_cons[:10])

df_con2_pos1 = reader("2","first","first")
df_con2_pos2 = reader("2","first","second")
df_con2_pos3 = reader("2","first","third")
df_con2_pos4 = reader("2","first","fourth")
df_con2_pos5 = reader("2","first","fifth")
df_con2_pos6 = reader("2","first","sixth")
df_con2_pos7 = reader("2","first","seventh")
df_con2_pos8 = reader("2","first","eight")
df_con2_pos9 = reader("2","first","ninth")

second_cons = [df_con2_pos1,df_con2_pos2,df_con2_pos3,df_con2_pos4,df_con2_pos5,df_con2_pos6,df_con2_pos7,df_con2_pos8,df_con2_pos9]

df_con3_pos1 = reader("3","first","first")
df_con3_pos2 = reader("3","first","second")
df_con3_pos3 = reader("3","first","third")
df_con3_pos4 = reader("3","first","fourth")
df_con3_pos5 = reader("3","first","fifth")
df_con3_pos6 = reader("3","first","sixth")
df_con3_pos7 = reader("3","first","seventh")
df_con3_pos8 = reader("3","first","eight")
df_con3_pos9 = reader("3","first","ninth")

third_cons = [df_con3_pos1,df_con3_pos2,df_con3_pos3,df_con3_pos4,df_con3_pos5,df_con3_pos6,df_con3_pos7,df_con3_pos8,df_con3_pos9]

df_con4_pos1 = reader("4","first","first")
df_con4_pos2 = reader("4","first","second")
df_con4_pos3 = reader("4","first","third")
df_con4_pos4 = reader("4","first","fourth")
df_con4_pos5 = reader("4","first","fifth")
df_con4_pos6 = reader("4","first","sixth")
df_con4_pos7 = reader("4","first","seventh")
df_con4_pos8 = reader("4","first","eight")
df_con4_pos9 = reader("4","first","ninth")

fourth_cons = [df_con4_pos1,df_con4_pos2,df_con4_pos3,df_con4_pos4,df_con4_pos5,df_con4_pos6,df_con4_pos7,df_con4_pos8,df_con4_pos9]

df_con5_pos1 = reader("5","first","first")
df_con5_pos2 = reader("5","first","second")
df_con5_pos3 = reader("5","first","third")
df_con5_pos4 = reader("5","first","fourth")
df_con5_pos5 = reader("5","first","fifth")
df_con5_pos6 = reader("5","first","sixth")
df_con5_pos7 = reader("5","first","seventh")
df_con5_pos8 = reader("5","first","eight")
df_con5_pos9 = reader("5","first","ninth")

fifth_cons = [df_con5_pos1,df_con5_pos2,df_con5_pos3,df_con5_pos4,df_con5_pos5,df_con5_pos6,df_con5_pos7,df_con5_pos8,df_con5_pos9]

first_con_data = pd.DataFrame({})

#print(first_seg_data.iloc[:,0])

counter = 0
for i, d in enumerate(first_cons):
        if i == 0:
                first_con_data.insert(loc=0,column="read_id",value=d["id"])
                first_con_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                first_con_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                first_con_data.insert(loc=3,column="pos1_base",value=d["base"])
                first_con_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                first_con_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        #if i == 1:
        #       first_seg_data.insert(loc=5,column="pos2_avg",value=d["current_avg"])
        #       first_seg_data.insert(loc=6,column="pos2_stdev",value=d["current_stdev"])
        #       first_seg_data.insert(loc=7,column="pos2_base",value=d["base"])
        #       first_seg_data.insert(loc=8,column="pos2_dwell",value=d["dwell"])
        #       counter+=4
        else:
                first_con_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                first_con_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                first_con_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                first_con_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                first_con_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5

first_con_data.insert(loc=counter,column="construct label",value="construct1")
first_con_data.insert(loc=counter,column="ino label",value="T")
print(len(first_con_data))
print(first_con_data.shape)

second_con_data = pd.DataFrame({})

counter = 0
for i, d in enumerate(second_cons):
        if i == 0:
                second_con_data.insert(loc=0,column="read_id",value=d["id"])
                second_con_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                second_con_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                second_con_data.insert(loc=3,column="pos1_base",value=d["base"])
                second_con_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                second_con_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        else:
                second_con_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                second_con_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                second_con_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                second_con_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                second_con_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5


second_con_data.insert(loc=counter,column="construct label",value="construct2")
second_con_data.insert(loc=counter,column="ino label",value="F")
print(len(second_con_data))
print(second_con_data.shape)

third_con_data = pd.DataFrame({})

counter = 0
for i, d in enumerate(third_cons):
        if i == 0:
                third_con_data.insert(loc=0,column="read_id",value=d["id"])
                third_con_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                third_con_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                third_con_data.insert(loc=3,column="pos1_base",value=d["base"])
                third_con_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                third_con_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        else:
                third_con_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                third_con_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                third_con_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                third_con_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                third_con_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5


third_con_data.insert(loc=counter,column="construct label",value="construct3")
third_con_data.insert(loc=counter,column="ino label",value="F")
print(len(third_con_data))
print(third_con_data.shape)

fourth_con_data = pd.DataFrame({})

counter = 0
for i, d in enumerate(fourth_cons):
        if i == 0:
                fourth_con_data.insert(loc=0,column="read_id",value=d["id"])
                fourth_con_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                fourth_con_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                fourth_con_data.insert(loc=3,column="pos1_base",value=d["base"])
                fourth_con_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                fourth_con_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        else:
                fourth_con_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                fourth_con_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                fourth_con_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                fourth_con_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                fourth_con_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5


fourth_con_data.insert(loc=counter,column="construct label",value="construct4")
fourth_con_data.insert(loc=counter,column="ino label",value="F")
print(len(fourth_con_data))
print(fourth_con_data.shape)

fifth_con_data = pd.DataFrame({})

counter = 0
for i, d in enumerate(fifth_cons):
        if i == 0:
                fifth_con_data.insert(loc=0,column="read_id",value=d["id"])
                fifth_con_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                fifth_con_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                fifth_con_data.insert(loc=3,column="pos1_base",value=d["base"])
                fifth_con_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                fifth_con_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        else:
                fifth_con_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                fifth_con_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                fifth_con_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                fifth_con_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                fifth_con_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5


fifth_con_data.insert(loc=counter,column="construct label",value="construct5")
fifth_con_data.insert(loc=counter,column="ino label",value="F")
print(len(fifth_con_data))
print(fifth_con_data.shape)

######################################################################
#
# Second Construct
#
######################################################################

df_con1_pos1 = reader("1","second","first")
df_con1_pos2 = reader("1","second","second")
df_con1_pos3 = reader("1","second","third")
df_con1_pos4 = reader("1","second","fourth")
df_con1_pos5 = reader("1","second","fifth")
df_con1_pos6 = reader("1","second","sixth")
df_con1_pos7 = reader("1","second","seventh")
df_con1_pos8 = reader("1","second","eight")
df_con1_pos9 = reader("1","second","ninth")

first_cons_second_seg = [df_con1_pos1,df_con1_pos2,df_con1_pos3,df_con1_pos4,df_con1_pos5,df_con1_pos6,df_con1_pos7,df_con1_pos8,df_con1_pos9]

#print(first_cons_second_seg[:10])

df_con2_pos1 = reader("2","second","first")
df_con2_pos2 = reader("2","second","second")
df_con2_pos3 = reader("2","second","third")
df_con2_pos4 = reader("2","second","fourth")
df_con2_pos5 = reader("2","second","fifth")
df_con2_pos6 = reader("2","second","sixth")
df_con2_pos7 = reader("2","second","seventh")
df_con2_pos8 = reader("2","second","eight")
df_con2_pos9 = reader("2","second","ninth")

second_cons_second_seg = [df_con2_pos1,df_con2_pos2,df_con2_pos3,df_con2_pos4,df_con2_pos5,df_con2_pos6,df_con2_pos7,df_con2_pos8,df_con2_pos9]

df_con3_pos1 = reader("3","second","first")
df_con3_pos2 = reader("3","second","second")
df_con3_pos3 = reader("3","second","third")
df_con3_pos4 = reader("3","second","fourth")
df_con3_pos5 = reader("3","second","fifth")
df_con3_pos6 = reader("3","second","sixth")
df_con3_pos7 = reader("3","second","seventh")
df_con3_pos8 = reader("3","second","eight")
df_con3_pos9 = reader("3","second","ninth")

third_cons_second_seg = [df_con3_pos1,df_con3_pos2,df_con3_pos3,df_con3_pos4,df_con3_pos5,df_con3_pos6,df_con3_pos7,df_con3_pos8,df_con3_pos9]

df_con4_pos1 = reader("4","second","first")
df_con4_pos2 = reader("4","second","second")
df_con4_pos3 = reader("4","second","third")
df_con4_pos4 = reader("4","second","fourth")
df_con4_pos5 = reader("4","second","fifth")
df_con4_pos6 = reader("4","second","sixth")
df_con4_pos7 = reader("4","second","seventh")
df_con4_pos8 = reader("4","second","eight")
df_con4_pos9 = reader("4","second","ninth")

fourth_cons_second_seg = [df_con4_pos1,df_con4_pos2,df_con4_pos3,df_con4_pos4,df_con4_pos5,df_con4_pos6,df_con4_pos7,df_con4_pos8,df_con4_pos9]

df_con5_pos1 = reader("5","second","first")
df_con5_pos2 = reader("5","second","second")
df_con5_pos3 = reader("5","second","third")
df_con5_pos4 = reader("5","second","fourth")
df_con5_pos5 = reader("5","second","fifth")
df_con5_pos6 = reader("5","second","sixth")
df_con5_pos7 = reader("5","second","seventh")
df_con5_pos8 = reader("5","second","eight")
df_con5_pos9 = reader("5","second","ninth")

fifth_cons_second_seg = [df_con5_pos1,df_con5_pos2,df_con5_pos3,df_con5_pos4,df_con5_pos5,df_con5_pos6,df_con5_pos7,df_con5_pos8,df_con5_pos9]

first_con_sec_seg_data = pd.DataFrame({})

#print(first_seg_data.iloc[:,0])

counter = 0
for i, d in enumerate(first_cons_second_seg):
        if i == 0:
                first_con_sec_seg_data.insert(loc=0,column="read_id",value=d["id"])
                first_con_sec_seg_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                first_con_sec_seg_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                first_con_sec_seg_data.insert(loc=3,column="pos1_base",value=d["base"])
                first_con_sec_seg_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                first_con_sec_seg_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        else:
                first_con_sec_seg_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                first_con_sec_seg_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                first_con_sec_seg_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                first_con_sec_seg_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                first_con_sec_seg_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5

first_con_sec_seg_data.insert(loc=counter,column="construct label",value="construct1")
first_con_sec_seg_data.insert(loc=counter,column="ino label",value="T")
print(len(first_con_sec_seg_data))
print(first_con_sec_seg_data.shape)

second_con_sec_seg_data = pd.DataFrame({})

counter = 0
for i, d in enumerate(second_cons_second_seg):
        if i == 0:
                second_con_sec_seg_data.insert(loc=0,column="read_id",value=d["id"])
                second_con_sec_seg_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                second_con_sec_seg_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                second_con_sec_seg_data.insert(loc=3,column="pos1_base",value=d["base"])
                second_con_sec_seg_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                second_con_sec_seg_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        else:
                second_con_sec_seg_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                second_con_sec_seg_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                second_con_sec_seg_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                second_con_sec_seg_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                second_con_sec_seg_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5


second_con_sec_seg_data.insert(loc=counter,column="construct label",value="construct2")
second_con_sec_seg_data.insert(loc=counter,column="ino label",value="F")
print(len(second_con_sec_seg_data))
print(second_con_sec_seg_data.shape)

third_con_sec_seg_data = pd.DataFrame({})

counter = 0
for i, d in enumerate(third_cons_second_seg):
        if i == 0:
                third_con_sec_seg_data.insert(loc=0,column="read_id",value=d["id"])
                third_con_sec_seg_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                third_con_sec_seg_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                third_con_sec_seg_data.insert(loc=3,column="pos1_base",value=d["base"])
                third_con_sec_seg_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                third_con_sec_seg_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        else:
                third_con_sec_seg_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                third_con_sec_seg_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                third_con_sec_seg_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                third_con_sec_seg_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                third_con_sec_seg_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5


third_con_sec_seg_data.insert(loc=counter,column="construct label",value="construct3")
third_con_sec_seg_data.insert(loc=counter,column="ino label",value="F")
print(len(third_con_sec_seg_data))
print(third_con_sec_seg_data.shape)

fourth_con_sec_seg_data = pd.DataFrame({})

counter = 0
for i, d in enumerate(fourth_cons_second_seg):
        if i == 0:
                fourth_con_sec_seg_data.insert(loc=0,column="read_id",value=d["id"])
                fourth_con_sec_seg_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                fourth_con_sec_seg_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                fourth_con_sec_seg_data.insert(loc=3,column="pos1_base",value=d["base"])
                fourth_con_sec_seg_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                fourth_con_sec_seg_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        else:
                fourth_con_sec_seg_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                fourth_con_sec_seg_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                fourth_con_sec_seg_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                fourth_con_sec_seg_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                fourth_con_sec_seg_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5


fourth_con_sec_seg_data.insert(loc=counter,column="construct label",value="construct4")
fourth_con_sec_seg_data.insert(loc=counter,column="ino label",value="F")
print(len(fourth_con_sec_seg_data))
print(fourth_con_sec_seg_data.shape)

fifth_con_sec_seg_data = pd.DataFrame({})

counter = 0
for i, d in enumerate(fifth_cons_second_seg):
        if i == 0:
                fifth_con_sec_seg_data.insert(loc=0,column="read_id",value=d["id"])
                fifth_con_sec_seg_data.insert(loc=1,column="pos1_avg",value=d["current_avg"])
                fifth_con_sec_seg_data.insert(loc=2,column="pos1_stdev",value=d["current_stdev"])
                fifth_con_sec_seg_data.insert(loc=3,column="pos1_base",value=d["base"])
                fifth_con_sec_seg_data.insert(loc=4,column="pos1_dwell",value=d["dwell"])
                fifth_con_sec_seg_data.insert(loc=5,column="pos1_score",value=d["score"])
                counter+=6
        else:
                fifth_con_sec_seg_data.insert(loc=counter,column="pos"+str(i+1)+"_avg",value=d["current_avg"])
                fifth_con_sec_seg_data.insert(loc=counter+1,column="pos"+str(i+1)+"_stdev",value=d["current_stdev"])
                fifth_con_sec_seg_data.insert(loc=counter+2,column="pos"+str(i+1)+"_base",value=d["base"])
                fifth_con_sec_seg_data.insert(loc=counter+3,column="pos"+str(i+1)+"_dwell",value=d["dwell"])
                fifth_con_sec_seg_data.insert(loc=counter+4,column="pos"+str(i+1)+"_score",value=d["score"])
                counter+=5


fifth_con_sec_seg_data.insert(loc=counter,column="construct label",value="construct5")
fifth_con_sec_seg_data.insert(loc=counter,column="ino label",value="F")
print(len(fifth_con_sec_seg_data))
print(fifth_con_sec_seg_data.shape)

all_constructs = pd.concat([first_con_data,second_con_data,third_con_data,fourth_con_data,fifth_con_data])
all_constructs_hf = h2o.H2OFrame(all_constructs)

all_constructs_sec_seg = pd.concat([first_con_sec_seg_data,second_con_sec_seg_data,third_con_sec_seg_data,fourth_con_sec_seg_data,fifth_con_sec_seg_data])
all_constructs_sec_seg_hf = h2o.H2OFrame(all_constructs_sec_seg)

print("First:\n")
print(all_constructs.head())
print("Second:\n")
print(all_constructs_sec_seg.head())

#all_constructs_sec_seg.describe()
#all_constructs_sec_seg_hf.describe()
'''
predictors = ["pos1_avg","pos1_base","pos1_dwell", "pos1_score",
              "pos2_avg","pos2_base","pos2_dwell", "pos2_score",
              "pos3_avg","pos3_base","pos3_dwell", "pos3_score",
              "pos4_avg","pos4_base","pos4_dwell", "pos4_score",
              "pos5_avg","pos5_base","pos5_dwell", "pos5_score",
              "pos6_avg","pos6_base","pos6_dwell", "pos6_score",
              "pos7_avg","pos7_base","pos7_dwell", "pos7_score",
              "pos8_avg","pos8_base","pos8_dwell", "pos8_score",
              "pos9_avg","pos9_base","pos9_dwell", "pos9_score"]
'''

predictors = ["pos1_avg","pos1_dwell", "pos1_score",
              "pos2_avg","pos2_dwell", "pos2_score",
              "pos3_avg","pos3_dwell", "pos3_score",
              "pos4_avg","pos4_dwell", "pos4_score",
              "pos5_avg","pos5_dwell", "pos5_score",
              "pos6_avg","pos6_dwell", "pos6_score",
              "pos7_avg","pos7_dwell", "pos7_score",
              "pos8_avg","pos8_dwell", "pos8_score",
              "pos9_avg","pos9_dwell", "pos9_score"]

#predictors = ["pos5_avg","pos5_base","pos5_dwell","pos5_score"]
response = "construct label"

# Split the dataset into a train and valid set:
train, valid = all_constructs_hf.split_frame(ratios=[.8], seed=1234)
sec_train, sec_valid = all_constructs_sec_seg_hf.split_frame(ratios=[.8], seed=1234)

# Build and train the model:
first_drf = H2ORandomForestEstimator(ntrees=20,
                                    max_depth=7,
                                    min_rows=10,
                                    calibrate_model=False,
				    calibration_frame=valid,
                                    balance_classes=True,
                                    binomial_double_trees=True)
first_drf.train(x=predictors,
               y=response,
               training_frame=train,
               validation_frame=sec_valid)

perf = first_drf.model_performance()
print(perf)
