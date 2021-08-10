import csv
import pandas as pd
from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, Whisker, SingleIntervalTicker

################################################################################
#
# This script reads data regarding the dwell time for each of the  9 positions 
# we are interested. This is specified for the first of the three segments. 
# These are then plotted in a line graph with bokeh.
#
################################################################################

con1_first = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/firstData.csv")
con1_second = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/secondData.csv")
con1_third = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/thirdData.csv")
con1_fourth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/fourthData.csv")
con1_fifth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/fifthData.csv")
con1_sixth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/sixthData.csv")
con1_seventh = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/seventhData.csv")
con1_eight = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/eightData.csv")
con1_ninth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct1/firstPos/ninthData.csv")

con2_first = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct2/firstPos/firstData.csv")
con2_second = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct2/firstPos/secondData.csv")
con2_third = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct2/firstPos/thirdData.csv")
con2_fourth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct2/firstPos/fourthData.csv")
con2_fifth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct2/firstPos/fifthData.csv")
con2_sixth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct2/firstPos/sixthData.csv")
con2_seventh = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct2/firstPos/seventhData.csv")
con2_eight = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct2/firstPos/eightData.csv")
con2_ninth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct2/firstPos/ninthData.csv")

con3_first = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct3/firstPos/firstData.csv")
con3_second = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct3/firstPos/secondData.csv")
con3_third = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct3/firstPos/thirdData.csv")
con3_fourth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct3/firstPos/fourthData.csv")
con3_fifth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct3/firstPos/fifthData.csv")
con3_sixth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct3/firstPos/sixthData.csv")
con3_seventh = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct3/firstPos/seventhData.csv")
con3_eight = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct3/firstPos/eightData.csv")
con3_ninth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct3/firstPos/ninthData.csv")

con4_first = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/firstPos/firstData.csv")
con4_second = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/firstPos/secondData.csv")
con4_third = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/firstPos/thirdData.csv")
con4_fourth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/firstPos/fourthData.csv")
con4_fifth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/firstPos/fifthData.csv")
con4_sixth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/firstPos/sixthData.csv")
con4_seventh = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/firstPos/seventhData.csv")
con4_eight = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/firstPos/eightData.csv")
con4_ninth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct4/firstPos/ninthData.csv")

con5_first = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct5/firstPos/firstData.csv")
con5_second = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct5/firstPos/secondData.csv")
con5_third = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct5/firstPos/thirdData.csv")
con5_fourth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct5/firstPos/fourthData.csv")
con5_fifth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct5/firstPos/fifthData.csv")
con5_sixth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct5/firstPos/sixthData.csv")
con5_seventh = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct5/firstPos/seventhData.csv")
con5_eight = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct5/firstPos/eightData.csv")
con5_ninth = pd.read_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/further_analysis/construct5/firstPos/ninthData.csv")

#print(first.head())
#print("Mean of first", first["current_avg"].mean())
#first["var"] = first["current_stdev"]**2
#print("Stdev of first", math.sqrt(first["var"].mean()))
#print("Mean dwell time",first["dwell"].mean())
#print(first.groupby("base").count())

con1_firstPoint = [1,con1_first["dwell"].mean()]
con1_secondPoint = [2,con1_second["dwell"].mean()]
con1_thirdPoint = [3,con1_third["dwell"].mean()]
con1_fourthPoint = [4,con1_fourth["dwell"].mean()]
con1_fifthPoint = [5,con1_fifth["dwell"].mean()]
con1_sixthPoint = [6,con1_sixth["dwell"].mean()]
con1_seventhPoint = [7,con1_seventh["dwell"].mean()]
con1_eightPoint = [8,con1_eight["dwell"].mean()]
con1_ninthPoint = [9,con1_ninth["dwell"].mean()]

con2_firstPoint = [1,con2_first["dwell"].mean()]
con2_secondPoint = [2,con2_second["dwell"].mean()]
con2_thirdPoint = [3,con2_third["dwell"].mean()]
con2_fourthPoint = [4,con2_fourth["dwell"].mean()]
con2_fifthPoint = [5,con2_fifth["dwell"].mean()]
con2_sixthPoint = [6,con2_sixth["dwell"].mean()]
con2_seventhPoint = [7,con2_seventh["dwell"].mean()]
con2_eightPoint = [8,con2_eight["dwell"].mean()]
con2_ninthPoint = [9,con2_ninth["dwell"].mean()]

con3_firstPoint = [1,con3_first["dwell"].mean()]
con3_secondPoint = [2,con3_second["dwell"].mean()]
con3_thirdPoint = [3,con3_third["dwell"].mean()]
con3_fourthPoint = [4,con3_fourth["dwell"].mean()]
con3_fifthPoint = [5,con3_fifth["dwell"].mean()]
con3_sixthPoint = [6,con3_sixth["dwell"].mean()]
con3_seventhPoint = [7,con3_seventh["dwell"].mean()]
con3_eightPoint = [8,con3_eight["dwell"].mean()]
con3_ninthPoint = [9,con3_ninth["dwell"].mean()]

con4_firstPoint = [1,con4_first["dwell"].mean()]
con4_secondPoint = [2,con4_second["dwell"].mean()]
con4_thirdPoint = [3,con4_third["dwell"].mean()]
con4_fourthPoint = [4,con4_fourth["dwell"].mean()]
con4_fifthPoint = [5,con4_fifth["dwell"].mean()]
con4_sixthPoint = [6,con4_sixth["dwell"].mean()]
con4_seventhPoint = [7,con4_seventh["dwell"].mean()]
con4_eightPoint = [8,con4_eight["dwell"].mean()]
con4_ninthPoint = [9,con4_ninth["dwell"].mean()]

con5_firstPoint = [1,con5_first["dwell"].mean()]
con5_secondPoint = [2,con5_second["dwell"].mean()]
con5_thirdPoint = [3,con5_third["dwell"].mean()]
con5_fourthPoint = [4,con5_fourth["dwell"].mean()]
con5_fifthPoint = [5,con5_fifth["dwell"].mean()]
con5_sixthPoint = [6,con5_sixth["dwell"].mean()]
con5_seventhPoint = [7,con5_seventh["dwell"].mean()]
con5_eightPoint = [8,con5_eight["dwell"].mean()]
con5_ninthPoint = [9,con5_ninth["dwell"].mean()]

con1_firstStd = [1,con1_first["dwell"].std()]
con1_secondStd = [2,con1_second["dwell"].std()]
con1_thirdStd = [3,con1_third["dwell"].std()]
con1_fourthStd = [4,con1_fourth["dwell"].std()]
con1_fifthStd = [5,con1_fifth["dwell"].std()]
con1_sixthStd = [6,con1_sixth["dwell"].std()]
con1_seventhStd = [7,con1_seventh["dwell"].std()]
con1_eightStd = [8,con1_eight["dwell"].std()]
con1_ninthStd = [9,con1_ninth["dwell"].std()]

con2_firstStd = [1,con2_first["dwell"].std()]
con2_secondStd = [2,con2_second["dwell"].std()]
con2_thirdStd = [3,con2_third["dwell"].std()]
con2_fourthStd = [4,con2_fourth["dwell"].std()]
con2_fifthStd = [5,con2_fifth["dwell"].std()]
con2_sixthStd = [6,con2_sixth["dwell"].std()]
con2_seventhStd = [7,con2_seventh["dwell"].std()]
con2_eightStd = [8,con2_eight["dwell"].std()]
con2_ninthStd = [9,con2_ninth["dwell"].std()]

con3_firstStd = [1,con3_first["dwell"].std()]
con3_secondStd = [2,con3_second["dwell"].std()]
con3_thirdStd = [3,con3_third["dwell"].std()]
con3_fourthStd = [4,con3_fourth["dwell"].std()]
con3_fifthStd = [5,con3_fifth["dwell"].std()]
con3_sixthStd = [6,con3_sixth["dwell"].std()]
con3_seventhStd = [7,con3_seventh["dwell"].std()]
con3_eightStd = [8,con3_eight["dwell"].std()]
con3_ninthStd = [9,con3_ninth["dwell"].std()]

con4_firstStd = [1,con4_first["dwell"].std()]
con4_secondStd = [2,con4_second["dwell"].std()]
con4_thirdStd = [3,con4_third["dwell"].std()]
con4_fourthStd = [4,con4_fourth["dwell"].std()]
con4_fifthStd = [5,con4_fifth["dwell"].std()]
con4_sixthStd = [6,con4_sixth["dwell"].std()]
con4_seventhStd = [7,con4_seventh["dwell"].std()]
con4_eightStd = [8,con4_eight["dwell"].std()]
con4_ninthStd = [9,con4_ninth["dwell"].std()]

con5_firstStd = [1,con5_first["dwell"].std()]
con5_secondStd = [2,con5_second["dwell"].std()]
con5_thirdStd = [3,con5_third["dwell"].std()]
con5_fourthStd = [4,con5_fourth["dwell"].std()]
con5_fifthStd = [5,con5_fifth["dwell"].std()]
con5_sixthStd = [6,con5_sixth["dwell"].std()]
con5_seventhStd = [7,con5_seventh["dwell"].std()]
con5_eightStd = [8,con5_eight["dwell"].std()]
con5_ninthStd = [9,con5_ninth["dwell"].std()]

con1_x = [con1_firstPoint[0], con1_secondPoint[0], con1_thirdPoint[0], con1_fourthPoint[0], con1_fifthPoint[0], con1_sixthPoint[0], con1_seventhPoint[0], con1_eightPoint[0], con1_ninthPoint[0]]
con1_y = [con1_firstPoint[1], con1_secondPoint[1], con1_thirdPoint[1], con1_fourthPoint[1], con1_fifthPoint[1], con1_sixthPoint[1], con1_seventhPoint[1], con1_eightPoint[1], con1_ninthPoint[1]]

con2_x = [con2_firstPoint[0], con2_secondPoint[0], con2_thirdPoint[0], con2_fourthPoint[0], con2_fifthPoint[0], con2_sixthPoint[0], con2_seventhPoint[0], con2_eightPoint[0], con2_ninthPoint[0]]
con2_y = [con2_firstPoint[1], con2_secondPoint[1], con2_thirdPoint[1], con2_fourthPoint[1], con2_fifthPoint[1], con2_sixthPoint[1], con2_seventhPoint[1], con2_eightPoint[1], con2_ninthPoint[1]]

con3_x = [con3_firstPoint[0], con3_secondPoint[0], con3_thirdPoint[0], con3_fourthPoint[0], con3_fifthPoint[0], con3_sixthPoint[0], con3_seventhPoint[0], con3_eightPoint[0], con3_ninthPoint[0]]
con3_y = [con3_firstPoint[1], con3_secondPoint[1], con3_thirdPoint[1], con3_fourthPoint[1], con3_fifthPoint[1], con3_sixthPoint[1], con3_seventhPoint[1], con3_eightPoint[1], con3_ninthPoint[1]]

con4_x = [con4_firstPoint[0], con4_secondPoint[0], con4_thirdPoint[0], con4_fourthPoint[0], con4_fifthPoint[0], con4_sixthPoint[0], con4_seventhPoint[0], con4_eightPoint[0], con4_ninthPoint[0]]
con4_y = [con4_firstPoint[1], con4_secondPoint[1], con4_thirdPoint[1], con4_fourthPoint[1], con4_fifthPoint[1], con4_sixthPoint[1], con4_seventhPoint[1], con4_eightPoint[1], con4_ninthPoint[1]]

con5_x = [con5_firstPoint[0], con5_secondPoint[0], con5_thirdPoint[0], con5_fourthPoint[0], con5_fifthPoint[0], con5_sixthPoint[0], con5_seventhPoint[0], con5_eightPoint[0], con5_ninthPoint[0]]
con5_y = [con5_firstPoint[1], con5_secondPoint[1], con5_thirdPoint[1], con5_fourthPoint[1], con5_fifthPoint[1], con5_sixthPoint[1], con5_seventhPoint[1], con5_eightPoint[1], con5_ninthPoint[1]]

con1_std_x = [con1_firstStd[0], con1_secondStd[0], con1_thirdStd[0], con1_fourthStd[0], con1_fifthStd[0], con1_sixthStd[0], con1_seventhStd[0], con1_eightStd[0], con1_ninthStd[0]]
con1_std_y = [con1_firstStd[1], con1_secondStd[1], con1_thirdStd[1], con1_fourthStd[1], con1_fifthStd[1], con1_sixthStd[1], con1_seventhStd[1], con1_eightStd[1], con1_ninthStd[1]]

con2_std_x = [con2_firstStd[0], con2_secondStd[0], con2_thirdStd[0], con2_fourthStd[0], con2_fifthStd[0], con2_sixthStd[0], con2_seventhStd[0], con2_eightStd[0], con2_ninthStd[0]]
con2_std_y = [con2_firstStd[1], con2_secondStd[1], con2_thirdStd[1], con2_fourthStd[1], con2_fifthStd[1], con2_sixthStd[1], con2_seventhStd[1], con2_eightStd[1], con2_ninthStd[1]]

con3_std_x = [con3_firstStd[0], con3_secondStd[0], con3_thirdStd[0], con3_fourthStd[0], con3_fifthStd[0], con3_sixthStd[0], con3_seventhStd[0], con3_eightStd[0], con3_ninthStd[0]]
con3_std_y = [con3_firstStd[1], con3_secondStd[1], con3_thirdStd[1], con3_fourthStd[1], con3_fifthStd[1], con3_sixthStd[1], con3_seventhStd[1], con3_eightStd[1], con3_ninthStd[1]]

con4_std_x = [con4_firstStd[0], con4_secondStd[0], con4_thirdStd[0], con4_fourthStd[0], con4_fifthStd[0], con4_sixthStd[0], con4_seventhStd[0], con4_eightStd[0], con4_ninthStd[0]]
con4_std_y = [con4_firstStd[1], con4_secondStd[1], con4_thirdStd[1], con4_fourthStd[1], con4_fifthStd[1], con4_sixthStd[1], con4_seventhStd[1], con4_eightStd[1], con4_ninthStd[1]]

con5_std_x = [con5_firstStd[0], con5_secondStd[0], con5_thirdStd[0], con5_fourthStd[0], con5_fifthStd[0], con5_sixthStd[0], con5_seventhStd[0], con5_eightStd[0], con5_ninthStd[0]]
con5_std_y = [con5_firstStd[1], con5_secondStd[1], con5_thirdStd[1], con5_fourthStd[1], con5_fifthStd[1], con5_sixthStd[1], con5_seventhStd[1], con5_eightStd[1], con5_ninthStd[1]]

####OPTIONALLY ADD WHISKERS
#base, lower, upper = [],[],[]

#for i,j,df in zip(con1_x,con1_y,[con1_first,con1_second,con1_third,con1_fourth,con1_fifth,con1_sixth,con1_seventh,con1_eight,con1_ninth]):
#	base.append(i)
#	mean = j
#	lower.append(mean-df["current_avg"].std())
#	upper.append(mean+df["current_avg"].std())
	
#source_error = ColumnDataSource(data=dict(base=base,lower=lower,upper=upper))

p = figure(plot_width=1000, plot_height=400,x_axis_label='Base',y_axis_label='Mean Dwell Time')

p.xaxis.ticker = SingleIntervalTicker(interval=1)
p.line(con1_x, con1_y, line_width=4, line_color='#173F5F', legend_label="Inosine")
p.line(con2_x, con2_y, line_width=2, line_color='#66ff00', legend_label="Adenine")
p.line(con3_x, con3_y, line_width=2, line_color='#3CAEA3', legend_label="Thymine")
p.line(con4_x, con4_y, line_width=2, line_color='#F6D55C', legend_label="Cytosine")
p.line(con5_x, con5_y, line_width=2, line_color='#ED553B', legend_label="Guanine")

p2 = figure(plot_width=1000, plot_height=400,x_axis_label='Base',y_axis_label='Dwell Time Standard Deviation')
p2.xaxis.ticker = SingleIntervalTicker(interval=1)

p2.line(con1_std_x, con1_std_y, line_width=2, line_color='#173F5F', legend_label="Inosine")
p2.line(con2_std_x, con2_std_y, line_width=2, line_color='#20639B', legend_label="Adenine")
p2.line(con3_std_x, con3_std_y, line_width=2, line_color='#3CAEA3', legend_label="Thymine")
p2.line(con4_std_x, con4_std_y, line_width=2, line_color='#F6D55C', legend_label="Cytosine")
p2.line(con5_std_x, con5_std_y, line_width=2, line_color='#ED553B', legend_label="Guanine")

#p.add_layout(Whisker(source=source_error, base="base", lower="lower", upper="upper"))

save(p)
#save(p2)

