from ont_fast5_api.fast5_interface import get_fast5_file
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pandas as pd
import numpy as np

################################################################################
#
# This script reads a set of .fast5 files. For each such file it aligns to a 
# segment of interest. It then finds the raw data corresponding to this segment
# and calculates and gathers descriptive statistics from this segment. In a 
# dataframe, this is stored and written out. 
#
################################################################################

def print_all_raw_data(path,out_df):
	fast5_filepath = path # This can be a single- or multi-read file
	with get_fast5_file(fast5_filepath, mode="r") as f5:
		for idx, read in enumerate(f5.get_reads()):
			#if idx == 0:
				#continue
			#if idx == 10:
			#	exit()
			if idx%80 == 0:
				print("On read",idx,"of 4000")
			f5_rh = f5.get_read(read.read_id)
			#if analysis_group:
			#	latest_analysis = analysis_group
			#else:
			latest_analysis = f5_rh.get_latest_analysis("Basecall_1D")
			if latest_analysis == None:
				continue
			seg_idx = latest_analysis[-3:]
			seg_path = f"Segmentation_{seg_idx}/Summary/segmentation"

			summary_data = f5_rh.get_summary_data(latest_analysis)
			basecall_template = f5_rh.get_chain(latest_analysis)[0][0].decode("utf-8") + "_template"
						
			raw_data = read.get_raw_data()
			summary = read.get_channel_info()
			seqInfo = read.get_analysis_dataset("Basecall_1D_000","BaseCalled_template")['Fastq'][()]
			moves = read.get_analysis_dataset("Basecall_1D_000","BaseCalled_template")['Move'][()]
			
			#Relevant read specific information
			stride = summary_data[basecall_template]["block_stride"]
			called_ev = summary_data[basecall_template]["called_events"]
			first_sample_template = f5_rh.get_analysis_attributes(seg_path)['first_sample_template']
			duration_template = f5_rh.get_analysis_attributes(seg_path)['duration_template']			
			num_events = f5_rh.get_analysis_attributes(seg_path)['num_events_template']
			fastq = str(seqInfo).split("\\n")
			offset = summary["offset"]	
			scale = summary["range"]/summary["digitisation"]
			
			#print("aligning")
			align = pairwise2.align.localms(fastq[1],"GTCAGTGCA",1,-1,-1,-1) #for first pos: GTCANTGCA, for second pos: AAGANGTAC, for third pos: TTCANCTGT, use GTTCGACTG for construct 5
			#print(format_alignment(*align[0]))
			
			if align[0].score < 9 or align[0].end-align[0].start != 9: #Set lower score to 9 for canonical base constructs, and to 7 for Inosine constructs
				#print(read.read_id)
				#print("skipping low score")
				#print("This is read sequence:",fastq[1], "it has length",len(fastq[1]))
				continue

			#print(read.read_id, len(raw_data))
			#print("This is read sequence:",fastq[1], "it has length",len(fastq[1]))
			#print("\nSome read info\n")
			#print("Stride:",stride,"\n")
			#print("Sample temp:",first_sample_template,"\n")
			#print("Duration temp:", duration_template,"\n")
			#print("Number of events:", num_events,"\n")
			#print("Called events:", called_ev,"\n")
			df = pd.DataFrame({"base": [], "move": [], "raw_index":[]})
			i_segment = pd.DataFrame({"base": [], "move": [], "raw_index":[]})
			moveSum = 0
			for idx, num in enumerate(moves):
				raw_idx = (first_sample_template+idx)*5
				if num == 1:
					if moveSum >= align[0].start and moveSum < align[0].end:
						moveSum+=1
						newRow = {'base': fastq[1][moveSum-1],'move': 1, "raw_index": raw_idx}
						i_segment = i_segment.append(newRow, ignore_index=True)
						df = df.append(newRow, ignore_index=True)
					#elif moveSum >= len(fastq[1]):
					#	newRow = {'base': fastq[1][moveSum-1],'move': 1, "raw_index": raw_idx}
					#	print("Have base:",fastq[1][moveSum-1])
					#	df = df.append(newRow, ignore_index=True)
					#	continue
					else:
						moveSum+=1
						newRow = {'base': fastq[1][moveSum-1],'move': 1, "raw_index": raw_idx}
						df = df.append(newRow, ignore_index=True)
				else:
					if moveSum > align[0].start and moveSum < align[0].end:
						newRow = {'base': fastq[1][moveSum-1],'move': 0, "raw_index": raw_idx}
						i_segment = i_segment.append(newRow, ignore_index=True)
					newRow = {'base': fastq[1][moveSum-1], 'move': 0, "raw_index": raw_idx}
					df = df.append(newRow, ignore_index=True)				
			#print(df.head())
			#print(df.tail())
			#print(i_segment)
			#print(align[0].end-align[0].start)
			#print(fastq[1][align[0].start:align[0].end])		
		
			startRawIndex = int(i_segment["raw_index"].loc[0])
			endRawIndex = int(i_segment["raw_index"].loc[len(i_segment)-1]+4)
			
			#print(raw_data[startRawIndex:endRawIndex])
			norm_data = scale * (raw_data[startRawIndex:endRawIndex] + offset)
			#print(len(norm_data))
			
			count = 0
			for index, row in i_segment.iterrows():
				if index == 0: 
					current_base = row["base"]
					start = row["raw_index"]
					offset = start
				if index == len(i_segment)-1:
					if current_base != row["base"]:
						stop = row["raw_index"]-1
						#print(stop,"and",start)
						#print("base", current_base,"has average of", np.mean(norm_data[int(start-offset):int(stop-offset)]), "and", count,"events.")
						newRow = {"read_id":read.read_id,"base":current_base,"signal_avg":np.mean(norm_data[int(start-offset):int(stop-offset)]),
								"signal_stdev":np.std(norm_data[int(start-offset):int(stop-offset)]),"dwell_time":count*stride}
						out_df = out_df.append(newRow, ignore_index=True)
						current_base = row["base"]
						start = row["raw_index"]
						count = 1
						#print("here with",current_base,"and",count)
					stop = row["raw_index"]+4
					#print(stop,"and",start)
					#print("base", current_base,"has average of", np.mean(norm_data[int(start-offset):int(stop-offset)]), "and", count,"events.")
					newRow = {"read_id":read.read_id,"base":current_base,"signal_avg":np.mean(norm_data[int(start-offset):int(stop-offset)]),
							"signal_stdev":np.std(norm_data[int(start-offset):int(stop-offset)]), "dwell_time":count*stride}
					out_df = out_df.append(newRow, ignore_index=True)
					#print("Last base is",row["base"])
					count = 0
					continue			
				elif current_base != row["base"] or row["move"] == 1.0 and index != 0:
					stop = row["raw_index"]-1
					#print(stop, "and", start)
					#print("base", current_base,"has average of", np.mean(norm_data[int(start-offset):int(stop-offset)]), "and", count,"events.")
					newRow = {"read_id":read.read_id,"base":current_base,"signal_avg":np.mean(norm_data[int(start-offset):int(stop-offset)]),
							"signal_stdev":np.std(norm_data[int(start-offset):int(stop-offset)]),"dwell_time":count*stride}
					out_df = out_df.append(newRow, ignore_index=True)
					#print("New base is",row["base"])
					current_base = row["base"]
					start = row["raw_index"]				
					count = 0
				count +=1
			#print(out_df)
			
		return out_df			

total_df = pd.DataFrame({"read_id":[],"base":[],"signal_avg":[],"signal_stdev":[],"dwell_time":[]})
with open("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/construct5/another_run/firstGuaBasePath.txt") as f:
	lines = f.readlines()
lines = [x.strip() for x in lines]
for idx,line in enumerate(lines):
	#if idx+1 >= 27:
	#	print("here")
	#	continue
	print("Processing file",idx+1,"out of",len(lines))
	#if idx == 1:
	#	break
	result_df = pd.DataFrame({"read_id":[],"base":[],"signal_avg":[],"signal_stdev":[],"dwell_time":[]})
	results_df = print_all_raw_data(line,result_df)
	total_df = total_df.append(results_df, ignore_index=True)

print(len(total_df))
total_df.to_csv("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/4_fast5_analysis/construct5/another_run/perfs/first_pos_data_out.csv")
