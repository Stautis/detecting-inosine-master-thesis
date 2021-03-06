import csv
import pandas as pd

################################################################################
#
# This script reads a specified .tsv file and edits the .tsv file to contain 
# only relevant information. In the context of the second experiment, this means
# only information regarding the last 69 bases of the construct. A modified 
# version of this script was also used to try to prune the eventalign_collapse
# output directly, but when that didn't work this approach was taken to prune
# the eventalign file prior to collapsing.
#
################################################################################

totalLines = 165032463
print("Opening in and out files.\n")
with open("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/3_event_align/construct1/A_in_reference/A_ref_eventalign.tsv","rt") as src, open("pruned.tsv","wt") as out:
	csrc = csv.reader(src, delimiter='\t')
	wout = csv.writer(out, delimiter='\t')
	#header = []
	cols = []
	newRead = False
	print("reading .tsv file line by line")
	for num, row in enumerate(csrc):
		if not (row[1].isdigit()):
			newRead = True
			#if len(str(row[0]).split("#")) > 1:
			#	header = row
			#else:
			cols = row
			continue
		if int(row[1]) >= 152:
			if newRead:
				newRead=False
				#wout.writerow(header)
				wout.writerow(cols)
			wout.writerow(row)	
		if (num%1000000)==0:
			print("read",num,"of",totalLines)
