import csv
import sys

################################################################################
#
# This script reads the output from ino_finder.py and extracts only the read ids
# included. These are all the reads we are interested in. 
#
################################################################################

csv.field_size_limit(sys.maxsize)

print("Opening in and out files.\n")
with open("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/3_tsv_analysis/construct2/third_adenines.tsv","rt") as src, open("construct2/thirdPerfAdenineReads.tsv","wt") as out:
	csrc = csv.reader(src, delimiter='\t')
	wout = csv.writer(out, delimiter='\t')
	newRead = False
	currentRead = ''
	for num, row in enumerate(csrc):
		if currentRead != row[0]:
			wout.writerow([row[0]])
			currentRead = row[0]
	print("Wrote all files to separate file")	
