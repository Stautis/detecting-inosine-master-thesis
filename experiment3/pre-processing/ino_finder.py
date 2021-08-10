import csv
import sys
import pandas as pd

################################################################################
#
# This script reads a specified .tsv file, then it looks through the reads 
# contained in the file. Looking for full segments matching the three 9-base
# segments of interest. For any given segment to be included in further analysis
# it must contain a match for the central variable positions and the four 
# constant bases on either side. Reference-relative indices are used to minimize
# the number of rows checked, thus these are specifically hard-coded for each of 
# the three variable sites in the five different constructs. 
#
################################################################################

def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

csv.field_size_limit(sys.maxsize)

print("Opening in and out files.\n")
with open("/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/2_tsv_conversion/construct2/construct2_test.tsv","rt") as src, open("construct2/third_adenines.tsv","wt") as out:
	csrc = csv.reader(src, delimiter='\t')
	wout = csv.writer(out, delimiter='\t')
	print("reading .tsv file line by line")
	newRead = False
	skipRead = False
	currentRead = ''
	currentFlag = ''
	readCount = 0
	rowsPrinted = 0
	readsPrinted = 0
	rowsToPrint = []
	totalLines = 114464213 
	for num, row in enumerate(csrc):
		if num == 0:
			continue
		if num%1000000 == 0:
			print("read",(num/totalLines)*100,"% of lines")
		if currentRead == row[0] and skipRead == True:
			continue
		if currentRead != row[0]:
			skipRead = False
			#print(rowsPrinted,"relevant rows in last read")
			if rowsPrinted < 9 and readCount > 0:
				#print("Last read did not have 9 entries between index 166 and 174, and this was",readCount)
				if rowsPrinted > 0:
					#print("rows in previous read",rowsPrinted)
					rowsToPrint = []
					currentRead = row[0]
					currentFlag = row[1]
					rowsPrinted = 0
					readCount += 1		
					continue
			elif rowsPrinted == 9 and currentFlag != '2048':
				#print("Printing rows from previous read to file\n")
				#print("Current read is",row[0])
				for r in rowsToPrint:
					#print(r)
					wout.writerow(r)
				readsPrinted += 1
			rowsToPrint = []	
			#print("Have a new read\n")
			currentRead = row[0]
			print("Here with",currentRead)
			currentFlag = row[1]
			rowsPrinted = 0
			readCount += 1
		elif (currentRead == row[0] and currentFlag != row[1]) or currentFlag == '2048' or len(row) <= 7:
			continue
		#if readCount == 3:
			#exit()
		if RepresentsInt(row[7]): 
			if (int(row[7]) >= 261 and int(row[7]) <= 269): #First Region: 166 - 174, Second Region: 180 - 188, Third Region: 195 - 203
				if row[9] != 'M':
					skipRead = True
					rowsToPrint = []
					continue
				rowsToPrint.append(row)
				rowsPrinted+=1
		else:
			continue
	print("printed a total of",readsPrinted,"reads")
