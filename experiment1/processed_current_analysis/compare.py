import sys

################################################################################
#
# This script prints data from two identical positions from eventalign outputs
# generated form the same set of reads but with different references. The printed
# results can then be compared. 
#
################################################################################

linesOne = open("run_ref_3/aoi.txt").read().splitlines()
linesTwo = open("run_ref_4/aoi.txt").read().splitlines()
currentId = ''

def searchSub(id, pos):
	#print("\nSearching for", id,"in run_ref_4")
	matchLine = ''
	for num, line in enumerate(linesTwo):
		if num > 6:
			if line == '':
				break
			listTwo = line.split("', '")
			if id == listTwo[3]:
				if pos == listTwo[1]:
					#print(line)
					matchLine = line
			else:
				continue
		else:
			continue
	return matchLine

for num, line in enumerate(linesOne):
	if num > 6:
		if line == '':
			break
		listOne = line.split("', '")
		if currentId == listOne[3]:
			continue
		else:
			currentId = listOne[3]
		match = searchSub(currentId, listOne[1])
		if match != '':
			print("\nFound match")
			print("\nin run_ref_3")   
			print(line)
			print("\nin run_ref_4")
			print(match)
		else:
			continue
			#print("\nCould not find match")
print("Procedure complete")
