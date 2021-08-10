import sys
import csv
import itertools

def restrict(low_bound, upper_bound):
	print("\nRestricting area")
	with open("run_ref_4/alkbh3_eventalign.tsv") as file:
		read = csv.reader(file, delimiter="\t")
		rows = []
		inosine_signal_sum = 0
		for num, row in enumerate(read):
			if num == 0:
				rows.append(row)
				continue
			if int(row[1]) > low_bound and int(row[1]) < upper_bound:
				rows.append(row)
				inosine_signal_sum+=float(row[6])
			#if num == 100000:
				#break
		print("The sum of event mean signals for the homopolymer between",low_bound,"and",upper_bound," is:", inosine_signal_sum)	
		return rows, inosine_signal_sum

def process(aoi):
	print("\nCalculating signal mean and dwell time per event")
	current_id = aoi[1][3]
	event_signal_sum = 0
	event_st_dev_sum = 0
	event_dwell_time = 0
	count=0
	rows = []
	for num, row in enumerate(aoi):
		if num == 0:
			rows.append(row)
			continue
		if row[3] == current_id:
			count+=1
			event_signal_sum+=float(row[6])
			event_st_dev_sum+=float(row[7])
			event_dwell_time+=float(row[8])
		else: 
			event = [current_id, event_signal_sum/count, event_st_dev_sum/count, event_dwell_time]
			rows.append(event)
			count=1
			event_signal_sum = float(row[6])
			event_st_dev_sum = float(row[7])
			event_dwell_time = float(row[8])
			current_id = row[3]
	return rows
			

print("Starting procedure")

#Calculate for Inosine
aoi, sum = restrict(1248, 1254)
print("\nAverage signal is: ", sum/len(aoi),"\n")
processed_rows = process(aoi)
print("\nProcessing complete")
per_event_signal_mean = 0
per_event_st_dev_mean = 0
per_event_length_mean=0
for num, row in enumerate(processed_rows):
	if num == 0:
		continue
	per_event_signal_mean+=float(row[1])
	per_event_st_dev_mean+=float(row[2])
	per_event_length_mean+=float(row[3])
print("\nFor recorded forward inosine events the following statistics are calculated:")
print("\nAverage event signal mean:", per_event_signal_mean/len(processed_rows))
print("\nAverage standard deviation mean:", per_event_st_dev_mean/len(processed_rows))
print("\nAverage event length:", per_event_length_mean/len(processed_rows)) 	

#Calculate for Thymine
aoi_t, sum = restrict(1273, 1278)
print("\nAverage signal is: ", sum/len(aoi_t),"\n")
processed_rows_t = process(aoi_t)
per_event_signal_mean = 0
per_event_st_dev_mean = 0
per_event_length_mean=0
for num, row in enumerate(processed_rows_t):
        if num == 0:
                continue
        per_event_signal_mean+=float(row[1])
        per_event_st_dev_mean+=float(row[2])
        per_event_length_mean+=float(row[3])
print("\nFor recorded forward thymine events the following statistics are calculated:")
print("\nAverage event signal mean:", per_event_signal_mean/len(processed_rows_t))
print("\nAverage standard deviation mean:", per_event_st_dev_mean/len(processed_rows_t))
print("\nAverage event length:", per_event_length_mean/len(processed_rows_t))

#Calculate for Adenine
aoi_a, sum = restrict(1291, 1296)
print("\nAverage signal is: ", sum/len(aoi_a),"\n")
processed_rows_a = process(aoi_a)
per_event_signal_mean = 0
per_event_st_dev_mean = 0
per_event_length_mean=0
for num, row in enumerate(processed_rows_a):
        if num == 0:
                continue
        per_event_signal_mean+=float(row[1])
        per_event_st_dev_mean+=float(row[2])
        per_event_length_mean+=float(row[3])
print("\nFor recorded forward adenine events the following statistics are calculated:")
print("\nAverage event signal mean:", per_event_signal_mean/len(processed_rows_a))
print("\nAverage standard deviation mean:", per_event_st_dev_mean/len(processed_rows_a))
print("\nAverage event length:", per_event_length_mean/len(processed_rows_a))

#Calculate for Guanine
aoi_g, sum = restrict(1310, 1315)
print("\nAverage signal is: ", sum/len(aoi_g),"\n")
processed_rows_g = process(aoi_g)
per_event_signal_mean = 0
per_event_st_dev_mean = 0
per_event_length_mean=0
for num, row in enumerate(processed_rows_g):
        if num == 0:
                continue
        per_event_signal_mean+=float(row[1])
        per_event_st_dev_mean+=float(row[2])
        per_event_length_mean+=float(row[3])
print("\nFor recorded forward guanine events the following statistics are calculated:")
print("\nAverage event signal mean:", per_event_signal_mean/len(processed_rows_g))
print("\nAverage standard deviation mean:", per_event_st_dev_mean/len(processed_rows_g))
print("\nAverage event length:", per_event_length_mean/len(processed_rows_g))

#Calculate for Cytosine
aoi_c, sum = restrict(1328, 1332)
print("\nAverage signal is: ", sum/len(aoi_c),"\n")
processed_rows_c = process(aoi_c)
per_event_signal_mean = 0
per_event_st_dev_mean = 0
per_event_length_mean=0
for num, row in enumerate(processed_rows_c):
        if num == 0:
                continue
        per_event_signal_mean+=float(row[1])
        per_event_st_dev_mean+=float(row[2])
        per_event_length_mean+=float(row[3])
print("\nFor recorded forward adenine events the following statistics are calculated:")
print("\nAverage event signal mean:", per_event_signal_mean/len(processed_rows_c))
print("\nAverage standard deviation mean:", per_event_st_dev_mean/len(processed_rows_c))
print("\nAverage event length:", per_event_length_mean/len(processed_rows_c))

print("Procedure completed for all areas of interest")
