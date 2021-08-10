import csv

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

aoi, sum = restrict(1248, 1254)
count=0
print("\nWriting area of interest\n")
for row in aoi:
	print(row)
	count+=1
print("\nThere are", count,"records detected in the currently defined area of interest\n")
