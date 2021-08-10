#!/bin/sh

BASECALLED=/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/1_basecalled/output/concat.fasta
INPUT=/export/valenfs/data/raw_data/minion/20181012_max_DNA_Inosine/20181012_1311_1/fast5/
SUMMARY=/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/1_basecalled/output/sequencing_summary.txt

nice -n 3 ~/nanopolish/./nanopolish index \
	-d $INPUT \
	--sequencing-summary $SUMMARY \
	$BASECALLED
	
