#!/bin/sh

REF_PATH=/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/2_aligned/construct1/A_in_reference/construct1.fasta
CONCAT_MERGE_FILE=/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/1_basecall/merged.fastq
OUT=/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/2_aligned/construct1/A_in_reference/output_alignment

SAM_OUT=$OUT/aln.sam
LOG_OUT=$OUT/aln.log

#--secondary=no
nice -n 3 minimap2 \
	-t 40 \
	-ax map-ont \
	-y \
	$REF_PATH $CONCAT_MERGE_FILE > $SAM_OUT \
	2> $LOG_OUT

