#!/bin/sh

REF_PATH=/Home/ii/thomasks/Ref/alkbh3_inosine_hompolymer_context.fasta
CONCAT_MERGE_FILE=/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/2_basecalled/output/concat.fastq
OUT=/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/2_aligned/rev_output

SAM_OUT=$OUT/aln.sam
LOG_OUT=$OUT/aln.log

nice -n 3 minimap2 \
	-t 40 \
	-ax map-ont \
	-y \
	--rev-only \
	--secondary=no \
	$REF_PATH $CONCAT_MERGE_FILE > $SAM_OUT \
	2> $LOG_OUT

