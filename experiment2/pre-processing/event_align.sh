#!/bin/sh

REF_PATH=/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/2_aligned/construct1/A_in_reference/construct1.fasta
CONCAT_MERGE_FILE=/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/1_basecall/merged.fastq
BAM=/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/2_aligned/construct1/A_in_reference/output_alignment/for_map_sort.bam
OUT=/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/3_event_align/construct1/A_in_reference/A_ref_eventalign.tsv

nice -n 3 nanopolish eventalign \
	-r $CONCAT_MERGE_FILE \
	-b $BAM \
	-g $REF_PATH \
	--progress \
        --threads 100 \
        --print-read-names \
	--scale-events \
	--samples > $OUT	
