#!/bin/sh

REF_PATH=/Home/ii/thomasks/Ref/alkbh3_inosine_hompolymer_context.fasta
CONCAT_MERGE_FILE=/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/alkbh3_analysis/run_ref_4/aln_sort.fastq
BAM=/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/alkbh3_analysis/run_ref_4/aln_sort.bam
OUT=/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/alkbh3_analysis/run_ref_4/alkbh3_eventalign.tsv

nice -n 3 ~/nanopolish/./nanopolish eventalign \
	-r $CONCAT_MERGE_FILE \
	-b $BAM \
	-g $REF_PATH \
	--progress \
        --threads 100 \
        --print-read-names \
	--scale-events \
	--samples > $OUT	
