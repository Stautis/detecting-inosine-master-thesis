#!/bin/sh

java -jar ~/jvarkit/dist/sam2tsv.jar \
				-R /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/1_alignment/construct5/construct5.fasta \
				-o /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/2_tsv_conversion/construct5/construct5_test.tsv \
				/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/1_alignment/construct5/for_map_sort.bam 
