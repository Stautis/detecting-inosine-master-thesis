#!/bin/sh

java -jar ~/picard/build/libs/picard.jar CreateSequenceDictionary \
	-R /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/1_alignment/construct5/construct5.fasta \
	-O /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/1_alignment/construct5/construct5.dict 
