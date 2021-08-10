#!/bin/sh

fast5_subset \
	-i /export/valenfs/data/raw_data/minion/20201016_DNA_Inosine-Pool1/20201016_1505_MN29576_FAO10324_b606dd53/fast5 \
	-s /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/3_tsv_analysis/construct4/another_run/thirdPosReads \
	-l /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/third_analysis/3_tsv_analysis/construct4/thirdPerfCytosineReads.tsv  \
	-t 128
