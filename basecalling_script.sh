#!/bin/sh

INPUT=/export/valenfs/data/raw_data/minion/2020_08_19_DNA_Inosine2/20200818_1417_MN29576_FAM96106_e8a61b8a/fast5
OUTPUT=/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/2_basecalled/output

nice -n 10 ~/guppy/guppy_341/ont-guppy-cpu/bin/./guppy_basecaller \
	--input_path $INPUT \
	--recursive \
	--save_path $OUTPUT \
	--config dna_r9.4.1_450bps_hac.cfg \
	--fast5_out \
	--trim_strategy none \
	--num_callers 11 \
	--cpu_threads_per_caller 6 \
	2>&1 | tee logfile.txt		
	
