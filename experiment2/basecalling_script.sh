#!/bin/sh

INPUT=/export/valenfs/data/raw_data/minion/20201016_DNA_Inosine-Pool1/20201016_1505_MN29576_FAO10324_b606dd53/fast5/
OUTPUT=/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/1_basecall/output/

nice -n 10 ~/ont-guppy-cpu/bin/./guppy_basecaller \
	--input_path $INPUT \
	--recursive \
	--save_path $OUTPUT \
	--config dna_r9.4.1_450bps_hac.cfg \
	--fast5_out \
	--trim_strategy none \
	--num_callers 11 \
	--cpu_threads_per_caller 6 \
	2>&1 | tee logfile.txt		
	
