#!/bin/sh

nanocompore sampcomp \
	--file_list1 /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/13_invariate_ref_event_collapsed/construct1/A_in_reference/out_eventalign_collapse.tsv \
	--file_list2 /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/13_invariate_ref_event_collapsed/construct2/out_eventalign_collapse.tsv \
	--label1 construct1_w_A_ref \
	--label2 construct2 \
	--fasta /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/14_invariate_ref_comp/A_vs_I/seg_ref_w_A.fasta	\
	--outpath /export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/19_6mer_comp_reprise/A_vs_I/results/ \
	--min_coverage 1 \
        --min_ref_length 50 \
	--max_invalid_kmers_freq 1 \
	--overwrite \
	--nthreads 100 \
    	--pvalue_thr 0.01 \
    	--comparison_methods GMM,KS,TT,MW \
    	--logit
