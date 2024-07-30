#!/bin/bash

# Download data from circadian gene expression atlas in mouse
SRA_accession_list=$(cat resources/A_circadian_gene_expression_atlas_in_mouse/SRR_Acc_List.txt)

for sra in ${SRA_accession_list[@]}
do
	if ! [ -e results/A_circadian_gene_expression_atlas_in_mouse/fastq/${sra}_1.fastq ] && ! [ -e results/A_circadian_gene_expression_atlas_in_mouse/fastq/${sra}_2.fastq ]
	then
		echo "Downloading $sra"
		fasterq-dump $sra --outdir results/A_circadian_gene_expression_atlas_in_mouse/fastq --threads 12 --mem 4000MB
	fi
done
