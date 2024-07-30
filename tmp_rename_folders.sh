#!/bin/bash

srr_to_samplename='resources/Zhang_PNAS_2014/SRR_per_SampleName.txt'
folder='results/Zhang_PNAS_2014/star'

N=$(wc -l < $srr_to_samplename)

for i in $(seq 1 $N)
do
    srr=$(awk -v i=$i 'NR==i {print $2}' $srr_to_samplename)
    sample_name=$(awk -v i=$i 'NR==i {print $1}' $srr_to_samplename)
    echo $srr " -> " $sample_name
    mv $folder/$srr $folder/$sample_name
done
