#!/bin/bash

dataset=$1

infile=resources/$dataset/SRR_per_SampleName.txt
infolder=resources/$dataset/fastq_SRR
outfolder=resources/$dataset/fastq

# get number of lines in file
N=$(wc -l < $infile)
for i in $(seq 2 $N)
do

    line=$(sed -n "${i}p" $infile)
    SampleName=$(echo $line | cut -d' ' -f 1)
    SRR=$(echo $line | cut -d ' ' -f 2)
    fq1="$infolder"/"${SampleName}"_R1.fastq.gz
    fq2="$infolder"/"${SampleName}"_R2.fastq.gz

    echo "SampleName: $SampleName"
    echo "SRR: $SRR"
    
    python scripts/concatenate.py --SampleName "$SampleName" --srr_list "$SRR" --infolder "$infolder" --fq1 "$fq1" --fq2 "$fq2"
    
done