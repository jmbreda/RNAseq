#!/bin/bash

dataset="Atger_PNAS_2015"
infile="resources/${dataset}/SRR_Acc_List.txt"
outfold="resources/${dataset}/fastq"
mapfile -t SRR < $infile

for srr in "${SRR[@]}"
do
    if [ ! -f "${outfold}/${srr}_1.fastq" ] || [ ! -f "${outfold}/${srr}_2.fastq" ]; then
        echo "Downloading $srr"
        stdout="resources/${dataset}/log/${srr}.out"
        stderr="resources/${dataset}/log/${srr}.err"
        fasterq-dump "$srr" --outdir $outfold --temp tmp --threads 6 --mem 10000MB 1> "$stdout" 2> "$stderr"
    fi
done