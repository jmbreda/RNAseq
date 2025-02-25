#!/bin/bash

DATASETS=(anderson_2023 mekbib_2022 astafev_2024 Weger_CellMetab_2019 Koike_Science_2012  Zhang_PNAS_2014)

for dataset in "${DATASETS[@]}"; do
    echo "Dataset: $dataset"
    
    srarun="../resources/${dataset}/SraRunTable.txt"
    gsm_sample="../resources/${dataset}/GSMID_SampleName.txt"
    output="../resources/${dataset}/SRR_per_SampleName.txt"

    if [ ! -e $output ] && [ -e $srarun ] && [ -e $gsm_sample ]; then
        echo $srarun
        echo $gsm_sample
        echo $output
        python get_SRR_SampleName.py --SraRunTable ${srarun} --GSMID_SampleName ${gsm_sample} --outfile ${output}
    fi
done
