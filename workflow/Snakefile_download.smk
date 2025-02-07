# Workflow RNA-seq processing

import pandas as pd

configfile: 'config/config.yaml'

Samples = {}
for dataset in config['Datasets']:
    Samples[dataset] = pd.read_csv("resources/"+dataset+"/SRR_Acc_List.txt", sep="\t", header=None)

rule all:
    input:
        expand("results/{dataset}/fastq/{sample}.out", dataset=config['Dataset'], sample=Samples[config['Dataset'][0]])
        #expand("resources/{dataset}/SRR_per_SampleName.txt", dataset=config['Datasets'])


##-------------------------##
##   Download fastq data   ##
##-------------------------##

rule download_fastq:
    input:
        acc_list = "resources/{dataset}/SRR_Acc_List.txt"
    output:
        stdout="log/{dataset}/fasterq-dump/{sample}.out"
    params:
        outfold="resources/{dataset}/fastq",
        stderr="log/{dataset}/fasterq-dump/{sample}.err"
    threads: 12
    shell:
        """
        fasterq-dump {wildcards.sample} --outdir {params.outfold} --temp tmp --threads {threads} --mem 10000MB 1> {output.stdout} 2> {params.stderr}
        """

##-------------------------##
##   Get SRR per sample    ##
##-------------------------##

rule srr_sample_table:
    input:
        srarun="resources/{dataset}/SraRunTable.txt",
        gsm_sample="resources/{dataset}/GSMID_SampleName.txt"
    output:
        "resources/{dataset}/SRR_per_SampleName.txt"
    shell:
        """
        python scripts/get_SRR_SampleName.py --SraRunTable {input.srarun} --GSMID_SampleName {input.gsm_sample} --outfile {output}
        """