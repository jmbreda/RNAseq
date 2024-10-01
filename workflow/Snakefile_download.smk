# Workflow RNA-seq processing

import pandas as pd

configfile: 'config/config.yaml'

Samples = {}
for dataset in config['Datasets']:
    Samples[dataset] = pd.read_csv(f'resources/{dataset}/SRR_Acc_List.txt', sep="\t", header=None)

rule all:
    input:
        #expand("results/{dataset}/fastq/{sample}.out", dataset=config['Dataset'], sample=Samples[config['Dataset'][0]])
        expand("resources/{dataset}/SRR_per_SampleName.txt", dataset=config['Datasets'])

rule download_fastq:
    output:
        stdout="results/{dataset}/fastq/{sample}.out",
        stderr="results/{dataset}/fastq/{sample}.err"
    params:
        outfold="results/A_circadian_gene_expression_atlas_in_mouse/fastq"
    threads: 12
    shell:
        """
        fasterq-dump {wildcards.sample} --outdir {params.outfold} --temp tmp --threads {threads} --mem 10000MB 1> {output.stdout} 2> {output.stderr}
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
        python script/get_SRR_SampleName.py --SraRunTable {input.srarun} --GSMID_SampleName {input.gsm_sample} --output {output}
        """