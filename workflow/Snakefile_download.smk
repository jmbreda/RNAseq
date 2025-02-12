# Workflow RNA-seq processing

import pandas as pd

configfile: 'config/Koike_Science_2012.yml'
dataset = config['Name']
Samples = pd.read_csv("resources/"+dataset+"/SRR_Acc_List.txt", sep="\t", header=None)[0]

rule all:
    input:
        expand("log/"+dataset+"/fasterq-dump/{sample}.done", sample=Samples),
        #"resources/"+dataset+"/SRR_per_SampleName.txt"


##-------------------------##
##   Download fastq data   ##
##-------------------------##

rule download_fastq:
    input:
        acc_list = "resources/"+dataset+"/SRR_Acc_List.txt"
    output:
        done="log/"+dataset+"/fasterq-dump/{sample}.done"
    params:
        outfold="resources/"+dataset+"/fastq",
        mem="72000MB"
    threads: 12
    shell:
        """
        fasterq-dump {wildcards.sample} --outdir {params.outfold} --temp tmp --threads {threads} --mem {params.mem}
        touch {output.done}
        """

##-------------------------##
##   Get SRR per sample    ##
##-------------------------##

rule srr_sample_table:
    input:
        srarun="resources/"+dataset+"/SraRunTable.txt",
        gsm_sample="resources/"+dataset+"/GSMID_SampleName.txt"
    output:
        "resources/"+dataset+"/SRR_per_SampleName.txt"
    shell:
        """
        python scripts/get_SRR_SampleName.py --SraRunTable {input.srarun} --GSMID_SampleName {input.gsm_sample} --outfile {output}
        """