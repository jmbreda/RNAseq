# RNAseq processing pipeline

# Prior to running the pipeline, the following files need to be prepared:
# SraRunTable.txt
# SRR_Acc_List.txt
# GSMID_SampleName.txt
# -> SRR_per_SampleName.txt generated with workflow/Snakefile_download.smk

# slurm submission: 
# snakemake -s workflow/Snakefile -j 999 --cluster-config config/cluster.json --cluster "sbatch --job-name {cluster.name} --output {cluster.stdout} --error {cluster.stderr} --qos {cluster.qos} --time {cluster.time} --mem {cluster.mem} --nodes {cluster.nodes} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus-per-task}"

import os

configfile: "config/Koike_Science_2012.yml"

def get_SampleNames(dataset=config['Name']):
    infile = "resources/"+dataset+"/SRR_per_SampleName.txt"
    with open(infile, 'r') as f:
        sample_names = [l.split('\t')[0] for l in f.read().splitlines()[1:]]
    
    return sample_names

# constants
wildcard_constraints:
    spec = config['Species'],
    sample = '|'.join(get_SampleNames())+'|test',

rule all:
    input:
        # get genome annotation
        #expand("resources/genome/{spec}/{file}",spec=config['Species'],file=config['GenomeFiles']),
        #expand("resources/genome/{spec}/ensembl.gtf.{splicing}.final.bed",spec=config['Species'],splicing=['introns','exons','genes']),
        #expand("resources/genome/{spec}/premrna.mrna.{transcriptome}.idx",spec=config['Species'],transcriptome=config['Transcriptome']),
        #
        # map to genome with STAR
        #"results/{dataset}/star/{sample}/Aligned.sortedByCoord.out.bam".format(dataset=config['Name'],sample=get_SampleNames()),
        #"results/Koike_Science_2012/star/liver_CT0/Aligned.sortedByCoord.out.bam",
        #expand("results/{dataset}/star/{sample}/Aligned.sortedByCoord.out.bam.bai",sample=get_SampleNames(),dataset=config['Name']),
        #
        # map to genome with bwa
        expand("results/{dataset}/bwa/{sample}_{strand}_coverage.bedgraph",dataset=config['Name'],sample=get_SampleNames(),strand=config['Strand']),
        #
        # map to transcriptome with kallisto
        #expand('results/{dataset}/kallisto/mrna_tpm_table.tab',dataset=config['Name']),

##---------------------------------------------------------------------------##
##   Download gtf, cDNA (coding transcripts), and DNA from Gencode/Ensembl   ##
##---------------------------------------------------------------------------##

rule download_gtf:
    output: 
        gtf = "resources/genome/{spec}/gene_annotation.gtf"
    params: 
        url = lambda wildcards: config['URLs']['gtf']
    shell:
        """
        wget -O {output.gtf}.gz {params.url}
        gunzip {output.gtf}.gz
        """

rule download_cdna:
    output: 
        cdna = "resources/genome/{spec}/cdna.fa"
    params: 
        url = lambda wildcards: config['URLs']['cDNA']
    shell:
        """
        wget -O {output.cdna}.gz {params.url}
        gunzip {output.cdna}.gz
        """

rule download_genome:
    output: 
        genome = "resources/genome/{spec}/genome.fa"
    params:
        url = lambda wildcards: config['URLs']['genome']
    shell:
        """
        wget -O {output.genome}.gz {params.url}
        gunzip {output.genome}.gz
        """

rule chrom_size:
    input:
        genome = "resources/genome/{spec}/genome.fa"
    output:
        chrom_size = "resources/genome/{spec}/chrom.sizes"
    shell:
        """
        if [ ! $(hostname -s) == "haas005" ]; then
            ml gcc/11.3.0 samtools
        fi
        samtools faidx {input.genome}
        cut -f1,2 {input.genome}.fai > {output.chrom_size}
        """

##-------------------##
## CSfastq to fastq  ##
##-------------------##

rule csfastq_to_fastq:
    input:
        csfastq = "resources/{dataset}/fastq/{sample}.csfastq.gz"
    output:
        fastq = "resources/{dataset}/fastq/{sample}.fastq.gz"
    shell:
        """
        scripts/csfq2fq.pl <( zcat {input.csfastq} ) | gzip > {output.fastq}
        """

##------------------------##
##  STAR : map to genome  ##
##------------------------##

rule star_index:
    input:
        genome = "resources/genome/{spec}/genome.fa",
        gtf = "resources/genome/{spec}/gene_annotation.gtf"
    output:
        genome_dir = directory("resources/genome/{spec}/star_index_{read_length}/")
    params:
        overhang = config['ReadLength'] - 1,
        mem = 160000000000
    threads: 24
    shell:
        """
        module load gcc/11.3.0 star
        mkdir -p {output.genome_dir}
        STAR --runMode genomeGenerate \
             --genomeDir {output.genome_dir} \
             --genomeFastaFiles {input.genome} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.overhang} \
             --limitGenomeGenerateRAM {params.mem} \
             --runThreadN {threads}
        """

rule star_align:
    input:
        fq='resources/{dataset}/fastq/{sample}.fastq.gz',
        gtf = "resources/genome/{spec}/gene_annotation.gtf".format(spec=config['Species']),
        genome_dir = "resources/genome/{spec}/star_index_{read_length}".format(spec=config['Species'],read_length=config['ReadLength'])
    output:
        bam="results/{dataset}/star/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        prefix="results/{dataset}/star/{sample}/",
        outsamtype = "BAM SortedByCoordinate",
        mem = 160000000000,
    threads: 24
    shell:
        """
        module load gcc/11.3.0 star
        if [ -e {params.prefix}tmp/ ]; then rm -r {params.prefix}tmp/; fi
        STAR --runMode alignReads \
             --outFileNamePrefix {params.prefix} \
             --outTmpDir {params.prefix}tmp/ \
             --runDirPerm All_RWX \
             --outSAMtype {params.outsamtype} \
             --readFilesIn {input.fq} \
             --readFilesCommand zcat \
             --genomeDir {input.genome_dir} \
             --sjdbGTFfile {input.gtf} \
             --limitBAMsortRAM {params.mem} \
             --runThreadN {threads}
        """

rule sam_index:
    input:
        bam="results/{dataset}/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        bam_index="results/{dataset}/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        ml gcc/11.3.0 samtools
        samtools index {input}
        """

rule coverage_bedgraph:
    input:
        bam="results/{dataset}/star/{sample}/Aligned.sortedByCoord.out.bam",
        bam_index="results/{dataset}/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        bg="results/{dataset}/star/{sample}/coverage.bedgraph"
    params:
        strand = lambda wildcards: '-' if wildcards.strand == "forward" else '+' if wildcards.strand == "reverse" else '',
        flag_pe = lambda wildcards: '81' if wildcards.strand == "forward" else '97' if wildcards.strand == "reverse" else '',
        flag_se = lambda wildcards: '16' if wildcards.strand == "forward" else '0' if wildcards.strand == "reverse" else ''
    shell:
        """
        if [ ! $(hostname -s) == "haas005" ]; then
            ml gcc/11.3.0 samtools bedtools2
        fi
        samtools view -h {input.bam} | awk '$1 ~ /^@/ || $2 == "{params.flag_pe}" || $2 == "{params.flag_se}"' | samtools view -b | bedtools genomecov -ibam stdin -bg -strand {params.strand} -5 | grep "^\<chr" > {output.bg}
        """

##----------------------##
## BWA : map to Genome  ##
##----------------------##

rule bwa_index:
    input:
        genome="resources/genome/{spec}/genome.fa"
    output:
        bwa_index="resources/genome/{spec}/genome.fa.bwt"
    shell:
        """
        bwa index {input.genome}
        """

rule bwa_aln_se:
    input:
        genome="resources/genome/{spec}/genome.fa".format(spec=config['Species']),
        genome_index="resources/genome/{spec}/genome.fa.bwt".format(spec=config['Species']),
        fq="resources/{dataset}/fastq/{sample}.fastq.gz",
    output:
        sai="results/{dataset}/bwa/{sample}.sai"
    threads: 12
    shell:
        """
        bwa aln -t {threads} -f {output.sai} {input.genome} {input.fq}
        """

rule bwa_samse:
    input:
        genome = "resources/genome/{spec}/genome.fa".format(spec=config['Species']),
        gemone_index = "resources/genome/{spec}/genome.fa.bwt".format(spec=config['Species']),
        sai = "results/{dataset}/bwa/{sample}.sai",
        fq = "resources/{dataset}/fastq/{sample}.fastq.gz",
    output:
        sam="results/{dataset}/bwa/{sample}.sam"
    threads: 3
    shell:
        """
        bwa samse -f {output.sam} {input.genome} {input.sai} {input.fq}
        """

rule mapped_uniq_MAPQ_sorted_se:
    input:
        sam="results/{dataset}/bwa/{sample}.sam"
    output:
        bam="results/{dataset}/bwa/{sample}.bam"
    params:
        MAPQ_min=config['bwa']['MAPQ_min'],
        mem=config['bwa']['mem']
    threads: 12
    shell:
        """
        if [ ! $(hostname -s) == "haas005" ]; then
            ml gcc/11.3.0 samtools
        fi
        samtools view -@ {threads} -h -q {params.MAPQ_min} -F 4 {input.sam} |\
            awk '$1 ~ /^@/ || $2 == "16" || $2 == "0"' |\
            grep -v 'XA:Z:'| grep -v 'SA:Z:' |\
            samtools sort -@ {threads} -m {params.mem} > {output.bam}
        """

rule bwa_bam_index:
    input:
        bam="results/{dataset}/bwa/{sample}.bam"
    output:
        bam_index="results/{dataset}/bwa/{sample}.bam.bai"
    shell:
        """
        ml gcc/11.3.0 samtools
        samtools index {input.bam}
        """

rule bwa_coverage_bedgraph:
    input:
        bam="results/{dataset}/bwa/{sample}.bam",
        bam_index="results/{dataset}/bwa/{sample}.bam.bai"
    output:
        bg="results/{dataset}/bwa/{sample}_{strand}_coverage.bedgraph"
    params:
        strand = lambda wildcards: '-' if wildcards.strand == "forward" else '+' if wildcards.strand == "reverse" else '',
        flag_pe = lambda wildcards: '81' if wildcards.strand == "forward" else '97' if wildcards.strand == "reverse" else '',
        flag_se = lambda wildcards: '16' if wildcards.strand == "forward" else '0' if wildcards.strand == "reverse" else ''
    shell:
        """
        if [ ! $(hostname -s) == "haas005" ]; then
            ml gcc/11.3.0 samtools bedtools2
        fi
        samtools view -h {input.bam} | awk '$1 ~ /^@/ || $2 == "{params.flag_pe}" || $2 == "{params.flag_se}"' | samtools view -b | bedtools genomecov -ibam stdin -bg -strand {params.strand} -5 | grep "^\<chr" > {output.bg}
        """

##----------------------------------##
##  Kallisto: map to transcriptome  ##
##----------------------------------##

rule kallisto_index:
    input:
        cdna = "resources/genome/{spec}/cdna.fa".format(spec=config['Species']),
    output: 
        index = "resources/genome/{spec}/cdna.idx".format(spec=config['Species']),
    threads: 24
    shell:
        """
        ml gcc/11.3.0 kallisto
        kallisto index -i {output.index} {input.cdna}
        """

rule kallisto_quant:
    input:
        index="resources/genome/{spec}/cdna.idx".format(spec=config['Species']),
        fq="resources/{dataset}/fastq/{sample}.fastq.gz"
    output:
        counts="results/{dataset}/kallisto/{sample}/abundance.tsv"
    params:
        dir=directory("results/{dataset}/kallisto/{sample}"),
        fragment_length_mean=config['ReadLength'],
        fragment_length_sd=config['kallisto']['fragment_length_sd']
    threads: 24
    shell:
        """
        module load gcc/11.3.0 kallisto
        kallisto quant \
                 -i {input.index} \
                 --single \
                 --fragment-length {params.fragment_length_mean} \
                 --sd {params.fragment_length_sd} \
                 --pseudobam \
                 --plaintext \
                 -b 50 \
                 -o {params.dir} \
                 -t {threads} \
                 {input.fq}
        """

##-----------------------------------------------------##
##  Make TPM and promoter tables from kallisto output  ##
##-----------------------------------------------------##

rule make_tpm_tables:
    input:
        expand('results/{dataset}/kallisto/{sample}/abundance.tsv', sample=get_SampleNames(), dataset=config['Name'])
    output:
        mrna='results/{dataset}/kallisto/mrna_tpm_table.tab'
    shell:
        """
        python scripts/make_tpm_table.py --intables {input} --out_mrna {output.mrna}
        """
#
#rule countreads:
#    input:
#        bam = "results/star/{sample}/Aligned.sortedByCoord.out.bam",
#        bam_index = "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
#        intron = "resources/genome/{spec}/ensembl.gtf.introns.final.bed".format(spec=config['Species']),
#        exon = "resources/genome/{spec}/ensembl.gtf.exons.final.bed".format(spec=config['Species'])
#    output:
#        "results/counting/{sample}/ie_ol.tsv"
#    shell:
#        """
#        ml gcc/11.3.0 samtools
#        perl scripts/Counting_IE_PairedEnd.pl {input.intron} {input.exon} {input.bam} > {output}
#        """
#
#rule merge_gene_counts:
#    input:
#        expand('results/star/{sample}_ReadsPerGene.out.tab', sample=get_acc_list(config['Name']))
#    output:
#        'results/star/GeneCount_table.tab'
#    shell:
#        """
#        outtable={output}
#        cut -f1 {input[0]} | tail -n+5 > $outtable
#        for i in {input}; do
#            paste $outtable <(cut -f2 $i | tail -n+5) > tmp.txt
#            mv tmp.txt $outtable
#        done
#        """
#
#rule merged_transcriptome:
#    input:
#        expand('results/star/{sample}_Aligned.toTranscriptome.out.bam', sample=get_acc_list(config['Name']))
#    output:
#        'results/star/merged_transcriptome.bam'
#    shell:
#        """
#        samtools sort -@ 12 -o {output} {input}
#        samtools index {output}
#        """