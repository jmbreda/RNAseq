# RNAseq processing pipeline

# Prior to running the pipeline, the following files need to be prepared:
# SraRunTable.txt
# SRR_Acc_List.txt
# GSMID_SampleName.txt
# -> SRR_per_SampleName.txt generated with workflow/Snakefile_download.smk

# slurm submission: 
# snakemake -s workflow/Snakefile -j 999 --cluster-config config/cluster.json --cluster "sbatch --job-name {cluster.name} --output {cluster.stdout} --error {cluster.stderr} --qos {cluster.qos} --time {cluster.time} --mem {cluster.mem} --nodes {cluster.nodes} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus-per-task}"

import os
import pandas as pd

configfile: "config/Koike_Science_2012.yml"


def get_SRR(dataset=config['Name']):
    infile = "resources/"+dataset+"/SRR_per_SampleName.txt"
    with open(infile, 'r') as f:
        srr_ids = [l.split('\t')[1] for l in f.read().splitlines()[1:]]
    
    return srr_ids

def get_SampleNames(dataset=config['Name']):
    infile = "resources/"+dataset+"/SRR_per_SampleName.txt"
    with open(infile, 'r') as f:
        sample_names = [l.split('\t')[0] for l in f.read().splitlines()[1:]]
    
    return sample_names

def get_srr(wildcards):
    infile = "resources/"+wildcards.dataset+"/SRR_per_SampleName.txt"
    df = pd.read_csv(infile, sep='\t',index_col=0)
    srr = df.at[wildcards.sample,'Run']
    
    return srr

# constants
wildcard_constraints:
    spec = config['Species'],
    srr = '|'.join(get_SRR(dataset=config['Name'])),
    sample = '|'.join(get_SampleNames(dataset=config['Name']))+'|test',
    strand = '|'.join(config['Strand']),
    chr = '|'.join(config['Chromosomes']),
    bin_size = '|'.join([str(b) for b in config['BinSize'] ]),
    dataset = config['Name'],
    genome_aligner = '|'.join(config['GenomeAligner']),

rule all:
    input:
        # get genome annotation
        #expand("resources/genome/{spec}/{file}",spec=config['Species'],file=config['GenomeFiles']),
        #expand("resources/genome/{spec}/ensembl.gtf.{splicing}.final.bed",spec=config['Species'],splicing=['introns','exons','genes']),
        #expand("resources/genome/{spec}/premrna.mrna.{transcriptome}.idx",spec=config['Species'],transcriptome=config['Transcriptome']),
        #
        # translate csfastq to fastq
        #expand("resources/{dataset}/fastq/{sample}.fastq.gz",dataset=config['Name'],sample=get_SampleNames(config['Name'])),
        #
        # map to genome with bwa
        #expand("results/{dataset}/bwa/coverage/{sample}_{strand}.bedgraph",dataset=config['Name'],sample=get_SampleNames(),strand=config['Strand']),
        #
        # map to transcriptome with kallisto
        #expand('results/{dataset}/kallisto/mrna_tpm_table.tab',dataset=config['Name']),
        #
        # map to genome with STAR
        #expand("results/{dataset}/star/{sample}/Aligned.sortedByCoord.out.bam",dataset=config['Name'],sample=get_SampleNames()),
        #"results/Koike_Science_2012/star/CT00/Aligned.sortedByCoord.out.bam",
        #expand("results/{dataset}/star/mapping/{sample}/Aligned.out.sam",dataset=config['Name'],sample=get_SampleNames()),
        #expand("results/{dataset}/star/{sample}/Aligned.sortedByCoord.out.bam.bai",sample=get_SampleNames(),dataset=config['Name']),
        #expand("results/{dataset}/star/coverage/{sample}_{strand}.bedgraph",dataset=config['Name'],sample=get_SampleNames(),strand=config['Strand']),
        #expand("results/{dataset}/star/coverage/{mapped}/{sample}_{strand}.bg",dataset=config['Name'],mapped=config['star']['Mapped'],sample=get_SampleNames(),strand=config['Strand']),
        #
        # genome coverage 
        #expand("results/{dataset}/{genome_aligner}/norm_coverage/{sample}_{strand}.bw",dataset=config['Name'],sample=get_SampleNames(),strand=config['Strand']),
        #expand("results/{dataset}/{genome_aligner}/binned_norm_coverage/{sample}_{strand}_bin{bin_size}bp.bw",dataset=config['Name'],sample=get_SampleNames(),strand=config['Strand'],bin_size=config['BinSize']),
        #expand("results/{dataset}/{genome_aligner}/log_binned_norm_coverage/{sample}_{strand}_bin{bin_size}bp.bw",dataset=config['Name'],sample=get_SampleNames(),strand=config['Strand'],bin_size=config['BinSize']),
        expand("results/{dataset}/star/binned_norm_coverage/{mapped}/expression_tables/bin_expression_table_bin{bin_size}bp.csv",dataset=config['Name'],mapped=config['star']['Mapped'],bin_size=config['BinSize'])
        #expand("results/{dataset}/star/binned_norm_coverage/Unique_rpm/expression_tables/bin_expression_table_bin10000bp.csv",dataset=config['Name'])


##-------------------------------##
##  Rename and CSfastq to fastq  ##tmux atta
##-------------------------------##

rule rename:
    params:
        srr = lambda wildcards: get_srr(wildcards),
        dataset = config['Name'],
    output:
        csfastq = "resources/{dataset}/fastq/{sample}.csfastq.gz"
    shell:
        """
        infile="resources/{params.dataset}/fastq/{params.srr}.csfastq.gz"
        mv $infile {output.csfastq}
        """

rule csfastq_to_fastq:
    input:
        csfastq = "resources/{dataset}/fastq/{sample}.csfastq.gz"
    output:
        fastq = "resources/{dataset}/fastq/{sample}.fastq.gz"
    shell:
        """
        scripts/csfq2fq.pl <( zcat {input.csfastq} ) | gzip > {output.fastq}
        """

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
        #if [ ! $(hostname -s) == "haas005" ]; then
        #    ml gcc/11.3.0 samtools
        #fi
        samtools faidx {input.genome}
        cut -f1,2 {input.genome}.fai > {output.chrom_size}
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
        sam="results/{dataset}/star/mapping/{sample}/Aligned.out.sam"
    params:
        outfolder = "results/{dataset}/star/mapping/{sample}/",
        prefix = "tmpdir/{dataset}/star/{sample}/"
    threads: 24
    shell:
        """
        module load gcc/11.3.0 star
        if [ -e {params.prefix}tmp/ ]; then rm -r {params.prefix}tmp/; fi
        STAR --runMode alignReads \
             --outFileNamePrefix {params.prefix} \
             --outTmpDir {params.prefix}tmp/ \
             --runDirPerm All_RWX \
             --readFilesIn {input.fq} \
             --readFilesCommand zcat \
             --genomeDir {input.genome_dir} \
             --sjdbGTFfile {input.gtf} \
             --runThreadN {threads}
        if [ ! -d {params.outfolder} ]; then mkdir -p {params.outfolder}; fi            
        mv {params.prefix}* {params.outfolder}
        """

rule sam_to_sorted_bam:
    input:
        sam="results/{dataset}/star/mapping/{sample}/Aligned.out.sam"
    output:
        bam="results/{dataset}/star/mapping/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        mem="6G"
    threads: 24
    shell:
        """
        if [ ! $(hostname -s) == "haas005" ]; then
            ml gcc/11.3.0 samtools
        fi
        samtools view -@ {threads} -h -F 4 {input.sam} |\
            awk '$1 ~ /^@/ || $2 == "16" || $2 == "0"' |\
            samtools sort -@ {threads} -m {params.mem} > {output.bam}
        """

rule sam_index:
    input:
        bam="results/{dataset}/star/mapping/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        bai="results/{dataset}/star/mapping/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        if [ ! $(hostname -s) == "haas005" ]; then
            ml gcc/11.3.0 samtools
        fi
        samtools index {input.bam}
        """

rule rpm_coverage:
    input:
        bam = "results/{dataset}/star/mapping/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/{dataset}/star/mapping/{sample}/Aligned.sortedByCoord.out.bam.bai",
    output:
        bg_u_fwd = "results/{dataset}/star/coverage/Unique_rpm/{sample}_forward.bg",
        bg_u_rev = "results/{dataset}/star/coverage/Unique_rpm/{sample}_reverse.bg",
        bg_um_fwd = "results/{dataset}/star/coverage/UniqueMultiple_rpm/{sample}_forward.bg",
        bg_um_rev = "results/{dataset}/star/coverage/UniqueMultiple_rpm/{sample}_reverse.bg",
    params:
        outfolder = "results/{dataset}/star/coverage/{sample}/",
        bg_u_str1 = "results/{dataset}/star/coverage/{sample}/Signal.Unique.str1.out.bg",
        bg_u_str2 = "results/{dataset}/star/coverage/{sample}/Signal.Unique.str2.out.bg",
        bg_um_str1 = "results/{dataset}/star/coverage/{sample}/Signal.UniqueMultiple.str1.out.bg",
        bg_um_str2 = "results/{dataset}/star/coverage/{sample}/Signal.UniqueMultiple.str2.out.bg"
    threads: 12
    shell:
        """
        module load gcc/11.3.0 star        
        mkdir -p {params.outfolder};
        chmod -R 777 {params.outfolder};
        STAR --runMode inputAlignmentsFromBAM \
             --runThreadN {threads} \
             --inputBAMfile {input.bam} \
             --outWigType bedGraph \
             --outFileNamePrefix {params.outfolder}
        grep "^chr" {params.bg_u_str1}  > {output.bg_u_fwd}
        grep "^chr" {params.bg_u_str2}  > {output.bg_u_rev}
        grep "^chr" {params.bg_um_str1} > {output.bg_um_fwd}
        grep "^chr" {params.bg_um_str2} > {output.bg_um_rev}
        """

#rule coverage_bedgraph:
#    input:
#        bam="results/{dataset}/star/mapping/{sample}/Aligned.sortedByCoord.out.bam",
#        bam_index="results/{dataset}/star/mapping/{sample}/Aligned.sortedByCoord.out.bam.bai",
#        chrom_size="resources/genome/{spec}/chrom.sizes".format(spec=config['Species'])
#    output:
#        bg="results/{dataset}/star/coverage/{sample}_{strand}.bedgraph"
#    params:
#        strand = lambda wildcards: '+' if wildcards.strand == "forward" else '-' if wildcards.strand == "reverse" else '.',
#    shell:
#        """
#        if [ ! $(hostname -s) == "haas005" ]; then
#            ml gcc/11.3.0 samtools bedtools2
#        fi
#        bedtools genomecov -ibam {input.bam} -bg -strand {params.strand} | grep "^\<chr" > {output.bg}
#        """
#
##-----------------##
##  Norm coverage  ##
##-----------------##

#rule get_total_count:
#    input:
#        bg="results/{dataset}/{genome_aligner}/coverage/{sample}_{strand}.bedgraph"
#    output:
#        counts="results/{dataset}/{genome_aligner}/coverage/{sample}_{strand}_total_counts.txt"
#    shell:
#        """
#        awk '{{sum+=$4}} END {{print sum}}' {input.bg} > {output.counts}
#        """
#
#rule normalize_sort_bedgraph:
#    input:
#        bg="results/{dataset}/{genome_aligner}/coverage/{sample}_{strand}.bedgraph",
#        total_count="results/{dataset}/{genome_aligner}/coverage/{sample}_{strand}_total_counts.txt",
#        total_counts=expand("results/{{dataset}}/{{genome_aligner}}/coverage/{sample}_{strand}_total_counts.txt",sample=get_SampleNames(config['Name']), strand=config['Strand'])
#    output:
#        norm_bg="results/{dataset}/{genome_aligner}/norm_coverage/{sample}_{strand}.bedgraph"
#    shell:
#        """
#        ./scripts/normalize_bedgraph.sh {input.bg} {input.total_count} {output.norm_bg} {input.total_counts}
#        """
#
#rule norm_coverage_bw:
#    input:
#        norm_bg="results/{dataset}/{genome_aligner}/norm_coverage/{sample}_{strand}.bedgraph",
#        chrom_size="resources/genome/{spec}/chrom.sizes".format(spec=config['Species'])
#    output:
#        bw="results/{dataset}/{genome_aligner}/norm_coverage/{sample}_{strand}.bw"
#    params:
#        tmp_sorted="results/{dataset}/{genome_aligner}/norm_coverage/{sample}_{strand}.sorted.tmp"
#    shell:
#        """
#        sort -k1,1 -k2,2n {input.norm_bg} > {params.tmp_sorted}
#        bedGraphToBigWig {params.tmp_sorted} {input.chrom_size} {output.bw}
#        rm {params.tmp_sorted}
#        """
#
##--------------------##
##  Bin-Log/coverage  ##
##--------------------##

rule make_bin_bed:
    input:
        chrom_size="resources/genome/{spec}/chrom.sizes".format(spec=config['Species'])
    output:
        bed="results/{dataset}/star/binned_norm_coverage/bin{bin_size}bp.bed"
    shell:
        """
        python scripts/make_bins.py --chrom_size {input.chrom_size} --bin_size {wildcards.bin_size} --output {output.bed}
        """

rule bin_coverage:
    input:
        bg="results/{dataset}/star/coverage/{mapped}/{sample}_{strand}.bg",
        bins="results/{dataset}/star/binned_norm_coverage/bin{bin_size}bp.bed"
    output:
        bg="results/{dataset}/star/binned_norm_coverage/{mapped}/{sample}_{strand}_bin{bin_size}bp.bg"
    resources:
        tmpdir = lambda wildcards: f"tmpdir/bin_coverage_{wildcards.mapped}_{wildcards.sample}_{wildcards.strand}_{wildcards.bin_size}"
    shell:
        """
        if [ ! $(hostname -s) == "haas005" ]; then
            ml gcc/11.3.0 bedtools2
        fi
        
        bedtools map -a {input.bins} -b {input.bg} -c 4 -o sum -null out | awk '$4 != "out"' > {output.bg}
        """

rule binned_norm_coverage_bw:
    input:
        bg="results/{dataset}/star/binned_norm_coverage/{mapped}/{sample}_{strand}_bin{bin_size}bp.bg",
        chrom_size="resources/genome/{spec}/chrom.sizes".format(spec=config['Species'])
    output:
        bw="results/{dataset}/star/binned_norm_coverage/{mapped}/{sample}_{strand}_bin{bin_size}bp.bw"
    params:
        tmp_sorted="results/{dataset}/star/binned_norm_coverage/{mapped}/{sample}_{strand}_bin{bin_size}bp.sorted.tmp"
    shell:
        """
        sort -k1,1 -k2,2n {input.bg} > {params.tmp_sorted}
        bedGraphToBigWig {params.tmp_sorted} {input.chrom_size} {output.bw}
        rm {params.tmp_sorted}
        """

rule bin_expression_table:
    input:
        bgs=expand("results/{{dataset}}/star/binned_norm_coverage/{{mapped}}/{sample}_{{strand}}_bin{{bin_size}}bp.bg",sample=get_SampleNames(config['Name']))
    output:
        table=temp("results/{dataset}/star/binned_norm_coverage/{mapped}/expression_tables/bin_expression_table_{chr}_{strand}_bin{bin_size}bp.csv")
    shell:
        """
        python scripts/make_expression_table_chr_strand.py --chr {wildcards.chr} \
                                                           --bedgraphs {input.bgs} \
                                                           --output {output.table}
        """
rule bin_expression_table_merge_strands:
    input:
        forward_table="results/{dataset}/star/binned_norm_coverage/{mapped}/expression_tables/bin_expression_table_{chr}_forward_bin{bin_size}bp.csv",
        reverse_table="results/{dataset}/star/binned_norm_coverage/{mapped}/expression_tables/bin_expression_table_{chr}_reverse_bin{bin_size}bp.csv"
    output:
        table=temp("results/{dataset}/star/binned_norm_coverage/{mapped}/expression_tables/bin_expression_table_{chr}_bin{bin_size}bp.csv")
    shell:
        """
        python scripts/make_expression_table_chr.py --forward_table {input.forward_table} \
                                                    --reverse_table {input.reverse_table} \
                                                    --output {output.table}
        """
rule bin_expression_table_merge_chr:
    input:
        tables=expand("results/{{dataset}}/star/binned_norm_coverage/{{mapped}}/expression_tables/bin_expression_table_{chr}_bin{{bin_size}}bp.csv",chr=config['Chromosomes'])
    output:
        table="results/{dataset}/star/binned_norm_coverage/{mapped}/expression_tables/bin_expression_table_bin{bin_size}bp.csv"
    threads: 12
    shell:
        """
        python scripts/make_expression_table.py --chr_tables {input.tables} \
                                                --output {output.table}
        """

rule log2_coverage:
    input:
        bg="results/{dataset}/star/binned_norm_coverage/{mapped}/{sample}_{strand}_bin{bin_size}bp.bg"
    output:
        bg="results/{dataset}/star/log_binned_norm_coverage/{mapped}/{sample}_{strand}_bin{bin_size}bp.bg"
    resources:
        tmpdir = lambda wildcards: f"tmpdir/log2_coverage_{wildcards.sample}_{wildcards.strand}_{wildcards.bin_size}.tmp"
    shell:
        """
        awk '{{print $1"\t"$2"\t"$3"\t"log($4+1)/log(2)}}' {input.bg} > {output.bg}
        """

rule log_coverage_bw:
    input:
        bg="results/{dataset}/star/log_binned_norm_coverage/{mapped}/{sample}_{strand}_bin{bin_size}bp.bg",
        chrom_size="resources/genome/{spec}/chrom.sizes".format(spec=config['Species'])
    output:
        bw="results/{dataset}/star/log_binned_norm_coverage/{mapped}/{sample}_{strand}_bin{bin_size}bp.bw"
    params:
        tmp_sorted="results/{dataset}/star/log_binned_norm_coverage/{mapped}/{sample}_{strand}_bin{bin_size}bp.sorted.tmp"
    shell:
        """
        sort -k1,1 -k2,2n {input.bg} > {params.tmp_sorted}
        bedGraphToBigWig {params.tmp_sorted} {input.chrom_size} {output.bw}
        rm {params.tmp_sorted}
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

rule make_tpm_tables:
    input:
        expand('results/{dataset}/kallisto/{sample}/abundance.tsv', sample=get_SampleNames(config['Name']), dataset=config['Name'])
    output:
        mrna='results/{dataset}/kallisto/mrna_tpm_table.tab'
    shell:
        """
        python scripts/make_tpm_table.py --intables {input} --out_mrna {output.mrna}
        """





##----------------------##
## BWA : map to Genome  ##
##----------------------##

#rule bwa_index:
#    input:
#        genome="resources/genome/{spec}/genome.fa"
#    output:
#        bwa_index="resources/genome/{spec}/genome.fa.bwt"
#    shell:
#        """
#        bwa index {input.genome}
#        """
#
#rule bwa_aln_se:
#    input:
#        genome="resources/genome/{spec}/genome.fa".format(spec=config['Species']),
#        genome_index="resources/genome/{spec}/genome.fa.bwt".format(spec=config['Species']),
#        fq="resources/{dataset}/fastq/{sample}.fastq.gz",
#    output:
#        sai="results/{dataset}/bwa/{sample}.sai"
#    threads: 12
#    shell:
#        """
#        bwa aln -t {threads} -f {output.sai} {input.genome} {input.fq}
#        """
#
#rule bwa_samse:
#    input:
#        genome = "resources/genome/{spec}/genome.fa".format(spec=config['Species']),
#        gemone_index = "resources/genome/{spec}/genome.fa.bwt".format(spec=config['Species']),
#        sai = "results/{dataset}/bwa/{sample}.sai",
#        fq = "resources/{dataset}/fastq/{sample}.fastq.gz",
#    output:
#        sam="results/{dataset}/bwa/{sample}.sam"
#    threads: 3
#    shell:
#        """
#        bwa samse -f {output.sam} {input.genome} {input.sai} {input.fq}
#        """
#
#rule mapped_uniq_MAPQ_sorted_se:
#    input:
#        sam="results/{dataset}/bwa/{sample}.sam"
#    output:
#        bam="results/{dataset}/bwa/{sample}.bam"
#    params:
#        MAPQ_min=config['bwa']['MAPQ_min'],
#        mem=config['bwa']['mem']
#    threads: 12
#    shell:
#        """
#        if [ ! $(hostname -s) == "haas005" ]; then
#            ml gcc/11.3.0 samtools
#        fi
#        samtools view -@ {threads} -h -q {params.MAPQ_min} -F 4 {input.sam} |\
#            awk '$1 ~ /^@/ || $2 == "16" || $2 == "0"' |\
#            grep -v 'XA:Z:'| grep -v 'SA:Z:' |\
#            samtools sort -@ {threads} -m {params.mem} > {output.bam}
#        """
#
#rule bwa_bam_index:
#    input:
#        bam="results/{dataset}/bwa/{sample}.bam"
#    output:
#        bam_index="results/{dataset}/bwa/{sample}.bam.bai"
#    shell:
#        """
#        ml gcc/11.3.0 samtools
#        samtools index {input.bam}
#        """
#
#rule bwa_coverage_bedgraph:
#    input:
#        bam="results/{dataset}/bwa/{sample}.bam",
#        bam_index="results/{dataset}/bwa/{sample}.bam.bai"
#    output:
#        bg="results/{dataset}/bwa/coverage/{sample}_{strand}.bedgraph"
#    params:
#        strand = lambda wildcards: '-' if wildcards.strand == "forward" else '+' if wildcards.strand == "reverse" else '',
#        flag_pe = lambda wildcards: '81' if wildcards.strand == "forward" else '97' if wildcards.strand == "reverse" else '',
#        flag_se = lambda wildcards: '16' if wildcards.strand == "forward" else '0' if wildcards.strand == "reverse" else ''
#    shell:
#        """
#        if [ ! $(hostname -s) == "haas005" ]; then
#            ml gcc/11.3.0 samtools bedtools2
#        fi
#        samtools view -h {input.bam} | awk '$1 ~ /^@/ || $2 == "{params.flag_pe}" || $2 == "{params.flag_se}"' | samtools view -b | bedtools genomecov -ibam stdin -bg -strand {params.strand} -5 | grep "^\<chr" > {output.bg}
#        """
#
