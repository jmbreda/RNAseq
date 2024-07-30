#!/bin/bash
# Generate intron-exon annotations for unstranded paired-end RNA-seq

GTF=$1
OUT_in=$2
OUT_ex=$3
OUT_ex_gene=$4


cat ${GTF} | grep 'transcript_type "protein_coding"' | grep 'transcript_support_level "1"\|transcript_support_level "2"' | awk -v OFS='\t' '$3 == "transcript" {print $1,$4,$5,$10,"1",$7}' > ${GTF}.gene.bed
cat ${GTF} | grep 'transcript_type "protein_coding"' | grep 'transcript_support_level "1"\|transcript_support_level "2"' | awk -v OFS='\t' '$3 == "exon" {print $1,$4,$5,$10,"1",$7}' > ${GTF}.exon.bed

sed 's/"\|;//g' ${GTF}.gene.bed > ${GTF}.gene.2.bed
sed 's/"\|;//g' ${GTF}.exon.bed > ${GTF}.exon.2.bed

sort -k1,1 -k4,4 -k2,2n ${GTF}.gene.2.bed   > ${GTF}.sorted.gene.2.bed

groupBy -g 1,4,5,6 -c 2,3 -o min,max -i ${GTF}.sorted.gene.2.bed |  awk -v OFS='\t' '{print $1, $5, $6, $2, $3, $4}' > ${GTF}.gene.2.merge.bed 

 #check overlapping region between genes coordinate and itself, remove the "self" intersection
bedtools intersect -wa -wb -a ${GTF}.gene.2.merge.bed -b ${GTF}.gene.2.merge.bed > ${GTF}_inter_wa_wb.bed
bedtools intersect  -a ${GTF}.gene.2.merge.bed -b ${GTF}.gene.2.merge.bed > ${GTF}_inter.bed
paste ${GTF}_inter_wa_wb.bed ${GTF}_inter.bed > ${GTF}_inter_all.bed
cat ${GTF}_inter_all.bed |  awk '($4!=$10)' | cut -f13,14,15,16,17,18 > ${GTF}.gene.2.bed_inter

# remove the overlapping regions to the default gene bed file
subtractBed -a ${GTF}.gene.2.merge.bed -b ${GTF}.gene.2.bed_inter -s > ${GTF}.gene.3.bed 

# generate exon and intron coordinates
subtractBed -a ${GTF}.exon.2.bed -b ${GTF}.gene.2.bed_inter -s > ${GTF}.exon.3.bed
subtractBed -a ${GTF}.gene.3.bed -b ${GTF}.exon.3.bed -s > ${GTF}.introns.bed
subtractBed -a ${GTF}.gene.3.bed -b ${GTF}.introns.bed -s > ${GTF}.exons_merge.bed

cp ${GTF}.introns.bed  $OUT_in
cp  ${GTF}.exons_merge.bed  $OUT_ex
cp ${GTF}.gene.3.bed $OUT_ex_gene


#GTF=ensembl.gtf
#bash /home/cgobet/270123_IE_GTEX/snakemake/script/Make_annot.sh /scratch/cgobet/270123_IE_GTEX_2/REF_2/Human/${GTF} \ /scratch/cgobet/270123_IE_GTEX_2/REF_2/Human/${GTF}.introns.final.bed /scratch/cgobet/270123_IE_GTEX_2/REF_2/Human/${GTF}.exons.final.bed /scratch/cgobet/270123_IE_GTEX_2/REF_2/Human/${GTF}.genes.final.bed /scratch/cgobet/270123_IE_GTEX_2/REF_2/Human/gencode.v43.pc_transcripts.fa /scratch/cgobet/270123_IE_GTEX_2/REF_2/Human/gencode.v43.pc_transcripts.fa.sub 