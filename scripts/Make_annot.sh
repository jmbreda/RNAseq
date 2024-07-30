#!/bin/bash
# Generate intron-exon annotations for unstranded paired-end RNA-seq

GTF=$1
OUT_intron=$2
OUT_exon=$3
OUT_ex_gene=$4

# get only protein coding transcripts and exons from gtf
grep 'transcript_type "protein_coding"' "${GTF}" \
    | grep 'transcript_support_level "1"\|transcript_support_level "2"' \
    | awk -v OFS='\t' '$3 == "transcript" {print $1,$4,$5,$10,"1",$7}' \
    | sed 's/"\|;//g' > "${GTF}".gene.bed
grep 'transcript_type "protein_coding"' "${GTF}" \
    | grep 'transcript_support_level "1"\|transcript_support_level "2"' \
    | awk -v OFS='\t' '$3 == "exon" {print $1,$4,$5,$10,"1",$7}' \
    | sed 's/"\|;//g' > "${GTF}".exon.bed

# sort by chr, transcript name and coordinate, then merge same transcript and take min-max coordinate
sort -k1,1 -k4,4 -k2,2n "${GTF}".gene.bed \
    | bedtools groupby -g 1,4,5,6 -c 2,3 -o min,max \
    | awk -v OFS='\t' '{print $1, $5, $6, $2, $3, $4}' > "${GTF}".gene.merge.bed

# check overlapping region between genes coordinate and itself, remove the "self" intersection
bedtools intersect -wa -wb -a "${GTF}".gene.merge.bed -b "${GTF}".gene.merge.bed > "${GTF}".inter_wa_wb.bed
bedtools intersect  -a "${GTF}".gene.merge.bed -b "${GTF}".gene.merge.bed > "${GTF}".inter_ab.bed
paste "${GTF}"_inter.wa_wb.bed "${GTF}".inter_ab.bed \
    | awk '($4!=$10)' \
    | cut -f13,14,15,16,17,18 > "${GTF}".gene.self_intersection.bed

# remove the self overlapping regions to the default gene bed file
bedtools subtract -a "${GTF}".gene.merge.bed -b "${GTF}".gene.self_intersection.bed -s > "$OUT_ex_gene"

# generate exon and intron coordinates
bedtools subtract -a "${GTF}".exon.bed -b "${GTF}".gene.self_intersection.bed -s > "${GTF}".exon.clean.bed
bedtools subtract -a "$OUT_ex_gene"    -b "${GTF}".exon.clean.bed -s             > "$OUT_intron"
bedtools subtract -a "$OUT_ex_gene"    -b "$OUT_intron" -s                       > "$OUT_exon"

# clean up
rm "${GTF}".gene.bed
rm "${GTF}".exon.bed
rm "${GTF}".gene.merge.bed
rm "${GTF}".inter*.bed
rm "${GTF}".gene.self_intersection.bed
rm "${GTF}".exon.clean.bed

