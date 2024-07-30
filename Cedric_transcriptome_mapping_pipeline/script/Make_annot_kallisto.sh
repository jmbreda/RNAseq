#!/bin/bash
 # Generate input cDNA fasta for kallisto

GTF=$1
CDNA=$2
DNA=$3
CDNA_OUT=$4

cat $GTF | grep 'transcript_type "protein_coding"' | grep 'transcript_support_level "1"\|transcript_support_level "2"' \
 | awk -v OFS='\t' '$3 == "transcript" {print $1,$4,$5,$12,"1",$7}' > select_transcripts.bed
sed 's/"\|;//g' select_transcripts.bed > select_transcripts_clean.bed

cut -f4 select_transcripts_clean.bed  > transcript_list.tsv 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $CDNA > cdna_linear.tsv
grep -i -E -f transcript_list.tsv  cdna_linear.tsv > cdna_linear_sub.tsv
tr "\t" "\n" <  cdna_linear_sub.tsv >  cdna_sub.fa

bedtools getfasta -fi $DNA -bed select_transcripts_clean.bed -fo premRNA.fa -name   

sed 's/>/>pre_/g' premRNA.fa  > premRNA_modif.fa 
cat premRNA_modif.fa  cdna_sub.fa > $CDNA_OUT


 
