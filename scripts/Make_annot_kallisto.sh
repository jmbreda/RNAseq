#!/bin/bash
 # Generate input cDNA fasta for kallisto

GTF=$1
CDNA=$2
DNA=$3
CDNA_OUT=$4

# get only protein coding transcripts with support level 1 and 2
grep 'transcript_type "protein_coding"' "${GTF}" \
    | grep 'transcript_support_level "1"\|transcript_support_level "2"' \
    | awk -v OFS='\t' '$3 == "transcript" {print $1,$4,$5,$12,"1",$7}' \
    | sed 's/"\|;//g' > "${GTF}".select_transcripts.bed

# get cDNA fasta with one transcript per line
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "$CDNA" > "$CDNA".linear.fa

# get cDNA fasta for selected transcripts
grep -i -E -f <( cut -f4 "${GTF}".select_transcripts.bed ) "$CDNA".linear.fa \
    | tr "\t" "\n" > "$CDNA".linear.selected.fa

# get pre-mRNA fasta for selected transcripts and add prefix to fasta headers
bedtools getfasta -fi "$DNA" -bed "${GTF}".select_transcripts.bed -name | sed 's/>/>pre_/g' > "$DNA".premRNA.fa

# Add pre-mRNA and cDNA fasta together to get input fasta for kallisto
cat "$DNA".premRNA.fa "$CDNA".linear.selected.fa > "$CDNA_OUT"

# clean up
#rm "${GTF}".select_transcripts.bed
#rm "${CDNA}".linear.fa
#rm "${CDNA}".linear.selected.fa
#rm "${DNA}".premRNA.fa
 
