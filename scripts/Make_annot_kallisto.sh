#!/bin/bash
 # Generate input cDNA fasta for kallisto

GTF=$1
CDNA=$2
DNA=$3
CDNA_OUT=$4
TRANSCRIPTOME=$5

if [ "$TRANSCRIPTOME" = "selected_transcripts" ]; then
    echo "Selected transcripts"

    # get only protein coding transcripts with support level 1 and 2
    grep 'transcript_type "protein_coding"' "${GTF}" \
        | grep 'transcript_support_level "1"\|transcript_support_level "2"' \
        | awk -v OFS='\t' '$3 == "transcript" {print $1,$4,$5,$12,"1",$7}' \
        | sed 's/"\|;//g' > "${GTF}"."${TRANSCRIPTOME}".bed

elif [ "$TRANSCRIPTOME" = "all_transcripts" ]; then
    echo "All transcripts"

    # get only transcripts bed file
    awk -v OFS='\t' '$3 == "transcript" {print $1,$4,$5,$12,"1",$7}' "${GTF}" \
        | sed 's/"\|;//g' > "${GTF}"."${TRANSCRIPTOME}".bed
else
    echo "Invalid option for TRANSCRIPTOME"
    exit 1
fi

# get cDNA fasta for selected transcripts
grep -i -E -f <( cut -f4 "${GTF}"."${TRANSCRIPTOME}".bed ) "$CDNA" \
    | tr "\t" "\n" > "${CDNA}"."${TRANSCRIPTOME}".fa

# get pre-mRNA fasta for selected transcripts and add prefix to fasta headers
bedtools getfasta -fi "$DNA" -bed  "${GTF}"."${TRANSCRIPTOME}".bed -name | sed 's/>/>pre_/g' > "$DNA"."${TRANSCRIPTOME}".premRNA.fa

# Add pre-mRNA and cDNA fasta together to get input fasta for kallisto
cat "$DNA"."${TRANSCRIPTOME}".premRNA.fa "${CDNA}"."${TRANSCRIPTOME}".fa > "$CDNA_OUT"

# clean up
#rm "${GTF}".select_transcripts.bed
#rm "${CDNA}".linear.fa
#rm "${CDNA}".linear.selected.fa
#rm "${DNA}".premRNA.fa