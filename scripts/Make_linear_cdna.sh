#!/bin/bash
 # Generate input cDNA fasta for kallisto

CDNA=$1
CDNA_OUT=$2


# get cDNA fasta with one transcript per line
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "$CDNA" > "$CDNA_OUT"
