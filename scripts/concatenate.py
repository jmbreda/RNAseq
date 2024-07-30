import argparse
import pandas as pd
import os

def parse_argument():
    parser = argparse.ArgumentParser(description='Get SRR, SampleName table from SRA metadata and GSMID_SampleName table.')
    parser.add_argument('--SampleName'
        ,required=True
        ,type=str
        ,help="Sample name")
    parser.add_argument('--srr_list'
        ,required=True
        ,type=str
        ,help="comma separated SRR list in one string")
    parser.add_argument('--infolder'
        ,required=True
        ,type=str
        ,help="input folder")
    parser.add_argument('--fq1'
        ,required=True
        ,type=str
        ,help="concatenated fastq file R1")
    parser.add_argument('--fq2'
        ,required=False
        ,type=str
        ,help="concatenated fastq file R2 (leave empty if single end)")
    
    return parser.parse_args()

if __name__ == '__main__':
    
    args = parse_argument()

    srr_list = args.srr_list.split(',')

    # concatenate fastq files
    if not os.path.exists(args.fq1):
        bashCommand = "cat " + ' '.join([args.infolder + '/' + srr + '_1.fastq' for srr in srr_list]) + " | gzip > " + args.fq1
        output = os.system(bashCommand)
        print(output)

    # concatenate fastq files for R2 if exists
    if args.fq2:
        if not os.path.exists(args.fq2):
            bashCommand = "cat " + ' '.join([args.infolder + '/' + srr + '_2.fastq' for srr in srr_list]) + " | gzip > " + args.fq2
            output = os.system(bashCommand)
            print(output)
