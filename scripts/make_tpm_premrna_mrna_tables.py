import numpy as np
import pandas as pd
import argparse

def parse_argument():
    parser = argparse.ArgumentParser(description='Get transcript to promoter file')
    parser.add_argument('--intables', type=str, nargs='+', help='Input tables')
    parser.add_argument('--out_premrna', type=str, help='Output pre-mRNA table',default='tmp_premrna.tsv')
    parser.add_argument('--out_mrna', type=str, help='Output mRNA table',default='tmp_mrna.tsv')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_argument()

    # concatenate all tables
    abundance = pd.DataFrame()
    for table in args.intables:
        tmp = pd.read_csv(table, sep='\t', index_col=0, usecols=['target_id', 'tpm'])
        tmp.columns = [table.split('/')[-2]]
        abundance = pd.concat([abundance, tmp], axis=1)
    del tmp

    # split between mRNA and pre-mRNA
    premrna_idx = [id[:11]=="pre_ENSMUST" for id in abundance.index]
    mrna_idx = [id[:7]=="ENSMUST" for id in abundance.index]

    # check that all transcripts are either pre-mRNA or mRNA
    assert sum(premrna_idx) + sum(mrna_idx) == len(abundance.index)
    premrna = abundance.loc[premrna_idx,:]
    mrna = abundance.loc[mrna_idx,:]

    # save tables
    premrna.to_csv(args.out_premrna, sep='\t')
    mrna.to_csv(args.out_mrna, sep='\t')
    