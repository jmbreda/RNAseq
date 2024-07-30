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

    # merge all tables
    premrna = pd.DataFrame()
    mrna = pd.DataFrame()
    for table in args.intables:
        print(table)
        tmp = pd.read_csv(table, sep='\t')
        N = tmp.shape[0]
        premrna = pd.concat([premrna, tmp.loc[:N//2-1,'tpm']], axis=1)
        mrna = pd.concat([mrna, tmp.loc[N//2:,'tpm']], axis=1)

    cols = [f.split('/')[-2] for f in args.intables]
    premrna.columns = cols
    mrna.columns = cols
    premrna.index = tmp.loc[:N//2-1,'target_id']
    mrna.index = tmp.loc[N//2:,'target_id']

    premrna.to_csv(args.out_premrna, sep='\t')
    mrna.to_csv(args.out_mrna, sep='\t')
            