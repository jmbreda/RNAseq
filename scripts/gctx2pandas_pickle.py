import argparse
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse


def parse_argument():
    parser = argparse.ArgumentParser(description='Save ChIP signals for tf in tensor (N_promoters x N_positions X N_experiments).')
    parser.add_argument('-i','--infile'
        ,required=True
        ,type=str
        ,help="input gctx file")
    parser.add_argument('-o','--outfile'
        ,required=True
        ,type=str
        ,help="Output pandas dataframe pickle file")
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_argument()
    df = parse(args.infile).data_df.to_pickle(args.outfile)


