import pandas as pd
import numpy as np
import zipfile
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Get transcript to promoter file')
    parser.add_argument('--promid_transcripts', type=str, help='Promoter gtf file')
    parser.add_argument('--premrna', type=str, help='pre-mRNA tmp trable')
    parser.add_argument('--mrna', type=str, help='mRNA tmp table')
    parser.add_argument('--output_premrna', type=str, help='Output prom-premrna table',default='tmp_prom_premrna.tsv')
    parser.add_argument('--output_mrna', type=str, help='Output prom-mrna table',default='tmp_prom_mrna.tsv')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    # Read in promoterome
    promid_transcripts = pd.read_csv(args.promid_transcripts, sep='\t', header=None)
    promid_transcripts.columns = ['prom','transcript']

    def remove_version(x):
        return x.split('.')[0]
    promid_transcripts.loc[:,'transcript'] = promid_transcripts.loc[:,['transcript']].map(remove_version)

    # fill in nan promoters
    for i in promid_transcripts.index:
        if pd.isnull(promid_transcripts.iloc[i, 0]):
            promid_transcripts.iloc[i, 0] = promid_transcripts.iloc[i-1, 0]

    # verify that all transcripts are unique
    assert len(promid_transcripts.transcript.unique()) == len(promid_transcripts)
    
    # make dictionary
    #my_dict = promid_transcripts.set_index('transcript').to_dict()
    my_dict = dict(zip(promid_transcripts.transcript.values, promid_transcripts.prom.values))

    # Read in kallisto tpm tables
    premrna = pd.read_csv(args.premrna, sep='\t', index_col=0)
    mrna = pd.read_csv(args.mrna, sep='\t', index_col=0)
    
    # get transcript ids as index
    premrna.index = [id.split('::')[0].split('_')[1] for id in premrna.index]
    mrna.index = [id.split('|')[0] for id in mrna.index]
    premrna.index = premrna.index.map(remove_version)
    mrna.index = mrna.index.map(remove_version)
    assert np.all(premrna.index == mrna.index)

    # check if how many transcripts are in the promoterome
    idx_in = np.array([id in my_dict for id in premrna.index])
    print('Fraction of transcripts in promoterome: ', np.sum(idx_in)/premrna.shape[0])
    
    # map to promoter and remove transcripts not in promoterome
    premrna.loc[:,'prom'] = premrna.index.map(my_dict)
    mrna.loc[:,'prom'] = mrna.index.map(my_dict)
    premrna = premrna.loc[idx_in,:]
    mrna = mrna.loc[idx_in,:]

    # rename promoter (remove species and version number)
    def rename_prom(p):
        return '_'.join(p.split('_')[2:])
    premrna.loc[:,'prom'] = premrna.loc[:,'prom'].map(rename_prom)
    mrna.loc[:,'prom'] = mrna.loc[:,'prom'].map(rename_prom)

    # sum same promoter
    premrna = premrna.groupby('prom').sum()
    mrna = mrna.groupby('prom').sum()

    # Count fraction of count mapped to promoter
    print('Fraction of count mapped to promoter: ', ( premrna.values.sum() + mrna.values.sum() )/( premrna.shape[1]*1e6 ) )

    # save to file
    premrna.to_csv(args.output_premrna, sep='\t')
    mrna.to_csv(args.output_mrna, sep='\t')


