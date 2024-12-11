import argparse
import pandas as pd

def parse_argument():
    parser = argparse.ArgumentParser(description='Get SRR, SampleName table from SRA metadata and GSMID_SampleName table.')
    parser.add_argument('--SraRunTable'
        ,required=True
        ,type=str
        ,help="SraRunTable.txt")
    parser.add_argument('--GSMID_SampleName'
        ,required=True
        ,type=str
        ,help="GSMID_SampleName.txt")
    parser.add_argument('--outfile'
        ,required=True
        ,type=str
        ,help="Output table file")
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_argument()

    # load tables
    SraRunTable = pd.read_csv(args.SraRunTable, sep=',')
    GSMID_SampleName = pd.read_csv(args.GSMID_SampleName, sep=' ',header=None,names=['GSMID','SampleName'])

    # keep only RNA-Seq data
    if 'Assay Type' in SraRunTable.columns:
        SraRunTable = SraRunTable.loc[SraRunTable['Assay Type'] == 'RNA-Seq']

    # rename GEO_Accession (exp) to GSMID
    if 'GEO_Accession (exp)' in SraRunTable.columns:
        SraRunTable = SraRunTable.rename(columns={'GEO_Accession (exp)':'GSMID'})
    elif 'Library Name' in SraRunTable.columns:
        SraRunTable = SraRunTable.rename(columns={'Library Name':'GSMID'})

    # merge tables
    SraRunTable = pd.merge(SraRunTable, GSMID_SampleName, left_on='GSMID', right_on='GSMID')
    SampleName_SRR = SraRunTable.loc[:,['Run','GSMID','SampleName']].groupby('SampleName')['Run'].apply(lambda x: ','.join(x)).reset_index()
    
    # sort by Run
    SampleName_SRR = SampleName_SRR.sort_values(by='Run')
    
    # save
    SampleName_SRR.to_csv(args.outfile, sep='\t', index=False)