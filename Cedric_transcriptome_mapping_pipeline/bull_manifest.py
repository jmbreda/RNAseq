#!/usr/bin/env python

import json
import pandas as pd
import os 

data = json.load(open('/home/cgobet/270123_IE_GTEX/script/file-manifest.json','r'))
df = pd.read_csv('/home/cgobet/270123_IE_GTEX/snakemake/Adipose_Visceral.tsv',sep="\t")
fn = "/home/cgobet/270123_IE_GTEX/snakemake/manifest/Adipose_Visceral.json"
def create_manifest_ID(data, df):
    my_json = []

    for v in df['ID']:
        for k in data:
            if((v in k['file_name']) and not k['file_name'].endswith("bai") and 'H3K27ac' not in k['file_name']):
                my_json.append(k)

    json.dump(my_json,open(fn,'w'), indent=4)


dd= create_manifest_ID(data, df)
