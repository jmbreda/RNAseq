#!/usr/bin/env python

import json
import pandas as pd
import os 

data = json.load(open('/home/cgobet/270123_IE_GTEX/script/file-manifest.json','r'))
df = pd.read_csv('/home/cgobet/270123_IE_GTEX/snakemake/ID.tsv',sep="\t")
homedir = "/home/cgobet/270123_IE_GTEX/snakemake/"
def create_manifest_ID(data, df):

    my_bam = []
    for v in df['ID']:
        my_json = []
        for k in data:
            if((v in k['file_name']) and not k['file_name'].endswith("bai") and 'H3K27ac' not in k['file_name']):
                mb = k['file_name']
                size = len(mb)
                mod_mb = mb[:size - 4]
                my_bam.append(mod_mb)
                fn = os.path.join(homedir + "manifest_2/file-manifest_") + k['file_name'].split('.')[0] + ".json"
                json.dump(k,open(fn,'w'), indent=4)
                with open(fn, "r") as file_object:
                    re=file_object.read()
                with open(fn, "w") as file_object:
                    file_object.write("[" + re + "]")


    return my_bam

dd= create_manifest_ID(data, df)