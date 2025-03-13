#!/usr/bin/env python

import xmltodict
import pandas as pd
import os
import sys

meme_fil = sys.argv[1]

program_lower = os.path.basename(meme_fil).replace('.xml', '')
program_upper = program_lower.upper()
outdir = os.path.dirname(meme_fil)

with open(meme_fil) as f:
    xml_content = f.read()
    xml_meme = xmltodict.parse(xml_content)

# parse xml files, different for meme vs streme
def meme_parse(xml_motif, rank, outdir):
    alph_arr = xml_motif['probabilities']["alphabet_matrix"]["alphabet_array"]
    # loop through each entry in alph_arr[i], access entry 'value', 
    # turn into a dataframe, transpose, and join with previous dataframes
    for j in range(len(alph_arr)):
        meme_df = pd.DataFrame(alph_arr[j]['value']).transpose()
        # use first row @letter_id as column names, remove first row and @letter_id column
        meme_df.columns = meme_df.iloc[0]
        meme_df = meme_df[1:]
        if j == 0:
            meme_df_all = meme_df   
        else:
            meme_df_all = pd.concat([meme_df_all, meme_df], axis=0)
    meme_df_all.to_csv(f"{outdir}/meme_pfm_{rank}.pfm", index=False,
    header=False, sep='\t')

def streme_parse(xml_motif, rank, outdir):
    df = pd.DataFrame(xml_motif['pos'])
    df.to_csv(f"{outdir}/streme_pfm_{rank}.pfm", index=False, 
    header=False, sep='\t')

n_motifs = len(xml_meme[program_upper]['motifs']['motif'])
if n_motifs > 1:
    # parse and write first three motifs
    for i in range(3):
        rank = i + 1
        xml_motif = xml_meme[program_upper]['motifs']['motif'][i]
        if program_lower == 'meme':
            meme_parse(xml_motif, rank, outdir)
        elif program_lower == 'streme':
            streme_parse(xml_motif, rank, outdir)
else:
    xml_motif = xml_meme[program_upper]['motifs']['motif']
    if program_lower == 'meme':
        meme_parse(xml_motif, 1, outdir)
    elif program_lower == 'streme':
        streme_parse(xml_motif, 1, outdir)