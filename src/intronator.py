#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" generate_introns.py

Usage: python intronator.py --gencode 27
___
--help | -h Display documentation

[description] 

This file can also be imported as a module and contains the following functions:
    * 

TO DO:  
    *
"""

from __future__ import absolute_import, division, print_function

import argparse, os, re, subprocess
import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd
from . import source

__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Production"


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-g', '--gencode', type=int, help='GENCODE version selected.')
    parser.add_argument('-f', '--file', type=str, help='Custom or not GENCODE gff file.')                        
    args = parser.parse_args()

    RAWDIR = 'data/raw'
    INTRONSDIR = 'data/interim'
    PROCESSEDDIR = 'data/processed'

    if args.gencode:
        path = source.Gencode(version=args.gencode, specie='human').download(outdir=RAWDIR, type='gff3')
    else:
        path = args.file

    gff_path = generate_introns(path, intronsdir=INTRONSDIR) # generating introns from gff

    df_introns = gff2pandas(gff_path) # converting in pandas data structure 
    df_introns.to_csv(os.path.join(INTRONSDIR, gff_path), index=None, sep='\t') 

    df_introns_annotated = annotate_introns(df_introns)
    df_introns_annotated.to_csv(os.path.join(PROCESSEDDIR, gff_path), index=None, sep='\t')


def generate_introns(inpath: str, intronsdir: str = '.') -> str:
    """This function generates introns file from gff*.gz file path.
    
    Arguments:
        inpath {str} -- The name of the annotation file to process.
    
    Keyword Arguments:
        intronsdir {str} -- Introns directory (default: {'.'})
    
    Returns:
        str -- Path where file has been stored.
    """    
    subprocess.call('which gt || apt install gt', shell=True)
    subprocess.call('which rg || apt install rg', shell=True)

    source.create_dir(intronsdir)
    intronsname = os.path.basename(re.sub(r'gff.*', 'introns.tsv', inpath))
    outpath = os.path.abspath(os.path.join(intronsdir, intronsname))

    # cmd_awk = r"awk -F'[\t=]' '{print $1,$4,$5,$7,$10}' OFS='\t'"
    # cmd_gt = f"zcat {inpath} | rg '(#|gene_type=protein_coding)' | gt gff3 -retainids -addintrons | rg 'intron' | {cmd_awk} > {outpath}"
    cmd_gt = f"zcat {inpath} | gt gff3 -tidy -retainids -addintrons > {outpath}" # tidy option described here http://genometools.org/pipermail/gt-users/2015-August/000794.html
    subprocess.call(cmd_gt, shell=True)
    print(f'Introns has been generated in: {os.path.dirname(outpath)}')
    return outpath


def gff2pandas(inpath: str) -> str:
    """Convert from gff* to pandas.
    
    Arguments:
        inpath {str} -- gff* path.
    
    Returns:
        df {list} -- gff* pandas DataFrame.
    """    
    df = gffpd.read_gff3(inpath).attributes_to_columns()
    df = df[['seq_id', 'source', 'type', 'start', 'end', 'strand', 'exon_id', 'exon_number', 'gene_name', 'gene_id', 'gene_type', 'transcript_id', 'transcript_type']]
    print(f'Introns has been processed and stored here: {os.path.abspath(inpath)}\nPlease, wait while our script is saving the data.')
    return df


# database processing
def add_introns_info(df: list) -> list:
    """Adding intron source, number and correct positions. (corrected by genome strand)
    
    Arguments:
        df {list} -- pandas DataFrame to input from initial introns annotation.
    
    Returns:
        list -- pandas DataFrame to output with introns info added.
    """    
    df.loc[df['type'] == 'intron', 'source'] = 'gt'
    df['intron_number'] = df['exon_number']
    intronf_list = ['gene_name','gene_id','gene_type','transcript_id','transcript_type']
    for feature in intronf_list:
        df.loc[df['type'].str.contains('exon|intron'), feature] = df.fillna(method='ffill')
    df.loc[(df['strand']=='+') & (df['type'].str.contains('exon|intron')), 'intron_number'] = df.fillna(method='ffill')
    df.loc[(df['strand']=='-') & (df['type'].str.contains('exon|intron')), 'intron_number'] = df.fillna(method='bfill')
    df.loc[~df['type'].str.contains('intron'), 'intron_number'] = '-'
    return df


def get_exons_cds(df: list) -> list:
    """Getting exons with cds annotations.
    
    Arguments:
        df {list} -- pandas DataFrame to input from initial introns annotation.
    
    Returns:
        list -- pandas DataFrame with exons annotations
    """    
    exon = df[df['type'] == 'exon'].drop(['intron_number', 'type'], axis=1).reset_index(drop=True)
    cds = df[df['type'] == 'CDS'].drop(['intron_number', 'type'], axis=1).reset_index(drop=True)

    df = pd.merge(exon, cds, how='left', 
                     on=list(df.drop(['start', 'end', 'type', 'intron_number'], axis=1).columns), 
                     indicator='cds_coverage')
    df['cds_coverage'] = df['cds_coverage'].replace('left_only', 'none')
    df['cds_coverage'] = df['cds_coverage'].replace('both', 'full')
    df.insert(2, 'type', 'exon')
    
    df = df.rename(columns={'start_x':'start', 
                            'end_x': 'end',
                            'start_y': 'start_cds',
                            'end_y': 'end_cds'})
    
    df.loc[df['start_cds'].isnull(), 'start_cds'] = df['start']
    df.loc[df['end_cds'].isnull(), 'end_cds'] = df['end']
    
    df['start_cds'] = df['start_cds'].astype(int) 
    df['end_cds'] = df['end_cds'].astype(int)
    
    df['cds_coverage'] = df['cds_coverage'].astype(str)
    df.loc[(df['start_cds'] != np.nan) & ((df['start_cds'] != df['start']) | (df['end_cds'] != df['end'])), 'cds_coverage'] = 'partial'
    return df


def annotate_introns(df: list) -> list:
    """Annotating introns cds positions merging exons and cds info.
    
    Arguments:
        df {list} -- pandas DataFrame to input from whole introns annotation.
    
    Returns:
        list -- pandas DataFrame to finally output.
    """    
    df_whole = add_introns_info(df)
    
    df_exons = get_exons_cds(df_whole)
    df_introns = df_whole[df_whole['type'] == 'intron'].reset_index(drop=True)

    df = pd.concat([df_exons, df_introns], axis=0).reset_index(drop=True).sort_values(by=['gene_name', 'transcript_id', 'start', 'end']).reset_index(drop=True)
    df['nexons'] = df.groupby('transcript_id')['type'].transform(lambda x: x.str.contains('exon').sum())
    df['ncds'] = df.groupby('transcript_id')['cds_coverage'].transform(lambda x: x.str.contains('full|partial').sum())
    df.loc[df['cds_coverage'].isnull(), 'l_intron_coverage'] = df.groupby('transcript_id')['cds_coverage'].shift(1)
    df.loc[df['cds_coverage'].isnull(), 'r_intron_coverage'] = df.groupby('transcript_id')['cds_coverage'].shift(-1)

    df.loc[(df['cds_coverage'].isnull()) & (df['l_intron_coverage'].str.contains('none') | df['r_intron_coverage'].str.contains('none')), 'cds_coverage'] = 'none'
    df.loc[(df['cds_coverage'].isnull()) & (df['l_intron_coverage'].str.contains('full|partial') | df['r_intron_coverage'].str.contains('full|partial')), 'cds_coverage'] = 'full' 
    df = df.drop(['l_intron_coverage', 'r_intron_coverage'], axis=1).fillna('-')
    return df


if __name__ == "__main__":
    main()
