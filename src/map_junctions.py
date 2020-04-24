#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" splice_junctions.py

DESCRIPTION 

This file can also be imported as a module and contains the following functions:
    * 

TO DO:  
    *
"""

from __future__ import absolute_import, division, print_function

import argparse, glob, os 
import pandas as pd
from . import star_sj

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
    # parser.add_argument('-f', '--file', type=str, help='Custom or not GENCODE gff file.')                        
    args = parser.parse_args()

    # concatenating several SJ.out after being annotated and parsed in same pandas DataFrame
    globdir = f'/media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g{args.gencode}/ER*.1/SJ.out.tab'
    df_sj = star_sj.concat_samples(globdir)

    # Getting max values per position and per tissue group sample
    df_sj['tissue_group'] = df_sj['tissue'].str.split('_').str[0]
    df_sj_max_position = df_sj.sort_values(by=['start', 'end', 'unique_reads'], ascending=[True, True, False]).drop_duplicates(subset=['start', 'end'], keep='first').reset_index(drop=True)
    df_sj_max_position.to_csv('../data/processed/sj_maxp.emtab2836.tsv.gz', index=None, sep='\t', compression='gzip')
    df_sj_max_tissue = df_sj.sort_values(by=['start', 'end', 'unique_reads'], ascending=[True, True, False]).drop_duplicates(subset=['start', 'end', 'tissue_group'], keep='first').reset_index(drop=True)
    df_sj_max_tissue.to_csv('../data/processed/sj_maxt.emtab2836.tsv.gz', index=None, sep='\t', compression='gzip')

    # Loading previously annotated introns
    df_introns_annotated = pd.read_csv(f'../data/processed/gencode.v{args.gencode}.annotation.complete.tsv.gz', sep='\t', compression='gzip')
    df_introns_annotated = df_introns_annotated[df_introns_annotated['type']=='intron'].reset_index(drop=True)

    # Mapping splice junctions read to introns positions
    df_qsj = pd.merge(df_introns_annotated, df_sj_max_position, how='left', on=['seq_id', 'start', 'end'])
    df_qsj.loc[df_qsj['unique_reads'].isnull(), 'unique_reads'] = 0
    df_qsj.loc[df_qsj['unique_reads'].isnull(), 'tissue'] = '-'
    df_qsj = df_qsj.sort_values(by=['seq_id', 'gene_id', 'transcript_id', 'start', 'end'], ascending=[True, True, True, True, True]).reset_index(drop=True).fillna('-')

    # Calculating means and score
    df_qsj = df_qsj[df_qsj['transcript_type'].str.contains('protein_coding|nonsense_mediated_decay')].reset_index(drop=True)
    df_qsj['gene_mean'] = df_qsj.groupby('gene_id')['unique_reads'].transform(lambda x:x.mean())
    df_qsj['gene_mean_cds'] = df_qsj.loc[df_qsj['cds_coverage'] == 'full'].groupby('gene_id')['unique_reads'].transform(lambda x:x.mean())

    df_qsj['RNA2sj'] = df_qsj['unique_reads']/df_qsj['gene_mean']
    df_qsj['norm_RNA2sj'] = df_qsj.groupby(['gene_id'])['RNA2sj'].transform(lambda x: (x-x.min()) / (x.max()-x.min()))

    df_qsj['RNA2sj_cds'] = df_qsj['unique_reads']/df_qsj['gene_mean_cds']
    df_qsj['norm_RNA2sj_cds'] = df_qsj.groupby(['gene_id'])['RNA2sj_cds'].transform(lambda x: (x-x.min()) / (x.max()-x.min()))
    df_qsj.to_csv(f'../data/processed/sj_maxp.emtab2836.v{args.gencode}.tsv.gz', sep='\t', index=None, compression='gzip')

    # Calculating mins and exporting qsplice
    df_qsj.loc[df_qsj['cds_coverage'] == 'none', 'unique_reads'] = df_qsj['unique_reads'].max()
    df_qsj = df_qsj.loc[df_qsj.groupby('transcript_id')['unique_reads'].idxmin()].sort_values(by=['seq_id', 'gene_id', 'transcript_id', 'start', 'end'], ascending=[True, True, True, True, True]).reset_index(drop=True).fillna('-')
    df_qsj = df_qsj[['seq_id','transcript_id','transcript_type','strand','gene_id','gene_name','gene_type','intron_number','start','end','nexons','ncds','unique_reads','tissue','tissue_group','gene_mean','gene_mean_cds','RNA2sj','norm_RNA2sj','RNA2sj_cds','norm_RNA2sj_cds']]
    df_qsj.to_csv(f'../data/processed/qsplice.emtab2836.v{args.gencode}.tsv.gz', sep='\t', index=None, compression='gzip')     

if __name__ == "__main__":
    main()
