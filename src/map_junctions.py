#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" map_junctions.py

Usage: 
python -m src.map_junctions --version g27 --custom --file /media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g27/SJ.out.tab.concat.gz
python -m src.map_junctions --version g33 --glob /media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g29/ER*.1/SJ.out.tab

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
    parser.add_argument('-v', '--version', type=str, help='Genome version selected. (g=GENCODE)')
    parser.add_argument('-g', '--globdir', type=str, help='Directory which contains the files to be globbed and concatenated', 
                        default='/media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g29')
    parser.add_argument('-c', '--custom', help='If you want to customize your file.', action='store_true', default=False)
    parser.add_argument('-f', '--file', type=str, help='Custom splice junctions file (in gzip)',
                        default='/media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g27/SJ.out.tab.concat.gz')
    args = parser.parse_args()

    if args.custom:
        # cat /media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g27/ERR*/SJ.out.tab > ./SJ.out.tab.concat && gzip SJ.out.concat 
        df_sj = pd.read_csv(args.file, compression='gzip', sep='\t', 
                            names=['seq_id', 'start', 'end', 'nstrand', 'unique_reads', 'tissue'])
    else:
        # concatenating several SJ.out after being annotated and parsed in same pandas DataFrame
        df_sj = star_sj.concat_samples(f'{args.globdir}/*/SJ.out.tab')
        # Getting max values per position and per tissue group sample
        df_sj['tissue'] = df_sj['tissue'].str.split('_').str[0]
        
    df_sj_max_position = df_sj.sort_values(by=['start', 'end', 'unique_reads'], ascending=[True, True, False]).drop_duplicates(subset=['start', 'end'], keep='first').reset_index(drop=True)
    df_sj_max_position.to_csv(f'data/processed/{args.version}/sj_maxp.emtab2836.{args.version}.tsv.gz', index=None, sep='\t', compression='gzip')
    df_sj_max_tissue = df_sj.sort_values(by=['start', 'end', 'unique_reads'], ascending=[True, True, False]).drop_duplicates(subset=['start', 'end', 'tissue'], keep='first').reset_index(drop=True)
    df_sj_max_tissue.to_csv(f'data/processed/{args.version}/sj_maxt.emtab2836.{args.version}.tsv.gz', index=None, sep='\t', compression='gzip')

    # Loading previously annotated introns
    df_introns_annotated = pd.read_csv(f'data/interim/{args.version}/gencode.v{args.version[1:]}.annotation.complete.tsv.gz', sep='\t', compression='gzip')
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
    df_qsj.to_csv(f'data/processed/{args.version}/sj_maxp.emtab2836.{args.version}.mapped.tsv.gz', sep='\t', index=None, compression='gzip')

    # Calculating mins and exporting qsplice
    df_qsj['cds_coverage'] = df_qsj.groupby('transcript_id')['cds_coverage'].transform(lambda x: 'null' if (x.str.contains('none')).all() else x)
    df_qsj.loc[df_qsj['cds_coverage'] == 'none', 'unique_reads'] = df_qsj['unique_reads'].max()
    df_qsj = df_qsj.loc[df_qsj.groupby('transcript_id')['unique_reads'].idxmin()].sort_values(by=['seq_id', 'gene_id', 'transcript_id', 'start', 'end'], ascending=[True, True, True, True, True]).reset_index(drop=True).fillna('-')
    df_qsj = df_qsj[['seq_id','transcript_id','transcript_type','strand','gene_id','gene_name','gene_type','intron_number','start','end','nexons','ncds','unique_reads','tissue','gene_mean','gene_mean_cds','RNA2sj','norm_RNA2sj','RNA2sj_cds','norm_RNA2sj_cds']]
    df_qsj['gene_id'] = df_qsj['gene_id'].str.split('.').str[0]
    df_qsj['transcript_id'] = df_qsj['transcript_id'].str.split('.').str[0]
    df_qsj.to_csv(f'data/processed/{args.version}/qsplice.emtab2836.{args.version}.tsv.gz', sep='\t', index=None, compression='gzip')

if __name__ == "__main__":
    main()
