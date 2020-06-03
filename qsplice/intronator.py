#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" src/intronator.py

Usage: python -m src.intronator --version g27
___
--help      | -h    Display documentation
--version   | -v    Genome annotation version. GENCODE: `g` + `nversion`
--file      | -f    Custom or not GENCODE gff file

This file can also be imported as a module and contains the following functions:
    * generate_introns
    * annotate_introns

TO DO:
    *
"""

from __future__ import absolute_import, division, print_function

import argparse, os, re, subprocess
from loguru import logger
import pandas as pd
import numpy as np
from genannpy.annotation import GFF

__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Production"


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-v', '--version', type=str, 
                        help='Genome annotation version. GENCODE: `g` + `nversion`')
    parser.add_argument('-f', '--file', type=str, 
                        help='Custom or not GENCODE (because it could be automatically downloaded) gff file.')
    args = parser.parse_args()
    
    logger.info(f"Program has been launched succesfully.")

    raw = os.path.join(os.path.dirname(__file__), '../data/raw' , f'{args.version}')
    interim = os.path.join(os.path.dirname(__file__), '../data/interim' , f'{args.version}')

    if args.version.lower().startswith('g'):
        filepath = source.Gencode(version=int(args.version[1:]), specie='human').download(outdir=raw, type='gff3')
    else:
        filepath = args.file

    gff_path = generate_introns(filepath, outdir=interim)
    logger.info(f'Introns generated.')

    df_annotations = load_annotations(gff_path, db=args.version[0], outdir=interim)
    logger.info(f'Annotations generated.')

    df_introns_annotated = annotate_introns(df_annotations)
    logger.info(f'CDS coverage and complete database generated.')

    complete_annotation_path = gff_path.replace('gz', 'complete.tsv.gz')
    df_introns_annotated.to_csv(complete_annotation_path, index=None, sep='\t', compression='gzip')

    introns_annotation_path = gff_path.replace('gz', 'introns.tsv.gz')
    df_introns = df_introns_annotated[df_introns_annotated['type'] == 'intron'].drop(['exon_id','exon_number','start_cds','end_cds'], axis=1)
    df_introns.to_csv(introns_annotation_path, index=None, sep='\t', compression='gzip')
    logger.info(f'Introns annotation generated: {df_introns.shape[0]}')
    logger.info(f'Introns coding coverage: {df_introns.cds_coverage.value_counts().to_dict()}')

    exons_annotation_path = gff_path.replace('gz', 'exons_cds.tsv.gz')
    df_exons = df_introns_annotated[~df_introns_annotated['type'].str.contains('intron')].drop(['intron_number'], axis=1)
    df_exons.to_csv(exons_annotation_path, index=None, sep='\t', compression='gzip')
    logger.info(f'Exons annotation generated: {df_introns.shape[0]}')
    logger.info(f'Exons coding coverage {df_exons.cds_coverage.value_counts().to_dict()}')


def load_annotations(filepath:str, db:str, outdir:str):
    df = GFF(filepath, db=db).load
    df = df[df['type'].str.contains('exon|intron|CDS')]
    if db.lower().startswith('r'): # RefSeq annotations
        df['exon_id'] = df['ID'].str.extract(r'^(?:exon|CDS)-(.*?)$')
        df['exon_number'] = df['exon_id'].str.split('-').str[1]
    df = df[['seqname', 'type', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type', 'transcript_id', 'exon_id', 'exon_number']]
    logger.info(f"{df.shape[0]} annotations loaded and stored with this columns:\n{';'.join(df.columns)}")
    return df


def annotate_introns(df: list) -> list:
    """Annotating introns cds positions merging exons and cds info.

    Arguments:
        df {list} -- pandas DataFrame to input from whole introns annotation.

    Returns:
        list -- pandas DataFrame to finally output.
    """
    df_whole = _add_strand_info(df)
    df_exons = _get_exons_cds(df_whole)
    df_introns = df_whole[df_whole['type'] == 'intron'].reset_index(drop=True)

    df = pd.concat([df_exons, df_introns], axis=0).reset_index(drop=True).sort_values(by=['gene_id', 'transcript_id', 'start', 'end'], ascending=[False, True, True, True]).reset_index(drop=True)
    df = df[df['gene_type'] == 'protein_coding']
    df['nexons'] = df.groupby('transcript_id')['type'].transform(lambda x: x.str.contains('exon').sum())
    df['ncds'] = df.groupby('transcript_id')['cds_coverage'].transform(lambda x: x.str.contains('full|partial').sum())
    df.loc[df['cds_coverage'].isnull(), 'l_intron_coverage'] = df.groupby('transcript_id')['cds_coverage'].shift(1)
    df.loc[df['cds_coverage'].isnull(), 'r_intron_coverage'] = df.groupby('transcript_id')['cds_coverage'].shift(-1)
    df.loc[(df['cds_coverage'].isnull()) & (df['l_intron_coverage'].str.contains('none') | df['r_intron_coverage'].str.contains('none')), 'cds_coverage'] = 'none'
    df.loc[(df['cds_coverage'].isnull()) & (df['l_intron_coverage'].str.contains('full|partial') | df['r_intron_coverage'].str.contains('full|partial')), 'cds_coverage'] = 'full'
    df = df.drop(['l_intron_coverage', 'r_intron_coverage'], axis=1)
    return df


def generate_introns(inpath: str, outdir: str) -> str:
    """This function generates introns file from gff*.gz file path.

    Arguments:
        inpath {str} -- The name of the annotation file to process.

    Keyword Arguments:
        outdir {str} -- Introns directory (default: {'.'})

    Returns:
        str -- Path where file has been stored.
    """
    _create_dir(outdir)
    outpath = os.path.abspath(os.path.join(outdir, os.path.basename(re.sub(r'gff', r'introns.gff', inpath))))  # used re if user want to add a non gz file as input
    cmd_gt = f"zcat {inpath} | gt gff3 -tidy -retainids -addintrons | gzip > {outpath}" # tidy option described here http://genometools.org/pipermail/gt-users/2015-August/000794.html
    subprocess.call(cmd_gt, shell=True)
    return outpath

def _add_strand_info(df: list) -> list:
    """Adding intron number and correct positions for CDS. (corrected by genome strand)

    Arguments:
        df {list} -- pandas DataFrame to input from initial introns annotation.

    Returns:
        list -- pandas DataFrame to output with introns info added.
    """
    df['intron_number'] = df['exon_number']
    intronf_list = ['gene_name','gene_id','gene_type','transcript_id']
    for feature in intronf_list:
        df.loc[df['type'].str.contains('exon|intron'), feature] = df.fillna(method='ffill')
    df.loc[(df['strand']=='+') & (df['type'].str.contains('exon|intron')), 'intron_number'] = df.fillna(method='ffill')
    df.loc[(df['strand']=='-') & (df['type'].str.contains('exon|intron')), 'intron_number'] = df.fillna(method='bfill')
    df.loc[~df['type'].str.contains('intron'), 'intron_number'] = '-'
    df.loc[(df['strand']=='+') & (df['type'].str.contains('exon|CDS')), 'exon_number'] = df.fillna(method='ffill')
    df.loc[(df['strand']=='-') & (df['type'].str.contains('exon|CDS')), 'exon_number'] = df.fillna(method='bfill')
    df.loc[(df['strand']=='+') & (df['type'].str.contains('exon|CDS')), 'exon_id'] = df.fillna(method='ffill')
    df.loc[(df['strand']=='-') & (df['type'].str.contains('exon|CDS')), 'exon_id'] = df.fillna(method='bfill')
    return df


def _create_dir(dirpath: str) -> str: 
    """mkdir -p python equivalent
    
    Arguments:
        dirpath {str} -- Path to create the new folder
    
    Returns:
        absolute_path {str} -- Absolute path of the new folder
    """    
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    absolute_path = os.path.abspath(dirpath)


def _get_exons_cds(df: list) -> list:
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
    df.insert(1, 'type', 'exon')
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
    df = df.drop_duplicates(['start', 'end', 'seqname', 'gene_id', 'transcript_id'])
    return df


if __name__ == "__main__":
    main()
