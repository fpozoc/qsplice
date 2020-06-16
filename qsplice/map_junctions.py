#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" qsplice/map_junctions.py

Usage: 
python -m src.map_junctions --version g27 --custom --file /media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g27/SJ.out.tab.concat.gz
python -m src.map_junctions --version g33 --glob /media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g29
___
--help      | -h    Display documentation
--version   | -v    Genome annotation version. GENCODE: `g` + `nversion`
--file      | -f    Customized SJ file
--glob      | -g    Directory which contains the files to be globbed and concatenated

TO DO:
    *
"""

from __future__ import absolute_import, division, print_function

import argparse, glob, os
import pandas as pd
from loguru import logger

__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Production"


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-v', '--version', type=str, help='Genome annotation version. GENCODE: `g` + `nversion`')
    parser.add_argument('-g', '--globdir', type=str, help='Directory which contains the files to be globbed and concatenated', 
                        default='/media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g29')
    parser.add_argument('-c', '--custom', help='Customized SJ file.', action='store_true', default=False)
    parser.add_argument('-f', '--file', type=str, help='Custom splice junctions file (in gzip)',
                        default='/media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g27/SJ.out.tab.concat.gz')
    args = parser.parse_args()

    logger.info(f"Program has been launched succesfully.")

    processed = os.path.join(os.path.dirname(__file__), '../data/processed' , f'{args.version}')
    interim = os.path.join(os.path.dirname(__file__), '../data/interim' , f'{args.version}')
    _create_dir(processed)

    if args.custom:
        # cat /media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g27/ERR*/SJ.out.tab > ./SJ.out.tab.concat && gzip SJ.out.concat 
        df_sj = pd.read_csv(args.file, compression='gzip', sep='\t', 
                            names=['seqname', 'start', 'end', 'nstrand', 'unique_reads', 'tissue'])
        df_sj = df_sj.drop('nstrand', axis=1)
    else:
        df_sj = concat_samples(f'{args.globdir}/*/SJ.out.tab') # concatenating several SJ.out after being annotated and parsed in same pandas DataFrame
        df_sj['tissue'] = df_sj['tissue'].str.split('_').str[0]  # Getting max values per position and per tissue group sample

    df_sj_max_position, df_sj_max_tissue =  map_junctions_positions(df_sj)

    sj_maxp_path = os.path.join(processed, f'sj_maxp.emtab2836.{args.version}.tsv.gz')
    df_sj_max_position.to_csv(sj_maxp_path, index=None, sep='\t', compression='gzip')

    sj_maxt_path = os.path.join(processed, f'sj_maxt.emtab2836.{args.version}.tsv.gz')
    df_sj_max_tissue.to_csv(sj_maxt_path, index=None, sep='\t', compression='gzip')

    df_introns = pd.read_csv(glob.glob(os.path.join(interim, '*introns.tsv.gz'))[0], sep='\t', compression='gzip')

    df_junction_score = score_per_junction(df_introns, df_sj_max_position)
    df_junction_score_path = os.path.join(processed, f'sj_maxp.emtab2836.{args.version}.mapped.tsv.gz')
    df_junction_score.to_csv(df_junction_score_path, sep='\t', index=None, compression='gzip')

    df_qsplice = score_per_transcript(df_junction_score)
    qsplice_path = os.path.join(processed, f'qsplice.emtab2836.{args.version}.tsv.gz')
    df_qsplice.to_csv(qsplice_path, sep='\t', index=None, compression='gzip')
    

def concat_samples(indir: str) -> []:
    '''Concatenate several SJ.out.tab files already processed in same pandas DataFrame.

    Args:
        indir (str): annotation file directory.

    Returns:
        df (list): pandas DataFrame with all SJ.out.tab concatenated with tissue annotations.
    '''
    globbed_dir = glob.glob(f'{indir}')
    annotation_dict = _parse_emtab(os.path.join(os.path.dirname(__file__), '../data/external/E-MTAB-2836.sdrf.txt'))
    df = pd.concat([_process_sj(filepath, annotation_dict) for filepath in globbed_dir]).reset_index(drop=True)
    logger.info(f'{len(globbed_dir)} RNA-seq samples loaded and concatenated.')
    return df


def map_junctions_positions(df:list):
    """This function use a pandas DataFrame with concatenated RNA-seq (SJ.out.tab)
    concatenated to extract the highest coverage per position and per tissue and 
    position.

    Arguments:
        df {list} -- pandas DataFrame splice junctions and reads concatenated.

    Returns:
        df_sj_max_position {list} -- pandas DataFrame with maximum coverage per junction position.
        df_sj_max_tissue {list} -- pandas DataFrame with maximum coverage per junction position and tissue.
    """    
    df_sj_max_position = df.sort_values(
        by=['start', 'end', 'unique_reads'], 
        ascending=[True, True, False]).drop_duplicates(subset=['start', 'end'], keep='first').reset_index(drop=True)
    logger.info(f"Mean unique reads per position: {df_sj_max_position['unique_reads'].mean().round(3)}")

    df_sj_max_tissue = df.sort_values(
        by=['start', 'end', 'unique_reads'], 
        ascending=[True, True, False]).drop_duplicates(subset=['start', 'end', 'tissue'], keep='first').reset_index(drop=True)
    logger.info(f"Mean unique reads per tissue and position: {df_sj_max_tissue['unique_reads'].mean().round(3)}")
    return df_sj_max_position, df_sj_max_tissue


def score_per_junction(df_introns:list, df_sj_max_position:list)->list:
    """Calculating means and score per junction.

    Arguments:
        df_introns {list} -- pandas DataFrame with introns.
        df_sj_max_postion {list} -- pandas DataFrame with junction read positions.

    Returns:
        list -- pandas DataFrame with score per exon, gene and transcript.
    """    
    df = pd.merge(df_introns, df_sj_max_position, how='left', on=['seqname', 'start', 'end'])
    df.loc[df['unique_reads'].isnull(), 'unique_reads'] = 0
    df.loc[df['unique_reads'].isnull(), 'tissue'] = '-'
    df['gene_mean'] = df.groupby('gene_id')['unique_reads'].transform(lambda x:x.mean())
    df['gene_mean'] = df['gene_mean'].fillna(df.groupby('gene_id')['gene_mean'].transform('mean'))
    df['gene_mean_cds'] = df.loc[df['cds_coverage'] == 'full'].groupby('gene_id')['unique_reads'].transform(lambda x:x.mean())
    df['gene_mean_cds'] = df['gene_mean_cds'].fillna(df.groupby('gene_id')['gene_mean_cds'].transform('mean'))
    df['RNA2sj'] = df['unique_reads']/df['gene_mean']
    df['RNA2sj_cds'] = df['unique_reads']/df['gene_mean_cds']
    df['norm_RNA2sj'] = df.groupby(['gene_id'])['RNA2sj'].transform(lambda x: (x-x.min()) / (x.max()-x.min()))
    df['norm_RNA2sj_cds'] = df.groupby(['gene_id'])['RNA2sj_cds'].transform(lambda x: (x-x.min()) / (x.max()-x.min()))
    df[df.columns[-6:]] = df[df.columns[-6:]].round(4)
    df = df.sort_values(by=['seqname', 'gene_id', 'transcript_id', 'start', 'end'])
    logger.info(f'{df.shape[0]} junctions from protein-coding genes loaded and quantified.')
    return df


def score_per_transcript(df:list)->list:
    """Calculating mins and exporting qsplice

    Arguments:
        df {list} -- pandas DataFrame with score per junction.

    Returns:
        list -- pandas DataFrame with scores per transcript (qsplice output).
    """    
    df['cds_coverage'] = df.groupby('transcript_id')['cds_coverage'].transform(
        lambda x: 'null' if (x.str.contains('none')).all() else x)
    df.loc[df['cds_coverage'] == 'none', 'unique_reads'] = df['unique_reads'].max()
    df = df.loc[df.groupby('transcript_id')['unique_reads'].idxmin()].sort_values(
        by=['seqname', 'gene_id', 'transcript_id', 'start', 'end'], 
        ascending=[True, True, True, True, True]).reset_index(drop=True)
    df = df.drop(['type', 'cds_coverage', 'start', 'end', 'strand', 'intron_number'], axis=1)
    df = df.sort_values(by=['seqname', 'gene_id', 'transcript_id', 'norm_RNA2sj_cds'], ascending=[True, True, True, False])
    logger.info(f'{df.shape[0]} transcripts from protein-coding genes loaded and quantified.')
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


def _parse_emtab(inpath: str) -> {}:
    '''Parsing experiment data from:
    https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2836/samples/
    It takes a table as input and returns a dictionary with identifier as key and tissue as value.

    Args:
        inpath (str): annotation file path.

    Returns:
        samples_ids (dict): experiments identifiers with tissues as values.
    '''
    df = pd.read_csv(inpath, sep ='\t')
    df['Comment[ENA_RUN]'] = df['Comment[ENA_RUN]'].apply(lambda x: x + '.1') # Adds '.1' to sample
    df  = df[['Comment[ENA_RUN]', 'Source Name']].drop_duplicates(subset='Comment[ENA_RUN]').sort_values(by='Comment[ENA_RUN]').reset_index(drop=True) # Dropping duplicates names
    samples_ids = dict(zip(df['Comment[ENA_RUN]'], df['Source Name']))
    return samples_ids


def _process_sj(inpath: str, samples_ids: {} = None) -> []:    
    '''SJ.out.tab processor. Removing non mapped regions and unannotated maps.

    STAR manual (4.4) https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

    Args:
        inpath (str): SJ.out.tab infile.
        samples_ids (dictionary): Select if you want to add a tissue annotation for that SJ file. In this case introduce a dictionary with dirname identifier as keys and tissue names as values.

    Returns:
        df (list): pandas DataFrame with file processed.
    '''
    df = pd.read_csv(inpath, sep='\t', names=['seqname', 'start', 'end', 'nstrand', 'intron-motif', 'annotated', 'unique_reads', 'multimapping_reads', 'overhang'])
    # df = df[df['seq_id'].str.contains('chr')][df['annotated'] == 1][['seq_id', 'start', 'end', 'nstrand', 'unique_reads']]
    df = df[df['annotated'] == 1][['seqname', 'start', 'end', 'nstrand', 'unique_reads']]
    if samples_ids:
        df['tissue'] = samples_ids[os.path.basename(os.path.dirname(inpath))] #adding tissue
    df = df.drop('nstrand', axis=1)
    return df


if __name__ == "__main__":
    main()
