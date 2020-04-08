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

__author__ = "Fernando Pozo"
__copyright__ = "Copyright 2020"
__license__ = "GNU General Public License"
__version__ = "1.0.1"
__maintainer__ = "Fernando Pozo"
__email__ = "fpozoc@cnio.es"
__status__ = "Production"


def main():
    return 'a'
    # max_position = df_sjs.drop_duplicates()
    # max_tissue = df_sjs.drop_duplicates([''])

    # introns_path = '/media/hdd1/fpozoc/projects/appris_rna/data/gencode.v27.annotation.introns.tsv'
    # introns = pd.read_csv(introns_path, sep='\t', names =['chr', 'start', 'end', 'strand', 'transcript_id'])


def process_sj(inpath: str, samples_ids: {} = None) -> []:
    '''SJ.out.tab processor. Removing non mapped regions and unannotated maps.

    STAR manual (4.4) https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

    Args:
        inpath (str): SJ.out.tab infile.
        samples_ids (dictionary): Select if you want to add a tissue annotation for that SJ file. In this case introduce a dictionary with dirname identifier as keys and tissue names as values.

    Returns:
        df (list): pandas DataFrame with file processed.
    '''
    df = pd.read_csv(inpath, sep='\t', names=['chr', 'start', 'end', 'nstrand', 'intron-motif', 'annotated', 'unique_reads', 'multimapping_reads', 'overhang'])
    df = df[df['chr'].str.contains('chr')][df['annotated'] == 1][['chr', 'start', 'end', 'nstrand', 'unique_reads']]
    if samples_ids:
        df['tissue'] = samples_ids[os.path.basename(os.path.dirname(inpath))] #adding tissue
    return df


# experiment data 
def parse_emtab(inpath: str) -> {}:
    '''Parsing experiment data from:
    https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2836/samples/
    It takes a table as input and returns a dictionary with identifier as key and tissue as value.

    Args:
        inpath (str): annotation file path.

    Returns:
        samples_ids (dict): experiments identifiers with tissues as values.
    '''
    df = pd.read_csv(inpath, sep ='\t') # read table
    df['Comment[ENA_RUN]'] = df['Comment[ENA_RUN]'].apply(lambda x: x + '.1') # Adds '.1' to sample
    # df['Source Name'] = df['Source Name'].str.split('_').str[0] # Split by _
    df  = df[['Comment[ENA_RUN]', 'Source Name']].drop_duplicates(subset='Comment[ENA_RUN]').sort_values(by='Comment[ENA_RUN]').reset_index(drop=True) # Dropping duplicates names
    samples_ids = dict(zip(df['Comment[ENA_RUN]'], df['Source Name'])) # Creates a dict
    return samples_ids


def concat_samples(indir: str) -> []:
    '''Concatenate several SJ.out.tab files already processed in same pandas DataFrame.

    Args:
        indir (str): annotation file directory.

    Returns:
        df (list): pandas DataFrame with all SJ.out.tab concatenated with tissue annotations.
    '''
    experiment_pattern = 'ER*.1'
    experiment_filepath = 'E-MTAB-2836.sdrf.txt'
    glob_path = f'{indir}/{experiment_pattern}/SJ.out.tab'
    df = pd.concat([process_sj(filepath, parse_emtab(experiment_filepath)) for filepath in glob.glob(glob_path)]).reset_index(drop=True)
    return df
