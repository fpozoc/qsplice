#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" generate_introns.py

Usage: python intronator.py --gencode 27
___
--help | -h Display documentation

DESCRIPTION 

This file can also be imported as a module and contains the following functions:
    * create_dir - returns the path of the directory created.
    * generate_introns - returns the path of the introns file.

TO DO:  
    *
"""

from __future__ import absolute_import, division, print_function

import argparse, os, re, subprocess
import source

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

    if args.gencode:
        path = source.Gencode(version=args.gencode, specie='human').download(outdir='../data/genomes', type='gff3')
    else:
        path = args.file

    generate_introns(path, intronsdir='../data/introns')


def create_dir(dirpath: str) -> str: 
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    return os.path.abspath(dirpath)


def generate_introns(inpath: str, intronsdir: str = '.') -> str:
    '''This function generates introns file from gff*.gz file path.

    Args:
       inpath (str):  The name of the annotation file to process.

    Returns:
        outpath (str): Path where file has been stored.
    '''
    subprocess.call('which gt || apt install gt', shell=True)
    subprocess.call('which rg || apt install rg', shell=True)

    create_dir(intronsdir)
    intronsname = os.path.basename(re.sub(r'gff.*', 'introns.gff', inpath))
    outpath = os.path.abspath(os.path.join(intronsdir, intronsname))

    # cmd_awk = r"awk -F'[\t=]' '{print $1,$4,$5,$7,$10}' OFS='\t'"
    # cmd_gt = f"zcat {inpath} | rg '(#|gene_type=protein_coding)' | gt gff3 -retainids -addintrons | rg 'intron' | {cmd_awk} > {outpath}"
    cmd_gt = f"zcat {inpath} | gt gff3 -tidy -retainids -addintrons > {outpath}" # tidy option described here http://genometools.org/pipermail/gt-users/2015-August/000794.html
    subprocess.call(cmd_gt, shell=True)

    print(f'{outpath} has been generated.')
    return outpath
        

if __name__ == "__main__":
    main()
