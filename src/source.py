#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" src/source.py

This file can also be imported as a module and contains the following functions:
    * create_dir
    * Gencode

TO DO:  
    *
"""

from __future__ import absolute_import, division, print_function

import ftplib, os, re, subprocess

def create_dir(dirpath: str) -> str: 
    """mkdir -p python equivalent
    
    Arguments:
        dirpath {str} -- Path to create the new folder
    
    Returns:
        absolute_path {str} -- Absolute path of the new folder
    """    
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    absolute_path = os.path.abspath(dirpath)
    return absolute_path


class Gencode():
    
    def __init__(self, version: int, specie: str):
        self.version = version
        self.specie = specie
        self.url = f'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{specie}/release_{version}'
        self.ftpfiles = []
        self.filepath = ''
        
    def connect_ftp(self):    # show last path of file
        url_site = self.url.split('ftp://')[1].split('/')[0]
        url_path = self.url.split(url_site)[1]
        fc = ftplib.FTP(url_site)
        try:
            fc.login()
            print(f'Conexion to {self.url} done.')
        except RuntimeError:
            print('Conexion cannot be established.')
            return
        fc.cwd(url_path)
        self.ftpfiles = fc.nlst()
    
    def show_files(self):
        print(self.ftpfiles)

    def regexes(self) -> []:
        catdict = {
            'all': '.', 
            'transcripts': r'gencode.v\d\d.pc_transcripts.fa.gz',
            'translation': r'gencode.v\d\d.pc_translations.fa.gz',
            'gff3': r'gencode.v\d\d.annotation.gff3.gz',
            'gtf': r'gencode.v\d\d.annotation.gtf.gz',
            'metadata': r'metadata'
            }
        return catdict

    def download(self, outdir: str = '.', type: str = None, pattern: str = None):
        self.connect_ftp()
        create_dir(outdir)
        catdict = self.regexes()
        if type == None and pattern == None:
            print(f'Please, select which file you want to download from {self.ftpfiles}')
        filename = [p for p in self.ftpfiles if re.search(catdict[type] if type else pattern, p)]
        urlname = os.path.join(self.url, filename[0])
        try:
            subprocess.call(f'wget -P {outdir} -nc --quiet {urlname}', shell=True)
            filepath = os.path.join(outdir, os.path.basename(urlname))
            print(f'{filepath} downloaded.')
            return filepath
        except subprocess.CalledProcessError:
            print("Failed calling wget - Aborting")
            return None
