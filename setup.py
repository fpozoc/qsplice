#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='QSplice',
      version='1.0.0',
      description='Quantifying splice junctions coverage from SJ.out.tab released by STAR mapping it to genome positions.',
      author='Fernando Pozo',
      author_email='fpozoc@cnio.es',
      url='https://gitlab.com/fpozoc/qsplice.git',
      download_url='https://gitlab.com/fpozoc/qsplice/tarball/0.1.0',
      license='GNU General Public License',
      install_requires=
            ['pandas', 
             'numpy',
             'loguru'],
      packages=find_packages()
     )