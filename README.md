# QSplice

Quantifying splice junctions coverage from SJ.out.tab released by STAR mapping it to genome positions.

### Table of Contents

**[Installation Instructions](#installation-instructions)**<br>
**[Usage](#model-reproducibility)**<br>
**[Author information and license](#author-information-and-license)**<br>
**[Release History](#release-history)**<br>
**[Contributing](#contributing)**<br>
**[License](#license)**

## Installation Instructions

Run the silent installation of Miniconda in case you don't have this software in your Linux Environment

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
```

Once you have installed Miniconda/Anaconda, create a Python environment. Clone this repository and install inside your recently created conda environment.

```sh
conda create --name qsplice python=3.7
conda activate qsplice
pip install gffpandas

git clone https://gitlab.com/fpozoc/qsplice.git
cd qsplice
python setup.py install
```

## Usage

At the moment, this package has been implemented for [GENCODE](https://www.gencodegenes.org/human/) genome annotation versions. 
User has to specify the desired version typing `g` + version number, e.g. `g27`.

First, user has to generate introns the genome annotation with `intronator.py` . Download will be performed automatically.

```sh
python -m src.intronator --version g27
```

To map splice junctions previously generated (see STAR) and calculate `QSplice` scores per transcript `map_junctions.py`.

If user wants to add a custom file:

```sh
python -m src.map_junctions --version g27 --custom --file /media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g27/SJ.out.tab.concat.gz
```

In the case the user aims to concat a set of experiments:

```sh
python -m src.map_junctions --version g33 --glob /media/hdd2/fpozoc/projects/rnaseq/out/E-MTAB-2836/GRCh38/STAR/g29/ER*.1/SJ.out.tab
```

## Reference files

```sh
$ head -n 3 SJ.out.tab

chr1    14830   14969   2       2       0       0       3       44
chr1    15039   15795   2       2       1       1       9       40
chr1    15948   16606   2       2       1       0       11      47
```

From [STAR manual](http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf) documentation:

```sh
SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format. The columns
have the following meaning:
column 1: chromosome
column 2: first base of the intron (1-based)
column 3: last base of the intron (1-based)
column 4: strand (0: undefined, 1: +, 2: -)
column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:
AT/AC, 6: GT/AT
column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
column 7: number of uniquely mapping reads crossing the junction
column 8: number of multi-mapping reads crossing the junction
column 9: maximum spliced alignment overhang
```

## Author information and license

Fernando Pozo ([@fpozoca](https://twitter.com/fpozoca) – fpozoc@cnio.es)

Distributed under the GNU General Public License. See ``LICENSE`` for more information.

## Release History

* 1.0.0.

## Contributing

1. Fork it (<https://gitlab.com/fpozoc/qsplice>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

## License

See `LICENSE` [file](LICENSE).
