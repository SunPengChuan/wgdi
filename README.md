# WGDI

![Latest PyPI version](https://img.shields.io/pypi/v/wgdi.svg) [![Downloads](https://pepy.tech/badge/wgdi/month)](https://pepy.tech/project/wgdi) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/wgdi/README.html)

| | |
| --- | --- |
| Author  | Pengchuan Sun ([sunpengchuan](https//github.com/sunpengchuan)) |
| Email   | <sunpengchuan@gmail.com> |
| License | [BSD](http://creativecommons.org/licenses/BSD/) |

## Description

**WGDI (Whole-Genome Duplication Integrated analysis)** is a Python-based command-line tool designed to simplify the analysis of whole-genome duplications (WGD) and cross-species genome alignments. It offers three main workflows that enhance the detection and study of WGD events:

## Key Features

### 1. Polyploid Inference
- Identifies and confirms polyploid events with high accuracy.

### 2. Genomic Homology Inference
- Traces the evolutionary history of duplicated regions across species, with a focus on distinguishing subgenomes. 

### 3. Ancestral Karyotyping
- Reconstructs protochromosomes and traces common chromosomal rearrangements to understand chromosome evolution. 


## Installation

Python package and command line interface (IDLE) for the analysis of whole genome duplications (WGDI). WGDI can be deployed in Windows, Linux, and Mac OS operating systems and can be installed via pip and conda.

#### Bioconda

```
conda install -c bioconda  wgdi
```

#### Pypi

```
pip3 install wgdi
```

Documentation for installation along with a user tutorial, a default parameter file, and test data are provided. please consult the docs at <http://wgdi.readthedocs.io/en/latest/>.

## Tips

Here are some videos with simple examples of WGDI.

###### [WGDI的简单使用（一）](https://www.bilibili.com/video/BV1qK4y1U7eK) or https://youtu.be/k-S6FVcBIQw

###### [WGDI的简单使用（二）](https://www.bilibili.com/video/BV195411P7L1) or https://youtu.be/QiZYFYGclyE

chatting group QQ : 966612552

## Citating WGDI

If you use wgdi in your work, please cite:

> Sun P., Jiao B., Yang Y., Shan L., Li T., Li X., Xi Z., Wang X., and Liu J. (2022). WGDI: A user-friendly toolkit for evolutionary analyses of whole-genome duplications and ancestral karyotypes. Mol. Plant. doi: https://doi.org/10.1016/j.molp.2022.10.018.

## News

## 0.74
* Improved the the fusion positions dataset (-fpd).
* Fixed some issues (-pc).

## 0.7.1
* Added extract the fusion positions dataset (-fpd).
* Added determine whether these fusion events occur in other genomes (-fd).
* Improved the karyotype_mapping (-km) effect.
* Fixed the problem caused by the Python version, now it is compatible with version 3.12.


## 0.6.5
* Fixed some issues (-sf).
* Added new tips to avoid some errors.

## 0.6.4
* Fixed the problem caused by the Python version, now it is compatible with version 3.11.3.

## 0.6.3
* Fixed some issues (-ks, -sf).

## 0.6.2
* Added find shared fusions between species (-sf).

## 0.6.1

* Fixed issue with alignment (-a). Only version 0.6.0 has this bug.

## 0.6.0

* Fixed issue with improved collinearity (-icl).
* Added a parameter 'tandem_ratio' to blockinfo (-bi).

## 0.5.9

* Update the improved collinearity (-icl). Faster than before, but lower than MCscanX, JCVI.
* Fixed issue with ancestral karyotype repertoire (-akr).

## 0.5.8

* Fixed issue with gene names (-ks).

## 0.5.7
- Fixed issue with chromosome order (-ak).
- Fixed issue with gene names (-ks).  This version is not fixed, please install the latest version.

## 0.5.5 and 0.5.6
* Add ancestral karyotype (-ak)
* Add ancestral karyotype repertoire (-akr)

## 0.5.4
* Improved the karyotype_mapping (-km) effect.
* little change (-at).

## 0.5.3
* Fixed legend issue with (-kf).
* Fixed calculate Ks issue with (-ks).
* Improved the karyotype_mapping (-km) effect.
* Improved the alignmenttrees (-at) effect.

## 0.5.2
* Fixed some bugs.

## 0.5.1
* Fixed the error of the command (-conf).
* Improved the karyotype_mapping (-km) effect.
* Added the available data set of alignmenttree (-at). Low copy data set (for example, single-copy_groups.tsv of sonicparanoid2 software).

## 0.4.9
* The latest version adds karyotype_mapping (-km) and karyotype (-k) display.
* The latest version changes the calculation of extracting pvalue from collinearity (-icl), making this parameter more sensitive. Therefore, it is recommended to set to 0.2 instead of 0.05.
* The latest version has also changed the drawing display of ksfigure (-kf) to make it more beautiful.
