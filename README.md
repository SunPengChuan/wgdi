# WGDI

![Latest PyPI version](https://img.shields.io/pypi/v/wgdi.svg) [![Downloads](https://pepy.tech/badge/wgdi/month)](https://pepy.tech/project/wgdi)

| | |
| --- | --- |
| Author  | Pengchuan Sun ([sunpengchuan](https//github.com/sunpengchuan)) |
| Email   | <sunpengchuan@gmail.com> |
| License | [BSD](http://creativecommons.org/licenses/BSD/) |

## Description

WGDI (Whole-Genome Duplication Integrated analysis), a Python-based command-line tool that facilitates comprehensive analysis of recursive polyploidizations and cross-species genome alignments.

WGDI supports three main workflows (polyploid inference, hierarchical inference of genomic homology, and ancestral chromosomal karyotyping) that can improve detection of WGD and characterization of related events. It incorporates a more sensitive and accurate collinearity detection algorithm than previous softwares, and can accelerate WGD-related karyotype research.

WGDI outperforms similar tools in terms of efficiency, flexibility and scalability.

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

###### [WGDI的简单使用（一）](https://www.bilibili.com/video/BV1qK4y1U7eK)

###### [WGDI的简单使用（二）](https://www.bilibili.com/video/BV195411P7L1)

## Citating WGDI

If you use wgdi in your work, please cite:

> Sun P, Jiao B, Yang Y, et al. WGDI: A user-friendly toolkit for evolutionary analyses of whole-genome duplications and ancestral karyotypes[J]. bioRxiv, 2021. **doi:** https://doi.org/10.1101/2021.04.29.441969

## News

* The latest version adds karyotype_mapping (-km) and karyotype (-k) display.
* The latest version changes the calculation of extracting pvalue from collinearity (-icl), making this parameter more sensitive. Therefore, it is recommended to set to 0.2 instead of 0.05.
* The latest version has also changed the drawing display of ksfigure (-kf) to make it more beautiful.


