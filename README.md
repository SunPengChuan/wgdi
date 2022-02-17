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
## 0.5.3


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

## Already cited WGDI articles
1. Genomic analysis of Medicago ruthenica provides insights into its tolerance to abiotic stress and demographic history
2. Chromosomal‐scale genome assembly of Eleutherococcus senticosus provides insights into chromosome evolution in Araliaceae
3. The Corylus mandshurica genome provides insights into the evolution of Betulaceae genomes and hazelnut breeding
4. An ancient whole-genome duplication event and its contribution to flavor compounds in the tea plant (Camellia sinensis)
5. A tetraploidization event shaped the Aquilaria sinensis genome and contributed to the ability of sesquiterpenes synthesis
6. High-quality genome assembly of Cinnamomum burmami (chvar. Borneol) provides insights into the natural borneol biosynthesis
7. The genome sequence provides insights into salt tolerance of Achnatherum splendens (Gramineae), a constructive species of alkaline grassland
8. Chromosome-level assembly of the common vetch reference genome (Vicia sativa)
9. A chromosome-level genome assembly of an alpine plant Crucihimalaya lasiocarpa provides insights into high-altitude adaptation
10. Chromosome-scale genome assembly of the diploid oat Avena longiglumis reveals the landscape of repetitive sequences, genes and chromosome evolution in grasses
11. The chromosome-level rambutan genome reveals a significant role of segmental duplication in the expansion of resistance genes
12. A chromosome-level genome assembly for the tertiary relict plant Tetracentron sinense oliv. (trochodendraceae)
13. Multi-omics reveal differentiation and maintenance of dimorphic flowers in an alpine plant on the Qinghai-Tibet Plateau
14. The Chromosome-Level Genome of Miracle Fruit (Synsepalum dulcificum) Provides New Insights Into the Evolution and Function of Miraculin.
15. A chromosome-level reference genome of Ensete glaucum gives insight into diversity, chromosomal and repetitive sequence evolution in the Musaceae
16. High-quality genome assembly of Cinnamomum burmannii (chvar. Borneol) provides insights into the natural borneol biosynthesis
17. The Chloranthus sessilifolius genome provides insight into early diversification of angiosperms