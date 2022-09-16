import argparse
import os
import shutil
import sys

import pandas as pd

import wgdi
import wgdi.base as base
from wgdi.align_dotplot import align_dotplot
from wgdi.block_correspondence import block_correspondence
from wgdi.block_info import block_info
from wgdi.block_ks import block_ks
from wgdi.circos import circos
from wgdi.dotplot import dotplot
from wgdi.karyotype import karyotype
from wgdi.karyotype_mapping import karyotype_mapping
from wgdi.ks import ks
from wgdi.ks_peaks import kspeaks
from wgdi.ksfigure import ksfigure
from wgdi.peaksfit import peaksfit
from wgdi.pindex import pindex
from wgdi.polyploidy_classification import polyploidy_classification
from wgdi.retain import retain
from wgdi.run_colliearity import mycollinearity
from wgdi.trees import trees
from wgdi.ancestral_karyotype import ancestral_karyotype
from wgdi.ancestral_karyotype_repertoire import ancestral_karyotype_repertoire

parser = argparse.ArgumentParser(
    prog='wgdi', usage='%(prog)s [options]', epilog="", formatter_class=argparse.RawDescriptionHelpFormatter,)
parser.description = '''\
WGDI(Whole-Genome Duplication Integrated analysis):  A user-friendly toolkit for evolutionary analyses of whole-genome duplications and ancestral karyotypes.

    https://wgdi.readthedocs.io/en/latest/
    -------------------------------------- '''
parser.add_argument("-v", "--version", action='version', version='0.6.1')
parser.add_argument("-d", dest="dotplot",
                    help="Show homologous gene dotplot")
parser.add_argument("-icl", dest="improvedcollinearity",
                    help="Improved version of ColinearScan ")
parser.add_argument("-ks", dest="calks",
                    help="Calculate Ka/Ks for homologous gene pairs by YN00")
parser.add_argument("-bk", dest="blockks",
                    help="Show Ks of blocks in a dotplot")
parser.add_argument("-bi", dest="blockinfo",
                    help="Collinearity and Ks speculate whole genome duplication")
parser.add_argument("-c", dest="correspondence",
                    help="Extract event-related genomic alignment")
parser.add_argument("-kp", dest="kspeaks",
                    help="A simple way to get ks peaks")
parser.add_argument("-kf", dest="ksfigure",
                    help="A simple way to draw ks distribution map")
parser.add_argument("-pf", dest="peaksfit",
                    help="Gaussian fitting of ks distribution")
parser.add_argument("-pc", dest="polyploidy_classification",
                    help="Polyploid distinguish among subgenomes")
parser.add_argument("-a", dest="alignment",
                    help="Show event-related genomic alignment in a dotplot")
parser.add_argument("-k", dest="karyotype",
                    help="Show genome evolution from reconstructed ancestors")
parser.add_argument("-ak", dest="ancestral_karyotype",
                    help="Generation of ancestral karyotypes from chromosomes that retain same structures in genomes")
parser.add_argument("-akr", dest="ancestral_karyotype_repertoire",
                    help="Incorporate genes from collinearity blocks into the ancestral karyotype repertoire")
parser.add_argument("-km", dest="karyotype_mapping",
                    help="Mapping from the known karyotype result to this species")                    
parser.add_argument("-at", dest="alignmenttrees",
                    help="Collinear genes construct phylogenetic trees")
parser.add_argument("-p", dest="pindex",
                    help="Polyploidy-index characterize the degree of divergence between subgenomes of a polyploidy")
parser.add_argument("-r", dest="retain",
                    help="Show subgenomes in gene retention or genome fractionation")
parser.add_argument("-ci", dest="circos",
                    help="A simple way to run circos")
parser.add_argument("-conf", dest="configure",
                    help="Display and modify the environment variable")
args = parser.parse_args()


def run_subprogram(program, conf, name):
    options = base.load_conf(conf, name)
    r = program(options)
    r.run()


def run_configure():
    base.rewrite(args.configure, 'ini')


def module_to_run(argument, conf):
    switcher = {
        'dotplot': (dotplot, conf, 'dotplot'),
        'correspondence': (block_correspondence, conf, 'correspondence'),
        'alignment': (align_dotplot, conf, 'alignment'),
        'retain': (retain, conf, 'retain'),
        'blockks': (block_ks, conf, 'blockks'),
        'blockinfo': (block_info, conf, 'blockinfo'),
        'calks': (ks, conf, 'ks'),
        'circos': (circos, conf, 'circos'),
        'kspeaks': (kspeaks, conf, 'kspeaks'),
        'peaksfit': (peaksfit, conf, 'peaksfit'),
        'ksfigure': (ksfigure, conf, 'ksfigure'),
        'pindex': (pindex, conf, 'pindex'),
        'alignmenttrees': (trees, conf, 'alignmenttrees'),
        'improvedcollinearity': (mycollinearity, conf, 'collinearity'),
        'configure': run_configure,
        'polyploidy_classification': (polyploidy_classification, conf, 'polyploidy classification'),
        'karyotype': (karyotype, conf, 'karyotype'),
        'ancestral_karyotype': (ancestral_karyotype, conf, 'ancestral_karyotype'),
        'karyotype_mapping': (karyotype_mapping, conf, 'karyotype_mapping'),
        'ancestral_karyotype_repertoire': (ancestral_karyotype_repertoire, conf, 'ancestral_karyotype_repertoire'),
    }
    if argument == 'configure':
        run_configure()
    else:
        program, conf, name = tuple(switcher.get(argument))
        run_subprogram(program, conf, name)


def main():
    path = wgdi.__path__[0]
    options = {'dotplot': 'dotplot.conf',
               'correspondence': 'corr.conf',
               'alignment': 'align.conf',
               'retain': 'retain.conf',
               'blockks': 'blockks.conf',
               'blockinfo': 'blockinfo.conf',
               'calks': 'ks.conf',
               'circos': 'circos.conf',
               'kspeaks': 'kspeaks.conf',
               'ksfigure': 'ksfigure.conf',
               'pindex': 'pindex.conf',
               'alignmenttrees': 'alignmenttrees.conf',
               'peaksfit': 'peaksfit.conf',
               'configure': 'conf.ini',
               'improvedcollinearity': 'collinearity.conf',
               'polyploidy_classification': 'polyploidy_classification.conf',
               'karyotype': 'karyotype.conf',
               'ancestral_karyotype': 'ancestral_karyotype.conf',
               'ancestral_karyotype_repertoire': 'ancestral_karyotype_repertoire.conf',
               'karyotype_mapping': 'karyotype_mapping.conf',
               }
    for arg in vars(args):
        value = getattr(args, arg)
        if value is not None:
            if value in ['?', 'help', 'example']:
                f = open(os.path.join(path, 'example', options[arg]))
                print(f.read())
                if arg == 'ksfigure':
                    if not os.path.exists('ks_fit_result.csv'):
                        shutil.copy2(os.path.join(
                            wgdi.__path__[0], 'example/ks_fit_result.csv'), os.getcwd())
            elif not os.path.exists(value):
                print(value+' not exits')
                sys.exit(0)
            else:
                module_to_run(arg, value)
