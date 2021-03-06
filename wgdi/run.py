import argparse
import os
import shutil
import sys

import pandas as pd

import wgdi
import wgdi.base as base
from wgdi.align_dotplot import align_dotplot
from wgdi.block_correspondence import block_correspondence
from wgdi.polyploidy_classification import polyploidy_classification
from wgdi.block_info import block_info
from wgdi.block_ks import block_ks
from wgdi.circos import circos
from wgdi.dotplot import dotplot
from wgdi.ks import ks
from wgdi.ks_peaks import kspeaks
from wgdi.ksfigure import ksfigure
from wgdi.peaksfit import peaksfit
from wgdi.pindex import pindex
from wgdi.retain import retain
from wgdi.run_colliearity import mycollinearity
from wgdi.trees import trees

parser = argparse.ArgumentParser(
    prog='wgdi', usage='%(prog)s [options]', epilog="", formatter_class=argparse.RawDescriptionHelpFormatter,)
parser.description = '''\
This is a gold standard for complex genomic analysis,including the construction 
of homologous gene dotplot, event-related genomic alignment, and synonymous 
substitutions, and differences in different evolution rates, etc.

    https://wgdi.readthedocs.io/en/latest/
    -------------------------------------- '''
parser.add_argument("-v", "--version", action='version', version='0.4.7')
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


def run_dotplot():
    options = base.load_conf(args.dotplot, 'dotplot')
    dot = dotplot(options)
    dot.run()


def run_block_info():
    options = base.load_conf(args.blockinfo, 'blockinfo')
    blockinfo = block_info(options)
    blockinfo.run()


def run_circos():
    options = base.load_conf(args.circos, 'circos')
    cir = circos(options)
    cir.run()


def run_peaksfit():
    options = base.load_conf(args.peaksfit, 'peaksfit')
    pf = peaksfit(options)
    pf.run()


def run_align_dotplot():
    options = base.load_conf(args.alignment, 'alignment')
    align_dot = align_dotplot(options)
    align_dot.run()


def run_align_correspondence():
    options = base.load_conf(args.correspondence, 'correspondence')
    align_cor = block_correspondence(options)
    align_cor.run()


def run_retain():
    options = base.load_conf(args.retain, 'retain')
    retained = retain(options)
    retained.run()


def run_block_ks():
    options = base.load_conf(args.blockks, 'blockks')
    blockks = block_ks(options)
    blockks.run()


def run_kspeaks():
    options = base.load_conf(args.kspeaks, 'kspeaks')
    kp = kspeaks(options)
    kp.run()


def run_ksfigure():
    options = base.load_conf(args.ksfigure, 'ksfigure')
    kf = ksfigure(options)
    kf.run()

def run_polyploidy_classification():
    options = base.load_conf(args.polyploidy_classification, 'polyploidy classification')
    pc = polyploidy_classification(options)
    pc.run()


def run_pindex():
    options = base.load_conf(args.pindex, 'pindex')
    p = pindex(options)
    p.run()


def run_trees():
    options = base.load_conf(args.alignmenttrees, 'alignmenttrees')
    t = trees(options)
    t.run()


def run_cal_ks():
    options = base.load_conf(args.calks, 'ks')
    calks = ks(options)
    calks.run()


def run_collinearity():
    options = base.load_conf(args.improvedcollinearity, 'collinearity')
    col = mycollinearity(options)
    col.run()


def run_configure():
    base.rewrite(args.configure, 'ini')


def module_to_run(argument):
    switcher = {
        'dotplot': run_dotplot,
        'correspondence': run_align_correspondence,
        'alignment': run_align_dotplot,
        'retain': run_retain,
        'blockks': run_block_ks,
        'blockinfo': run_block_info,
        'calks': run_cal_ks,
        'circos': run_circos,
        'kspeaks': run_kspeaks,
        'peaksfit': run_peaksfit,
        'ksfigure': run_ksfigure,
        'pindex': run_pindex,
        'alignmenttrees': run_trees,
        'improvedcollinearity': run_collinearity,
        'configure': run_configure,
        'polyploidy_classification':run_polyploidy_classification,
    }
    return switcher.get(argument)()


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
               'polyploidy_classification':'polyploidy_classification.conf',
               }
    for arg in vars(args):
        value = getattr(args, arg)
        if value is not None:
            if value in ['?', 'help', 'example']:
                f = open(os.path.join(path, 'example', options[arg]))
                print(f.read())
                if arg == 'ksfigure':
                    if not os.path.exists('ks_fit_result.csv'):
                        shutil.copy2(os.path.join(wgdi.__path__[0], 'example/ks_fit_result.csv'), os.getcwd())
            elif not os.path.exists(value):
                print(value+' not exits')
                sys.exit(0)
            else:
                module_to_run(arg)
