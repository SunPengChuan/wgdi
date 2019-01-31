import argparse
import os
import sys

import pandas as pd
import wgdi
import wgdi.base as base
from wgdi.align_dotplot import align_dotplot
from wgdi.block_correspondence import block_correspondence
from wgdi.block_ks import block_ks
from wgdi.colinearscan import colinearscan
from wgdi.dotplot import dotplot
from wgdi.ks import ks
from wgdi.retain import retain

parser = argparse.ArgumentParser(
    prog='wgdi', usage='%(prog)s [options]', epilog="", formatter_class=argparse.RawDescriptionHelpFormatter,)
parser.description = '''\
This is a gold standard for complex genomic analysis,including the construction 
of homologous gene dotplot, event-related genomic alignment, and synonymous 
substitutions, and differences in different evolution rates, etc.

    https://wgdi.readthedocs.io/en/latest/
    -------------------------------------- '''
parser.add_argument("-v", "--version", action='version', version='0.1.8')
parser.add_argument("-d", dest="dotplot",
                    help="Show homologous gene dotplot")
parser.add_argument("-c", dest="correspondence",
                    help="Extract event-related genomic alignment")
parser.add_argument("-a", dest="alignment",
                    help="Show event-related genomic alignment in a dotplot")
parser.add_argument("-r", dest="retain",
                    help="Show subgenomes in gene retention or genome fractionation")
parser.add_argument("-bk", dest="blockks",
                    help="Show Ks of blocks in a dotplot")
parser.add_argument("-ks", dest="calks",
                    help="Calculate Ka/Ks for homologous gene pairs by Comdel")
parser.add_argument("-cl", dest="collinearity",
                    help="A simple way to run ColinearScan")

args = parser.parse_args()


def run_dotplot():
    options = base.load_conf(args.dotplot, 'dotplot')
    dot = dotplot(options)
    dot.run()


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
    options = base.load_conf(args.blkks, 'blkks')
    blockks = block_ks(options)
    blockks.run()


def run_cal_ks():
    options = base.load_conf(args.calks, 'ks')
    calks = ks(options)
    calks.run()


def run_colinearscan():
    options = base.load_conf(args.collinearity, 'colinearscan')
    col = colinearscan(options)
    col.run()


def module_to_run(argument):
    switcher = {
        'dotplot': run_dotplot,
        'correspondence': run_align_correspondence,
        'alignment': run_align_dotplot,
        'retain': run_retain,
        'blkks': run_block_ks,
        'calks': run_cal_ks,
        'collinearity': run_colinearscan
    }
    return switcher.get(argument)()


def main():
    path = wgdi.__path__[0]
    options = {'dotplot': 'dotplot.conf',
               'correspondence': 'corr.conf',
               'alignment': 'align.conf',
               'retain': 'retain.conf',
               'blkks': 'blkks.conf',
               'calks': 'ks.conf',
               'collinearity': 'colinearscan.conf'}
    for arg in vars(args):
        value = getattr(args, arg)
        if value is not None:
            if value in ['?', 'help', 'example']:
                f = open(os.path.join(path, 'example', options[arg]))
                print(f.read())
            elif not os.path.exists(value):
                print(value+' not exits')
                sys.exit(0)
            else:
                module_to_run(arg)
