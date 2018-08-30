import pandas as pd
import sys
import os
import wgdi
import argparse
from wgdi.dotplot import dotplot
from wgdi.align_dotplot import align_dotplot
from wgdi.block_correspondence import block_correspondence
from wgdi.block_ks import block_ks
from wgdi.retain import retain
from wgdi.ks import ks
import wgdi.base as base

parser = argparse.ArgumentParser(prog='wgdi',usage='%(prog)s [options]',epilog="")
parser.description='''The sets of programs mainly packaged the most popular and professional methods
 to comparatively analyze the post Whole Genome Duplication, for revealing the regular of genome evolutionary. '''
parser.add_argument("-v", "--version",action='version', version='0.1.1')
parser.add_argument("-d", "--dotplot",help="""More information, try to use 'wgdi -d example'""")
parser.add_argument("-c", "--correspondence",help="More information, try to use 'wgdi -c example'")
parser.add_argument("-a", "--alignment",help="More information, try to use 'wgdi -a example'")
parser.add_argument("-r", "--retain",help="More information, try to use 'wgdi -r example'")
parser.add_argument("-bk", "--blkks",help="More information, try to use 'wgdi -bk example'")
parser.add_argument("-ks", "--calks",help="More information, try to use 'wgdi -ks example'")
args = parser.parse_args()

def run_dotplot():
    options=base.load_conf(args.dotplot,'dotplot')
    dot=dotplot(options)
    dot.run()

def run_align_dotplot():
    options=base.load_conf(args.alignment,'alignment')
    align_dot=align_dotplot(options)
    align_dot.run()

def run_align_correspondence():
    options=base.load_conf(args.correspondence,'correspondence')
    align_cor=block_correspondence(options)
    align_cor.run()

def run_retain():
    options=base.load_conf(args.retain,'retain')
    retained=retain(options)
    retained.run()

def run_block_ks():
    options=base.load_conf(args.blkks,'blkks')
    blockks=block_ks(options)
    blockks.run()

def run_cal_ks():
    options=base.load_conf(args.calks,'ks')
    calks=ks(options)
    calks.run()

def module_to_run(argument):
    switcher = {
        'dotplot' : run_dotplot,
        'correspondence' : run_align_correspondence,
        'alignment' : run_align_dotplot,
        'retain' : run_retain,
        'blkks' : run_block_ks,
        'calks' : run_cal_ks
    }
    return switcher.get(argument)()

def main():
    path=wgdi.__path__[0]
    options={'dotplot':'dotplot.conf',
            'correspondence': 'corr.conf',
            'alignment' : 'align.conf',
            'retain' : 'retain.conf',
            'blkks' : 'blkks.conf',
            'calks' : 'ks.conf'}
    for arg in vars(args):
        value=getattr(args, arg)
        if value is not None:
            if value in ['?','help','example']:
                f=open(os.path.join(path,'example\\'+options[arg]))
                print(f.read())
            elif not os.path.exists(value):
                print(value+' not exits')
                sys.exit(0)
            else:
                module_to_run(arg)