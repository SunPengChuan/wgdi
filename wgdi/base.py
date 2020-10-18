import configparser
import os
import re

import numpy as np
import pandas as pd
import wgdi
from Bio import Seq, SeqIO, SeqRecord


def config():
    conf = configparser.ConfigParser()
    conf.read(os.path.join(wgdi.__path__[0], 'conf.ini'))
    return conf.items('ini')


def load_conf(file, section):
    conf = configparser.ConfigParser()
    conf.read(file)
    return conf.items(section)


def read_colinearscan(file):
    data, b, flag, num = [], [], 0, 1
    with open(file) as f:
        for line in f.readlines():
            line = line.strip()
            if re.match(r"the", line):
                num = re.search('\d+', line).group()
                b = []
                flag = 1
                continue
            if re.match(r"\>LOCALE", line):
                flag = 0
                p = re.split(':', line)
                if len(b) > 0:
                    data.append([num, b, p[1]])
                b = []
                continue
            if flag == 1:
                a = re.split(r"\s", line)
                b.append(a)
    return data

def read_mcscanx(fn):
    f1=open(fn)   
    data,b=[],[]
    flag,num=0,0
    for line in f1.readlines():
        line=line.strip()
        if re.match(r"## Alignment",line):
            flag=1
            if len(b)==0:
                arr= re.findall(r"[\d+\.]+", line)[0]
                continue
            data.append([num,b,0])
            b=[]
            num = re.findall(r"\d+", line)[0]
            continue
        if flag==0:
            continue
        a=re.split(r"\:",line)
        c=re.split(r"\s+",a[1])
        b.append([c[1],c[1],c[2],c[2]])
    data.append([num,b,0])
    return data

def read_coliearity(fn):
    f1=open(fn)   
    data,b=[],[]
    flag,num=0,0
    for line in f1.readlines():
        line=line.strip()
        if re.match(r"# Alignment",line):
            flag=1
            if len(b)==0:
                arr =re.findall('[\.\d+]+',line)
                continue
            data.append([arr[0],b,arr[2]])
            b=[]
            arr =re.findall('[\.\d+]+',line)
            continue
        if flag==0:
            continue
        b.append(re.split(r"\s", line))
    data.append([arr[0],b,arr[2]])
    return data

def read_ks(file, col):
    ks = pd.read_csv(file, sep='\t')
    ks.drop_duplicates(subset=['id1', 'id2'], keep='first', inplace=True)
    ks[col] = ks[col].astype(float)
    ks = ks[ks[col] >= 0]
    ks.index = ks['id1']+','+ks['id2']
    return ks[col]


def get_median(data):
    if len(data) == 0:
        return 0
    data.sort()
    half = len(data) // 2
    return (data[half] + data[~half]) / 2


def cds_to_pep(cds_file, pep_file, fmt='fasta'):
    records = list(SeqIO.parse(cds_file, fmt))
    for k in records:
        k.seq = k.seq.translate()
    SeqIO.write(records, pep_file, 'fasta')
    return True


def newblast(file, score, evalue, gene_loc1, gene_loc2, reverse):
    blast = pd.read_csv(file, sep="\t", header=None)
    if reverse.upper() == 'TRUE':
        blast[[0, 1]] = blast[[1, 0]]
    blast = blast[(blast[11] >= score) & (
        blast[10] < evalue) & (blast[1] != blast[0])]
    blast = blast[(blast[0].isin(gene_loc1.index)) &
                  (blast[1].isin(gene_loc2.index))]
    blast.drop_duplicates(subset=[0, 1], keep='first', inplace=True)
    blast[0] = blast[0].astype(str)
    blast[1] = blast[1].astype(str)
    return blast


def newgff(file):
    gff = pd.read_csv(file, sep="\t", header=None, index_col=1)
    gff.rename(columns={0: 'chr', 2: 'start',
                        3: 'end', 4: 'stand', 5: 'order'}, inplace=True)
    gff['chr'] = gff['chr'].astype(str)
    gff['start'] = gff['start'].astype(int)
    gff['end'] = gff['end'].astype(int)
    gff['stand'] = gff['stand'].astype(str)
    gff['order'] = gff['order'].astype(int)
    return gff


def newlens(file, position):
    lens = pd.read_csv(file, sep="\t", header=None, index_col=0)
    lens.index = lens.index.astype(str)
    if position == 'order':
        lens = lens[2]
    if position == 'end':
        lens = lens[1]
    return lens


def gene_location(gff, lens, step, position):
    loc_gene, dict_chr, n = {}, {}, 0
    gff = gff[gff['chr'].isin(lens.index)].copy()
    dict_chr = dict(zip(lens.index, np.append(
        np.array([0]), lens.cumsum()[:-1].values)))
    gff.loc[:, 'loc'] = ''
    for name, group in gff.groupby(['chr']):
        gff.loc[group.index, 'loc'] = (dict_chr[name]+group[position])*step
    return gff


def dotplot_frame(fig, ax, lens1, lens2, step1, step2, genome1_name, genome2_name, arr):
    for k in lens1.cumsum()[:-1]*step1:
        ax.axhline(y=k, alpha=0.8, color='black', lw=0.5)
    for k in lens2.cumsum()[:-1]*step2:
        ax.axvline(x=k, alpha=0.8, color='black', lw=0.5)
    align = dict(family='Arial', style='normal',
                 horizontalalignment="center", verticalalignment="center")
    align1 = dict(family='Arial', style='normal',
                  horizontalalignment="right", verticalalignment="center")
    yticks = lens1.cumsum()*step1-0.5*lens1*step1
    ax.set_yticks(yticks)
    ax.set_yticklabels(lens1.index, fontsize=12, **align1)
    xticks = lens2.cumsum()*step2-0.5*lens2*step2
    ax.set_xticks(xticks)
    # ax.set_xticks([])
    # ax.set_yticks([])
    ax.set_xticklabels(lens2.index, fontsize=12, **align)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    if arr[0] < 0:
        ax.text(-0.065, 0.5, genome1_name, weight='semibold',
                fontsize=18, rotation=90, **align)
    elif arr[0] == 1:
        ax.text(-0.065, 0.5, genome1_name, weight='semibold',
                fontsize=18, rotation=90, **align)
    else:
        ax.text(-0.05, 0.5, genome1_name, weight='semibold',
                fontsize=18, rotation=90, **align)
    if arr[1] < 0:
        ax.text(0.5, -0.065, genome2_name,
                weight='semibold', fontsize=18, **align)
    else:
        ax.text(0.5, -0.05, genome2_name,
                weight='semibold', fontsize=18, **align)


def Bezier3(plist, t):
    p0, p1, p2 = plist
    return p0*(1-t)**2+2*p1*t*(1-t)+p2*t**2


def Bezier4(plist, t):
    p0, p1, p2, p3, p4 = plist
    return p0*(1-t)**4+4*p1*t*(1-t)**3+6*p2*t**2*(1-t)**2+4*p3*(1-t)*t**3+p4*t**4
