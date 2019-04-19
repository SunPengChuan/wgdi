import configparser
import os
import re

import pandas as pd
import numpy as np
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
    data, b, flag, num= [], [], 0, 1
    with open(file) as f:
        for line in f.readlines():
            line = line.strip()
            if re.match(r"the", line):
                num = re.search('\d+',line).group()
                b = []
                flag = 1
                continue
            if re.match(r"\>LOCALE", line):
                flag = 0
                p = re.split(':', line)
                if len(b)>0:
                    data.append([num,b,p[1]])
                b = []
                continue
            if flag == 1:
                a = re.split(r"\s", line)
                b.append(a)
    return data


def read_mcscanx(fn):
    f1 = open(fn)
    data, b = [], []
    flag = 0
    for line in f1.readlines():
        line = line.strip()
        if re.match(r"## Alignment", line):
            flag = 1
            if len(b) == 0:
                b.append([line])
                continue
            data.append(b)
            b = []
            b.append([line])
            continue
        if flag == 0:
            continue
        if re.match(r'#', line):
            continue
        a = re.split(r"\:", line)
        c = re.split(r"\s+", a[1])
        b.append([c[1], c[2]])
    data.append(b)
    return data


def read_ks(file):
    ks = pd.read_csv(file, sep='\t', header=None)
    ks = ks.drop_duplicates()
    ks = ks[ks[3] > 0]
    ks.index = ks[0]+','+ks[1]
    return ks


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


def tandem(chr1, chr2, loc1, loc2):
    if (chr1 == chr2) and (abs(float(loc1)-float(loc2)) < 200):
        return True
    return False


def newblast(file, score, evalue, gene_loc1, gene_loc2):
    blast = pd.read_csv(file, sep="\t", header=None)
    blast = blast[(blast[11] >= score) & (
        blast[10] < evalue) & (blast[1] != blast[0])]
    blast = blast[(blast[0].isin(gene_loc1.index)) & (blast[1].isin(gene_loc2.index))]
    blast.drop_duplicates(subset=[0, 1], keep='first', inplace=True)
    blast[0] = blast[0].astype(str)
    blast[1] = blast[1].astype(str)
    return blast


def newgff(file):
    gff = pd.read_csv(file, sep="\t", header=None,index_col=1)
    gff.rename(columns={0: 'chr', 2: 'start',
                        3: 'end', 4: 'stand', 5: 'order'}, inplace=True)
    gff['chr'] = gff['chr'].astype(str)
    gff['start'] = gff['start'].astype(float)
    gff['end'] = gff['end'].astype(float)
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
    gff = gff[gff['chr'].isin(lens.index)]
    dict_chr=dict(zip(lens.index, np.append(np.array([0]),lens.cumsum()[:-1].values)))
    gff.loc[:,'loc']=''
    for name, group in gff.groupby(['chr']):
        gff.loc[group.index,'loc'] = (dict_chr[name]+group[position])*step
    return gff


def dotplot_frame(fig, ax, lens1, lens2, step1, step2, genome1_name, genome2_name):
    for k in lens1.cumsum()[:-1]*step1:
        ax.axhline(y=k, alpha=0.8, color='black', lw=0.5)
    for k in lens2.cumsum()[:-1]*step2:
        ax.axvline(x=k, alpha=0.8, color='black', lw=0.5)
    align = dict(family='Times New Roman', style='normal',
                 horizontalalignment="center", verticalalignment="center")
    yticks = lens1.cumsum()*step1-0.5*lens1*step1
    ax.set_yticks(yticks)
    ax.set_yticklabels(lens1.index, fontsize=12, **align)
    xticks = lens2.cumsum()*step2-0.5*lens2*step2
    ax.set_xticks(xticks)
    ax.set_xticklabels(lens2.index, fontsize=12, **align)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.axis([0, 1, 1, 0])
    ax.set_ylabel(genome1_name, labelpad=8,
                  weight='semibold', fontsize=18, **align)
    fig.suptitle(genome2_name, weight='semibold', fontsize=18, **align)

# if __name__ == "__main__":
#     config()
