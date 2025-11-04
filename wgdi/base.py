import configparser
import hashlib
import os
import re

import matplotlib
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from Bio import SeqIO

import wgdi


def gen_md5_id(item):
    """Generate MD5 hash for the given item."""
    return hashlib.md5(item.encode('utf-8')).hexdigest()


def config():
    """Read configuration from the example conf.ini file."""
    conf = configparser.ConfigParser()
    conf.read(os.path.join(wgdi.__path__[0], 'example/conf.ini'))
    return conf.items('ini')


def load_conf(file, section):
    """Load configuration items from the specified section."""
    conf = configparser.ConfigParser()
    conf.read(file)
    return conf.items(section)


def rewrite(file, section):
    """Rewrite the configuration file to keep only the specified section."""
    conf = configparser.ConfigParser()
    conf.read(file)
    if conf.has_section(section):
        for k in conf.sections():
            if k != section:
                conf.remove_section(k)
        conf.write(open(os.path.join(wgdi.__path__[0], 'example/conf.ini'), 'w'))
        print('Option ini has been modified')
    else:
        print('Option ini no change')


def read_colinearscan(file):
    """Read colinearscan output and parse into data structure."""
    data, b, flag, num = [], [], 0, 1
    with open(file) as f:
        for line in f:
            line = line.strip()
            if re.match(r"the", line):
                num = re.search(r'\d+', line).group()
                b = []
                flag = 1
                continue
            if re.match(r"\>LOCALE", line):
                flag = 0
                p = re.split(':', line)
                if b:
                    data.append([num, b, p[1]])
                b = []
                continue
            if flag == 1:
                a = re.split(r"\s", line)
                b.append(a)
    if b:
        data.append([num, b, p[1]])
    return data


def read_mcscanx(fn):
    """Read mcscanx output and parse into data structure."""
    with open(fn) as f1:
        data, b = [], []
        flag, num = 0, 0
        for line in f1:
            line = line.strip()
            if re.match(r"## Alignment", line):
                flag = 1
                if not b:
                    arr = re.findall(r"[\d+\.]+", line)[0]
                    continue
                data.append([num, b, 0])
                b = []
                num = re.findall(r"\d+", line)[0]
                continue
            if flag == 0:
                continue
            a = re.split(r"\:", line)
            c = re.split(r"\s+", a[1])
            b.append([c[1], c[1], c[2], c[2]])
        if b:
            data.append([num, b, 0])
    return data


def read_jcvi(fn):
    """Read jcvi output and parse into data structure."""
    with open(fn) as f1:
        data, b = [], []
        num = 1
        for line in f1:
            line = line.strip()
            if re.match(r"###", line):
                if b:
                    data.append([num, b, 0])
                    b = []
                num += 1
                continue
            a = re.split(r"\t", line)
            b.append([a[0], a[0], a[1], a[1]])
        if b:
            data.append([num, b, 0])
    return data


def read_collinearity(fn):
    """Read collinearity output and parse into data structure."""
    with open(fn) as f1:
        data, b = [], []
        flag, arr = 0, []
        for line in f1:
            line = line.strip()
            if re.match(r"# Alignment", line):
                flag = 1
                if not b:
                    arr = re.findall(r'[\.\d+]+', line)
                    continue
                data.append([arr[0], b, arr[2]])
                b = []
                arr = re.findall(r'[\.\d+]+', line)
                continue
            if flag == 0:
                continue
            b.append(re.split(r"\s", line))
        if b:
            data.append([arr[0], b, arr[2]])
    return data


def read_ks(file, col):
    """Read KS values from file and select specified column."""
    ks = pd.read_csv(file, sep='\t')
    ks.drop_duplicates(subset=['id1', 'id2'], keep='first', inplace=True)
    ks[col] = ks[col].astype(float)
    ks = ks[ks[col] >= 0]
    ks.index = ks['id1'] + ',' + ks['id2']
    return ks[col]


def get_median(data):
    """Calculate the median of the data list."""
    if not data:
        return 0
    data_sorted = sorted(data)
    half = len(data_sorted) // 2
    return (data_sorted[half] + data_sorted[-(half + 1)]) / 2


def cds_to_pep(cds_file, pep_file, fmt='fasta'):
    """Translate CDS sequences to peptide sequences and write to file."""
    records = list(SeqIO.parse(cds_file, fmt))
    for rec in records:
        rec.seq = rec.seq.translate()
    SeqIO.write(records, pep_file, 'fasta')
    return True


def newblast(file, score, evalue, gene_loc1, gene_loc2, reverse):
    """Filter BLAST results based on score, evalue, and gene locations."""
    blast = pd.read_csv(file, sep="\t", header=None)
    
    if reverse == 'true':
        blast[[0, 1]] = blast[[1, 0]]
    blast = blast[(blast[11] >= score) & (blast[10] < evalue) & (blast[1] != blast[0])]
    blast = blast[(blast[0].isin(gene_loc1.index)) & (blast[1].isin(gene_loc2.index))]
    blast.drop_duplicates(subset=[0, 1], keep='first', inplace=True)
    blast[0] = blast[0].astype(str)
    blast[1] = blast[1].astype(str)
    return blast


def newgff(file):
    """Read GFF file and rename columns with appropriate data types."""
    gff = pd.read_csv(file, sep="\t", header=None, index_col=1)
    gff.rename(columns={0: 'chr', 2: 'start', 3: 'end', 4: 'strand', 5: 'order'}, inplace=True)
    gff['chr'] = gff['chr'].astype(str)
    gff['start'] = gff['start'].astype(np.int64)
    gff['end'] = gff['end'].astype(np.int64)
    gff['strand'] = gff['strand'].astype(str)
    gff['order'] = gff['order'].astype(int)
    return gff


def newlens(file, position):
    """Read lens file and select position based on 'order' or 'end'."""
    lens = pd.read_csv(file, sep="\t", header=None, index_col=0)
    lens.index = lens.index.astype(str)
    if position == 'order':
        lens = lens[2]
    elif position == 'end':
        lens = lens[1]
    return lens


def read_classification(file):
    """Read classification data and convert columns to appropriate types."""
    classification = pd.read_csv(file, sep="\t", header=None)
    classification[0] = classification[0].astype(str)
    classification[1] = classification[1].astype(int)
    classification[2] = classification[2].astype(int)
    classification[3] = classification[3].astype(str)
    classification[4] = classification[4].astype(int)
    return classification


def gene_location(gff, lens, step, position):
    """Calculate gene locations based on lens and step."""
    gff = gff[gff['chr'].isin(lens.index)].copy()
    if gff.empty:
        print('Stoped! \n\nChromosomes in gff file and lens file do not correspond.')
        exit(0)
    dict_chr = dict(zip(lens.index, np.append(np.array([0]), lens.cumsum()[:-1].values)))
    gff['loc'] = ''
    for name, group in gff.groupby('chr'):
        gff.loc[group.index, 'loc'] = (dict_chr[name] + group[position]) * step
    return gff


def dotplot_frame(fig, ax, lens1, lens2, step1, step2, genome1_name, genome2_name, arr, pad = 0):
    """Set up the dotplot frame with grid lines and labels."""
    for k in lens1.cumsum()[:-1] * step1:
        ax.axhline(y=k, alpha=0.8, color='black', lw=0.5)
    for k in lens2.cumsum()[:-1] * step2:
        ax.axvline(x=k, alpha=0.8, color='black', lw=0.5)
    align = dict(family='DejaVu Sans', style='italic', horizontalalignment="center", verticalalignment="center")
    yticks = lens1.cumsum() * step1 - 0.5 * lens1 * step1
    ax.set_yticks(yticks)
    ax.set_yticklabels(lens1.index, fontsize = 13, family='DejaVu Sans', style='normal')
    ax.tick_params(axis='y', which='major', pad = pad)
    ax.tick_params(axis='x', which='major', pad = pad)
    xticks = lens2.cumsum() * step2 - 0.5 * lens2 * step2
    ax.set_xticks(xticks)
    ax.set_xticklabels(lens2.index, fontsize = 13, family='DejaVu Sans', style='normal')
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    if arr[0] <= 0:
        ax.text(-0.06, 0.5, genome1_name, weight='semibold', fontsize=16, rotation=90, **align)
    else:
        ax.text(-0.06, 0.5, genome1_name, weight='semibold', fontsize=16, rotation=90, **align)
    if arr[1] < 0:
        ax.text(0.5, -0.06, genome2_name, weight='semibold', fontsize=16, **align)
    else:
        ax.text(0.5, -0.06, genome2_name, weight='semibold', fontsize=16, **align)

def Bezier3(plist, t):
    """Calculate Bezier curve of degree 3."""
    p0, p1, p2 = plist
    return p0 * (1 - t) ** 2 + 2 * p1 * t * (1 - t) + p2 * t ** 2


def Bezier4(plist, t):
    """Calculate Bezier curve of degree 4."""
    p0, p1, p2, p3, p4 = plist
    return p0 * (1 - t) ** 4 + 4 * p1 * t * (1 - t) ** 3 + 6 * p2 * t ** 2 * (1 - t) ** 2 + 4 * p3 * (1 - t) * t ** 3 + p4 * t ** 4


def Rectangle(ax, loc, height, width, color, alpha):
    """Draw a rectangle on the axes with specified properties."""
    p = mpatches.Rectangle(loc, width, height, edgecolor=None, facecolor=color, alpha=alpha)
    ax.add_patch(p)

def str_to_bool(s):
    if isinstance(s, bool):
        return s 
    return str(s).strip().lower() == 'true'