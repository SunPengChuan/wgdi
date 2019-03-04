# usr/bin/env python
# coding:utf-8

import configparser
import os
import re

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
    data, b, flag = [], [], 0
    with open(file) as f1:
        for line in f1.readlines():
            line = line.strip()
            if re.match(r"MAXIMUM GAP", line):
                continue
            if re.match(r"the", line):
                b = []
                flag = 1
                continue
            if flag == 0:
                continue
            if re.match(r"\>LOCALE", line):
                flag = 0
                if len(b) == 0:
                    continue
                p = re.split(':', line)
                data.append([b, p[1]])
                b = []
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
        return None
    data.sort()
    half = len(data) // 2
    return (data[half] + data[~half]) / 2


def cds_to_pep(cds_file, pep_file, fmt='fasta'):
    records = list(SeqIO.parse(cds_file, fmt))
    for k in records:
        k.seq = k.seq.translate()
    SeqIO.write(records, pep_file, 'fasta')
    return True


def tendem(chr1, chr2, loc1, loc2):
    if (chr1 == chr2) and (abs(float(loc1)-float(loc2)) < 200):
        return True
    return False

if __name__ == "__main__":
    config()
