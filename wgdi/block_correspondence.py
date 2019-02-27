import re
import sys

import pandas as pd
import wgdi.base as base


class block_correspondence():
    def __init__(self, options):
        self.block_len = 0
        self.correspondence =  'all'
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.homo =  [float(k) for k in self.homo.split(',')]
    def run(self):
        colinearity = base.read_colinearscan(self.colinearity)
        gff1 = pd.read_csv(self.gff1, sep='\t', header=None, index_col=1)
        gff2 = pd.read_csv(self.gff2, sep='\t', header=None, index_col=1)
        gff1[0] = gff1[0].astype(str)
        gff2[0] = gff2[0].astype(str)
        homopairs = self.deal_blast(gff1, gff2)
        lens_1 = pd.read_csv(self.lens1, sep="\t", header=None, index_col=0)
        lens_2 = pd.read_csv(self.lens2, sep="\t", header=None, index_col=0)
        if self.correspondence == 'all':

            cor = gff1
        else:
            
            





        cor[1] = cor[1].astype(str)
        cor[4] = cor[4].astype(str)
        cols = cor[0].drop_duplicates().values
        for k in cols:
            gff1[k] = ''
        gff1 = pd.concat([gff1, pd.DataFrame(columns=cols)], sort=True)
        align = self.colinearity_region(
            gff1, gff2, colinearity, cor, homopairs)
        align[cols].to_csv(self.savefile, sep='\t', header=None)

    def deal_blast(self, gff1, gff2):
        blast = pd.read_csv(self.blast, sep="\t", header=None)
        score, evalue, repnum = 200, 1e-5, 20
        blast = blast[(blast[11] >= score) & (
            blast[10] < evalue) & (blast[1] != blast[0])]
        blast = blast[(blast[0].isin(gff1.index.values)) &
                      (blast[1].isin(gff2.index.values))]
        homopairs = {}
        dupnum = int(self.wgd)
        hitnum = dupnum+4
        for name, group in blast.groupby([0])[1]:
            newgroup = group.values.tolist()[:repnum]
            for i, el in enumerate(newgroup, start=1):
                if i <= dupnum:
                    homopairs[name+","+el] = 1
                elif i <= hitnum:
                    homopairs[name+","+el] = 0
                else:
                    homopairs[name+","+el] = -1
        return homopairs

    def colinearity_region(self, gff1, gff2, colinearity, cor, homopairs):
        for k in colinearity:
            if len(k[0]) <= int(self.block_len):
                continue
            chr1, chr2 = gff1.loc[k[0][0][0], 0], gff2.loc[k[0][0][2], 0]
            array1, array2 = [float(i[1]) for i in k[0]], [
                float(i[3]) for i in k[0]]
            start1, end1 = min(array1), max(array1)
            start2, end2 = min(array2), max(array2)
            for name, group in cor.groupby([0]):
                group.columns = ['sub', 'chr1', 'start1',
                                 'end1', 'chr2', 'start2', 'end2']
                for index, row in group.iterrows():
                    if (str(chr1) == row['chr1']) and (int(row['start1']) <= start1) and \
                        (int(row['end1']) >= end1) and (str(chr2) == row['chr2']) and \
                            (int(row['start2']) <= start2) and (int(row['end2']) >= end2):
                        homo = 0
                        for block in k[0]:
                            if (block[0] not in gff1.index) or (block[2] not in gff2.index):
                                continue
                            if block[0]+","+block[2] in homopairs.keys():
                                homo += homopairs[block[0]+","+block[2]]
                        if homo >= float(self.homo[0]) or homo <= float(self.homo[1]):
                            continue
                        index = gff1[(gff1[0] == chr1) & (gff1[5] >= start1) & (
                            gff1[5] <= end1)].index
                        length1 = len(gff1[gff1.index.isin(
                            index) & gff1[name].str.match(r'\w+') == True])
                        if length1 > len(k[0]):
                            continue
                        gff1.loc[index, name] = '.'
                        new_index = [i[0] for i in k[0]]
                        gff1.loc[new_index, name] = [i[2] for i in k[0]]
        return gff1
