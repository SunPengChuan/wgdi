import re
import sys

import numpy as np
import pandas as pd
import wgdi.base as base


class block_correspondence():
    def __init__(self, options):
        self.block_len = 0
        self.correspondence = 'all'
        self.tandem = True
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.homo = [float(k) for k in self.homo.split(',')]

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
            cor = [[k, i, 0, lens_1.at[i, 2], j, 0, lens_2.at[j, 2]]
                   for k in range(1, int(self.wgd)+1) for i in lens_1.index for j in lens_2.index]
            cor = pd.DataFrame(
                cor, columns=['sub', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'])
            cor.to_csv('all.coor.txt', header=None, index=False)
        else:
            cor = pd.read_csv(self.correspondence, sep=',',
                              header=None, engine='python')
            cor.columns = ['sub', 'chr1', 'start1',
                           'end1', 'chr2', 'start2', 'end2']
        cor['chr1'] = cor['chr1'].astype(str)
        cor['chr2'] = cor['chr2'].astype(str)
        gff1 = pd.concat([gff1, pd.DataFrame(columns=list('L'+str(i)
                                                          for i in range(1, int(self.wgd)+1)))], sort=True)
        if self.tandem == True or self.tandem == 'true' or self.tandem == 1:
            colinearity = self.remove_tandem(colinearity, gff1, gff2)
        align = self.colinearity_region(
            gff1, gff2, colinearity, cor, homopairs)
        align[gff1.columns[-int(self.wgd):]
              ].to_csv(self.savefile, sep='\t', header=None)

    def deal_blast(self, gff1, gff2):
        blast = pd.read_csv(self.blast, sep="\t", header=None)
        score, evalue, repnum = 200, 1e-5, 20
        blast = blast[(blast[11] >= score) & (
            blast[10] < evalue) & (blast[1] != blast[0])]
        blast = blast[(blast[0].isin(gff1.index.values)) &
                      (blast[1].isin(gff2.index.values))]
        blast.drop_duplicates(subset=[0, 1], keep='first', inplace=True)
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

    def remove_tandem(self, colinearity, gff1, gff2):
        newcolinearity = []
        for k in colinearity:
            block = []
            chr1, chr2 = gff1.loc[k[0][0][0], 0], gff2.loc[k[0][0][2], 0]
            if chr1 != chr2:
                newcolinearity.append(k)
                continue
            for v in k[0]:
                if base.tendem(chr1, chr2, v[1], v[3]):
                    continue
                block.append(v)
            if len(block) == 0:
                continue
            k[0] = block
            newcolinearity.append(k)
        return newcolinearity

    def colinearity_region(self, gff1, gff2, colinearity, cor, homopairs):
        for k in colinearity:
            if len(k[0]) <= int(self.block_len):
                continue
            chr1, chr2 = gff1.loc[k[0][0][0], 0], gff2.loc[k[0][0][2], 0]
            array1, array2 = [float(i[1]) for i in k[0]], [
                float(i[3]) for i in k[0]]
            start1, end1 = min(array1), max(array1)
            start2, end2 = min(array2), max(array2)
            newcor = cor[(cor['chr1'] == str(chr1)) &
                         (cor['chr2'] == str(chr2))]
            group = newcor.drop_duplicates(
                subset=['start1', 'end1', 'start2', 'end2'], keep='first', inplace=False)
            for index, row in group.iterrows():
                if (int(row['start1']) <= start1) and (int(row['end1']) >= end1) and (int(row['start2']) <= start2) and (int(row['end2']) >= end2):
                    homo = 0
                    for block in k[0]:
                        if (block[0] not in gff1.index) or (block[2] not in gff2.index):
                            continue
                        if block[0]+","+block[2] in homopairs.keys():
                            homo += homopairs[block[0]+","+block[2]]
                    homo = homo/len(k[0])
                    if homo <= float(self.homo[0]) or homo >= float(self.homo[1]):
                        continue
                    index = gff1[(gff1[0] == chr1) & (gff1[5] >= start1) & (
                        gff1[5] <= end1)].index
                    new_index = [i[0] for i in k[0]]
                    for i in range(1, int(self.wgd)+1):
                        name = 'L'+str(i)
                        old_index = gff1[gff1.index.isin(
                            index) & gff1[name].str.match(r'\w+') == True].index
                        inters = np.intersect1d(old_index, new_index)
                        if len(inters)/len(new_index) > 0.2:
                            continue
                        gff1.loc[gff1.index.isin(
                            index) & gff1[name].isnull(), name] = '.'
                        gff1.loc[new_index, name] = [i[2] for i in k[0]]
                        break
        return gff1
