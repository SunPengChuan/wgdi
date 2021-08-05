
import sys

import matplotlib.pyplot as plt
import pandas as pd

import wgdi.base as base


class karyotype_mapping():
    def __init__(self, options):
        self.position = 'order'
        self.coverage = 0.1
        self.block_length = 5 
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)

    def karyotype_left(self, bkinfo, ancestor):
        data = []
        for name, group in bkinfo.groupby('chr1'):
            d1 = ancestor[ancestor[0] == name]
            for index1, row1 in group.iterrows():
                a, b = map(int, sorted([row1['start1'], row1['end1']]))
                block1 = row1['block1'].split('_')
                for index2, row2 in d1.iterrows():
                    c, d = sorted([row2[1], row2[2]])
                    h1 = len([k for k in range(a, b) if k in range(c, d)])/(b-a)
                    h2 = len([k for k in block1 if c <= int(k) <= d])
                    if h1 >= float(self.coverage) and h2 >= int(self.block_length):
                        a, b = map(int, sorted([row1['start2'], row1['end2']]))
                        data.append([row1['chr2'], a, b, row2[3], row2[4]])
        return data

    def karyotype_top(self, bkinfo, ancestor):
        data = []
        for name, group in bkinfo.groupby('chr2'):
            d1 = ancestor[ancestor[0] == name]
            for index1, row1 in group.iterrows():
                a, b = map(int, sorted([row1['start2'], row1['end2']]))
                block2 = row1['block2'].split('_')
                for index2, row2 in d1.iterrows():
                    c, d = sorted([row2[1], row2[2]])
                    h1 = len([k for k in range(a, b) if k in range(c, d)])/(b-a)
                    h2 = len([k for k in block2 if c <= int(k) <= d])
                    if h1 >= float(self.coverage) and h2 >= int(self.block_length):
                        a, b = map(int, sorted([row1['start1'], row1['end1']]))
                        data.append([row1['chr1'], a, b, row2[3], row2[4]])
        return data

    def karyotype_map(self, data, lens):
        df = pd.DataFrame(data)
        df = df.sort_values(by=[0, 1], ascending=[True, True])
        df.to_csv('out.csv')
        ancestor = []
        for chr, length in lens.items():
            group = df[df[0] == chr]
            color, classid = '', 1
            arr = []
            for index, row in group.iterrows():
                if color == row[3]:
                    arr.append(row[1])
                    arr.append(row[2])
                else:
                    if len(arr) >= 1:
                        ancestor.append(
                            [chr, min(arr), max(arr), color, classid])
                    arr = []
                    color = row[3]
                    classid = row[4]
                    arr.append(row[1])
                    arr.append(row[2])
            ancestor.append([chr, min(arr), max(arr), color, classid])
        ancestor = pd.DataFrame(ancestor)
        for chr, group in ancestor.groupby([0]):
            ancestor.loc[group.index[0], 1] = 1
            ancestor.loc[group.index[-1], 2] = lens[chr]
        return ancestor

    def run(self):
        bkinfo = pd.read_csv(self.blockinfo, index_col='id')
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo = bkinfo[bkinfo['length'] > int(self.block_length)]
        bkinfo['class1'] = ''
        bkinfo['col1'] = ''
        bkinfo['class2'] = ''
        bkinfo['col2'] = ''
        lens = base.newlens(self.the_other_lens, self.position)
        if hasattr(self, 'ancestor_top'):
            ancestor = base.read_calassfication(self.ancestor_top)
            data = self.karyotype_top(bkinfo, ancestor)
        if hasattr(self, 'ancestor_left'):
            ancestor = base.read_calassfication(self.ancestor_left)
            data = self.karyotype_left(bkinfo, ancestor)
        the_other_ancestor_file = self.karyotype_map(data, lens)
        the_other_ancestor_file .to_csv(
            self.the_other_ancestor_file, sep='\t', header=False, index=False)
