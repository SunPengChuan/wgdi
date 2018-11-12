import pandas as pd
import re
import sys
import wgdi.base as base


class block_correspondence():
    def __init__(self, options):
        self.block_len = 0
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def run(self):
        colinearity = base.read_colinearscan(self.colinearity)
        gff1 = pd.read_csv(self.gff1, sep='\t', header=None, index_col=1)
        gff2 = pd.read_csv(self.gff2, sep='\t', header=None, index_col=1)
        gff1[0] = gff1[0].astype(str)
        gff2[0] = gff2[0].astype(str)
        cor = pd.read_csv(self.correspondence, sep='\t|:|\-',
                          header=None, engine='python')
        cor[1] = cor[1].astype(str)
        cor[4] = cor[4].astype(str)
        cols = cor[0].drop_duplicates().values
        for k in cols:
            gff1[k] = ''
        gff1 = pd.concat([gff1, pd.DataFrame(columns=cols)], sort=True)
        align = self.colinearity_region(gff1, gff2, colinearity, cor)
        align[cols].to_csv(self.savefile, sep='\t', header=None)

    def colinearity_region(self, gff1, gff2, colinearity, cor):
        for k in colinearity:
            if len(k[0])<=self.block_len:
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
                        index = gff1[(gff1[0] == chr1) & (gff1[5] >= start1) & (
                            gff1[5] <= end1) & (gff1[name] == '')].index
                        gff1.loc[index, name] = '.'
                        new_index = [i[0] for i in k[0]]
                        gff1.loc[new_index, name] = [i[2] for i in k[0]]
        return gff1
