import re
import sys

import numpy as np
import pandas as pd
import wgdi.base as base


class block_correspondence():
    def __init__(self, options):
        self.block_len = 0
        self.tandem = True
        self.pvalue = 0.05
        self.position = 'order'
        self.block_length = 5
        self.tandem_length = 200
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.homo = [float(k) for k in self.homo.split(',')]

    def run(self):
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        bkinfo = pd.read_csv(self.blockinfo)
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo = bkinfo[(bkinfo['length'] > int(self.block_length)) & (bkinfo['chr1'].isin(
            lens1.index)) & (bkinfo['chr2'].isin(lens2.index)) & (bkinfo['pvalue'] < float(self.pvalue))]
        cor = [[k, i, 0, lens1[i], j, 0, lens2[j], float(self.homo[0]), float(self.homo[1])] for k in range(
            1, int(self.multiple)+1) for i in lens1.index for j in lens2.index]
        cor = pd.DataFrame(
            cor, columns=['sub', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'homo1', 'homo2'])
        cor['chr1'] = cor['chr1'].astype(str)
        cor['chr2'] = cor['chr2'].astype(str)
        if self.tandem == True or self.tandem == 'true' or self.tandem == 1:
            bkinfo = self.remove_tandem(bkinfo)
        arr = self.colinearity_region(cor, bkinfo)
        bkinfo.loc[bkinfo.index.isin(arr), :].to_csv(
            self.savefile, index=False)

    def remove_tandem(self, bkinfo):
        group = bkinfo[bkinfo['chr1'] == bkinfo['chr2']].copy()
        group.loc[:, 'start'] = group.loc[:, 'start1']-group.loc[:, 'start2']
        group.loc[:, 'end'] = group.loc[:, 'end1']-group.loc[:, 'end2']
        index = group[(group['start'].abs() < int(self.tandem_length)) | (
            group['end'].abs() < int(self.tandem_length))].index
        bkinfo = bkinfo.drop(index)
        return bkinfo

    def colinearity_region(self, cor, bkinfo):
        arr = []
        for chr1, group in bkinfo.groupby(['chr1']):
            group = group.sort_values(by=['start1'], ascending=[True])
            for index, row in group.iterrows():
                newcor = cor[(cor['chr1'] == row['chr1']) &
                             (cor['chr2'] == row['chr2'])]
                for sub, cor_row in newcor.groupby(['sub']):
                    if row['homo'+self.multiple] <= float(cor_row['homo1']) or row['homo'+self.multiple] >= float(cor_row['homo2']):
                        continue
                    arr.append(index)
        return arr
