import re

import numpy as np
import pandas as pd

import wgdi.base as base


class block_correspondence():
    def __init__(self, options):
        self.tandem = True
        self.pvalue = 0.2
        self.position = 'order'
        self.block_length = 5
        self.tandem_length = 200
        self.tandem_ratio = 1
        self.ks_hit = 0.5
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        if hasattr(self, 'ks_area'):
            self.ks_area = [float(k) for k in self.ks_area.split(',')]
        else:
            self.ks_area = [-1, 3]
        self.homo = [float(k) for k in self.homo.split(',')]
        self.tandem_ratio  = float(self.tandem_ratio)

    def run(self):
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        bkinfo = pd.read_csv(self.blockinfo)
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo = bkinfo[(bkinfo['length'] >= int(self.block_length)) & (bkinfo['chr1'].isin(
            lens1.index)) & (bkinfo['chr2'].isin(lens2.index)) & (bkinfo['pvalue'] <= float(self.pvalue))]
        if 'tandem_ratio' in bkinfo.columns:
            bkinfo = bkinfo[bkinfo['tandem_ratio']<=self.tandem_ratio]
        cor = [[k, i, 0, lens1[i], j, 0, lens2[j], float(self.homo[0]), float(self.homo[1])] for k in range(
            1, int(self.multiple)+1) for i in lens1.index for j in lens2.index]
        cor = pd.DataFrame(
            cor, columns=['sub', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'homo1', 'homo2'])
        cor['chr1'] = cor['chr1'].astype(str)
        cor['chr2'] = cor['chr2'].astype(str)
        if self.tandem == False or self.tandem.upper() == 'FALSE':
            bkinfo = self.remove_tandem(bkinfo)
        self.remove_ks_hit(bkinfo)
        arr = self.collinearity_region(cor, bkinfo, lens1)
        bkinfo.loc[bkinfo.index.isin(arr), :].to_csv(
            self.savefile, index=False)

    def remove_tandem(self, bkinfo):
        group = bkinfo[bkinfo['chr1'] == bkinfo['chr2']].copy()
        group.loc[:, 'start'] = group.loc[:, 'start1']-group.loc[:, 'start2']
        group.loc[:, 'end'] = group.loc[:, 'end1']-group.loc[:, 'end2']
        index = group[(group['start'].abs() <= int(self.tandem_length)) | (
            group['end'].abs() <= int(self.tandem_length))].index
        bkinfo = bkinfo.drop(index)
        return bkinfo

    def collinearity_region(self, cor, bkinfo, lens):
        arr = []
        for (chr1, chr2), group in bkinfo.groupby(['chr1', 'chr2']):
            group = group.sort_values(by=['length'], ascending=[False])
            df = pd.Series(0, index=range(1, int(lens[str(chr1)])+1))
            for index, row in group.iterrows():
                if row['homo'+self.multiple] < float(self.homo[0]) or row['homo'+self.multiple] > float(self.homo[1]):
                    continue
                b1 = row['block1'].split('_')
                df1 = df.copy()
                df1[[int(k) for k in b1]] += 1
                ratio = (len(df1[df1 > 0])-len(df[df > 0]))/len(b1)
                if ratio < 0.5:
                    continue
                df[[int(k) for k in b1]] += 1
                arr.append(index)
        return arr

    def remove_ks_hit(self, bkinfo):
        for index, row in bkinfo.iterrows():
            ks = row['ks'].split('_')
            if ks[0] == '':
                ks = list(map(float, ks[1:]))
            else:
                ks = list(map(float, ks))
            newks = [k for k in ks if self.ks_area[0] <= k <= self.ks_area[1]]
            percet = len(newks)/len(ks)
            if percet < float(self.ks_hit):
                bkinfo.drop(index, inplace=True)
        return bkinfo
