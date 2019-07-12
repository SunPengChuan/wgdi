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
        gff = base.newgff(self.gff)
        colinearity = base.read_colinearscan(self.colinearity)
        bkinfo = pd.read_csv(self.blockinfo)
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo = bkinfo[(bkinfo['length'] > int(self.block_length)) & (bkinfo['chr1'].isin(
            lens1.index)) & (bkinfo['chr2'].isin(lens2.index)) & (bkinfo['pvalue'] < float(self.pvalue))]
        if self.correspondence == 'all':
            cor = [[k, i, 0, lens1[i], j, 0, lens2[j], float(self.homo[0]), float(self.homo[1])]
                   for k in range(1, int(self.multiple)+1) for i in lens1.index for j in lens2.index]
            cor = pd.DataFrame(
                cor, columns=['sub', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'homo1', 'homo2'])
            cor.to_csv('all.coor.txt', header=True, index=False)
        else:
            cor = pd.read_csv(self.correspondence, sep=',',
                              header=None, engine='python')
            cor.columns = ['sub', 'chr1', 'start1',
                           'end1', 'chr2', 'start2', 'end2', 'homo1', 'homo2']
        cor['chr1'] = cor['chr1'].astype(str)
        cor['chr2'] = cor['chr2'].astype(str)
        if self.tandem == True or self.tandem == 'true' or self.tandem == 1:
            bkinfo = self.remove_tandem(bkinfo)
        for k in cor['sub'].drop_duplicates().values:
            gff['sub'+str(k)] = ''
        print(gff.head())
        align = self.colinearity_region(gff, colinearity, cor, bkinfo)
        # align[gff1.columns[-int(self.multiple):]
        #       ].to_csv(self.savefile, sep='\t', header=None)

    def remove_tandem(self, bkinfo):
        group = bkinfo[bkinfo['chr1'] == bkinfo['chr2']].copy()
        group.loc[:,'start'] = group.loc[:,'start1']-group.loc[:,'start2']
        group.loc[:,'end'] = group.loc[:,'end1']-group.loc[:,'end2']
        index = group[(group['start'].abs() < int(self.tandem_length)) | (
            group['end'].abs() < int(self.tandem_length))].index
        bkinfo = bkinfo.drop(index)
        return bkinfo

    def colinearity_region(self, gff, colinearity, cor, bkinfo):
        block_dict =dict(zip([int(block[0]) for block in colinearity], [block[1] for block in colinearity]))
        for chr1, group in bkinfo.groupby(['chr1']):
            group = group.sort_values(by=['length', 'start2'], ascending = [False,False])
            print(chr1,len(group))
            for index, row in group.iterrows():
                newcor = cor[(cor['chr1'] == row['chr1'])&(cor['chr2'] == row['chr2'])]
                for sub, cor_row in newcor.groupby(['sub']):
                    if row['homo'+self.multiple] <= float(cor_row['homo1']) or row['homo'+self.multiple] >= float(cor_row['homo2']):
                        continue
                    # print(index,row['length'])
                    index = gff[(gff[0] == chr1) & (gff1[5] >= start1) & (gff1[5] <= end1)].index
            break
            # if (end1-start1)/len(array1) <= 0.05 or (end2-start2)/len(array2) <= 0.05:
    #             continue
    #         if (end1-start1)/len(array1) <= 0.05 or (end2-start2)/len(array2) <= 0.05:
    #             continue
    #         for index, row in group.iterrows():
    #             if (int(row['start1']) <= start1) and (int(row['end1']) >= end1) and (int(row['start2']) <= start2) and (int(row['end2']) >= end2):
    #                 homo = 0
    #                 for block in k[1]:
    #                     if (block[0] not in gff1.index) or (block[2] not in gff2.index):
    #                         continue
    #                     if block[0]+","+block[2] in homopairs.keys():
    #                         homo += homopairs[block[0]+","+block[2]]
    #                 homo = homo/len(k[1])
    #                 if homo <= float(row['homo1']) or homo >= float(row['homo2']):
    #                     continue
    #                 index = gff1[(gff1[0] == chr1) & (gff1[5] >= start1) & (
    #                     gff1[5] <= end1)].index
    #                 new_index = [i[0] for i in k[1]]
    #                 for i in range(1, int(self.multiple)+1):
    #                     name = 'L'+str(i)
    #                     old_index = gff1[gff1.index.isin(
    #                         index) & gff1[name].str.match(r'\w+') == True].index
    #                     inters = np.intersect1d(old_index, new_index)
    #                     if len(inters)/len(new_index) > 0.2:
    #                         continue
    #                     gff1.loc[gff1.index.isin(
    #                         index) & gff1[name].isnull(), name] = '.'
    #                     gff1.loc[new_index, name] = [i[2] for i in k[1]]
    #                     break
    #     return gff1
