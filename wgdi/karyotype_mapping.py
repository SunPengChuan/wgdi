
import sys

import numpy as np
import pandas as pd

import wgdi.base as base


class karyotype_mapping():
    def __init__(self, options):
        self.blast_reverse = 'False'
        self.blockinfo_reverse = 'False'
        self.position = 'order'
        self.block_length = 5
        self.limit_length = 5
        self.repeat_number = 10
        self.score = 100
        self.evalue = 1e-5
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)

    def karyotype_left(self, pairs, ancestor, gff1, gff2):
        for index, row in ancestor.iterrows():
            (loc_min, loc_max) = sorted([row[1], row[2]])
            index1 = gff1[(gff1['chr'] == row[0]) & (
                gff1['order'] >= loc_min) & (gff1['order'] <= loc_max)].index
            gff1.loc[index1, 'color'] = row[3]
            gff1.loc[index1, 'classification'] = row[4]
        data = pd.merge(pairs, gff1, left_on=0,
                        right_on=gff1.index, how='left')
        data.drop_duplicates(subset=[1], inplace=True)
        data.index = data[1].values
        gff2.loc[data.index, 'color'] = data['color']
        gff2.loc[data.index, 'classification'] = data['classification']
        return gff2

    def karyotype_top(self, pairs, ancestor, gff1, gff2):
        for index, row in ancestor.iterrows():
            (loc_min, loc_max) = sorted([row[1], row[2]])
            index1 = gff2[(gff2['chr'] == row[0]) & (
                gff2['order'] >= loc_min) & (gff2['order'] <= loc_max)].index
            gff2.loc[index1, 'color'] = row[3]
            gff2.loc[index1, 'classification'] = row[4]
        gff2.to_csv('out1.csv', sep='\t')
        data = pd.merge(pairs, gff2, left_on=1,
                        right_on=gff2.index, how='left')
        data.drop_duplicates(subset=[0], inplace=True)
        data.index = data[0].values
        gff1.loc[data.index, 'color'] = data['color']
        gff1.loc[data.index, 'classification'] = data['classification']
        return gff1

    def karyotype_map(self, gff, lens):
        gff = gff[gff['chr'].isin(lens.index)]
        gff = gff[gff['color'].notnull()]
        ancestor = []
        for chr, group in gff.groupby(['chr']):
            color, classid, arr = '', 1, []
            for index, row in group.iterrows():
                if color == row['color'] and classid == row['classification']:
                    arr.append(row['order'])
                else:
                    if len(arr) >= int(self.limit_length):
                        ancestor.append(
                            [chr, min(arr), max(arr), color, classid, len(arr)])
                    arr = []
                    color = row['color']
                    classid = row['classification']
                    if len(ancestor) >= 1 and color == ancestor[-1][3] and classid == ancestor[-1][4] and chr == ancestor[-1][0]:
                        arr.append(ancestor[-1][1])
                        arr += np.random.randint(
                            ancestor[-1][1], ancestor[-1][2], size=ancestor[-1][5]-1).tolist()
                        ancestor.pop()
                    arr.append(row['order'])
            if len(arr) >= int(self.limit_length):
                ancestor.append(
                    [chr, min(arr), max(arr), color, classid, len(arr)])
        ancestor = pd.DataFrame(ancestor)
        for chr, group in ancestor.groupby([0]):
            ancestor.loc[group.index[0], 1] = 1
            ancestor.loc[group.index[-1], 2] = lens[chr]
        ancestor[4] = ancestor[4].astype(int)
        return ancestor[[0, 1, 2, 3, 4]]

    def colinear_gene_pairs(self, bkinfo, gff1, gff2):
        data = []
        bkinfo = bkinfo.sort_values(by=['length'], ascending=[True])
        for index, row in bkinfo.iterrows():
            b1 = list(map(int, row['block1'].split('_')))
            b2 = list(map(int, row['block2'].split('_')))
            newgff1 = gff1[(gff1['chr'] == row['chr1'])
                           & (gff1['order'].isin(b1))]
            newgff2 = gff2[(gff2['chr'] == row['chr2'])
                           & (gff2['order'].isin(b2))]
            for i in range(len(b1)):
                a, b = newgff1.loc[newgff1['order'] == b1[i]
                                   ].index[0], newgff2.loc[newgff2['order'] == b2[i]].index[0]
                data.append([a, b])
        data = pd.DataFrame(data)
        return data

    def new_ancestor(self, ancestor, gff1, gff2, blast):
        for index, row in ancestor.iterrows():
            if row[1] == 1:
                chr_old = row[0]
                color_old = row[3]
                class_old = row[4]
                area_mo = row[2]
                continue
            if (chr_old == row[0]) and (class_old == row[4]):
                area = row[1]-area_mo
                if area == 1:
                    area_mo = row[2]
                else:
                    if hasattr(self, 'ancestor_top'):
                        index1 = gff1[(gff1['chr'] == row[0]) & (
                            gff1['order'] >= area_mo) & (gff1['order'] <= row[1])].index
                        index2 = gff2[(gff2['color'] == color_old)].index
                        index3 = gff2[(gff2['color'] == row[3])].index
                    if hasattr(self, 'ancestor_left'):
                        index1 = gff1[(gff1['chr'] == row[0]) & (
                            gff1['order'] >= area_mo) & (gff1['order'] <= row[1])].index
                        index2 = gff2[(gff2['color'] == color_old)].index
                        index3 = gff2[(gff2['color'] == row[3])].index
                    newblast1 = blast[(blast[0].isin(index1))
                                      & (blast[1].isin(index2))]
                    newblast2 = blast[(blast[0].isin(index1))
                                      & (blast[1].isin(index3))]
                    if len(newblast1) >= len(newblast2):
                        ancestor.loc[index-1, 2] = row[1]-1
                    else:
                        ancestor.loc[index, 1] = area_mo+1
            else:
                chr_old = row[0]
                color_old = row[3]
                class_old = row[4]
            area_mo = row[2]
        return ancestor

    def run(self):
        bkinfo = pd.read_csv(self.blockinfo, index_col='id')
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        if self.blockinfo_reverse == True or self.blockinfo_reverse.upper() == 'TRUE':
            bkinfo[['chr1', 'chr2']] = bkinfo[['chr2', 'chr1']]
            bkinfo[['block1', 'block2']] = bkinfo[['block2', 'block1']]
        bkinfo = bkinfo[bkinfo['length'] > int(self.block_length)]
        bkinfo['class1'] = ''
        bkinfo['col1'] = ''
        bkinfo['class2'] = ''
        bkinfo['col2'] = ''
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        lens = base.newlens(self.the_other_lens, self.position)
        blast = base.newblast(self.blast, int(self.score), float(
            self.evalue), gff1, gff2, self.blast_reverse)
        index = [group.sort_values(by=11, ascending=False)[:int(
            self.repeat_number)].index.tolist() for name, group in blast.groupby([0])]
        blast = blast.loc[np.concatenate(
            np.array([k[:int(self.repeat_number)] for k in index], dtype=object)), [0, 1]]
        pairs = self.colinear_gene_pairs(bkinfo, gff1, gff2)
        data = []
        for k, v in lens.items():
            for i in range(1, v+1):
                data.append([k, i, 0])
        df = pd.DataFrame(data)
        gff1['color'] = np.nan
        gff2['color'] = np.nan
        gff1['classification'] = np.nan
        gff2['classification'] = np.nan
        if hasattr(self, 'ancestor_top'):
            ancestor = base.read_calassfication(self.ancestor_top)
            data = self.karyotype_top(pairs, ancestor, gff1, gff2)
        if hasattr(self, 'ancestor_left'):
            ancestor = base.read_calassfication(self.ancestor_left)
            data = self.karyotype_left(pairs, ancestor, gff1, gff2)
        the_other_ancestor_file = self.karyotype_map(data, lens)
        the_other_ancestor_file = self.new_ancestor(
            the_other_ancestor_file, gff1, gff2, blast)
        the_other_ancestor_file .to_csv(
            self.the_other_ancestor_file, sep='\t', header=False, index=False)
