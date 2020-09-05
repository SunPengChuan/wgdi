import os
import sys

import numpy as np
import pandas as pd
import wgdi.base as base


class pindex():
    def __init__(self, options):
        self.remove_delta = True
        self.position = 'order'
        self.retention = 0.05
        self.diff = 0.05
        self.gap = 50
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.gap = int(self.gap)
        self.retention = float(self.retention)
        self.diff = float(self.diff)

    def Pindex(self, sub1, sub2):
        r1 = self.retain(sub1)
        r2 = self.retain(sub2)
        r = []
        for i in range(len(r2)):
            if(r1[i] < self.retention or r2[i] < self.retention):
                r.append(0)
                continue
            d = (r1[i]-r2[i])/(r1[i]+r2[i])*0.5
            if d > self.diff:
                r.append(1)
            elif -d > self.diff:
                r.append(-1)
            else:
                r.append(0)
        a, b, c = len([i for i in r if i == 1]), len(
            [i for i in r if i == -1]), len([i for i in r if i == 0])
        return [a, -b, c, len(r)]

    def retain(self, arr):
        a = []
        for i in range(0, len(arr), 2*self.gap):
            start, end = i-self.gap, i+self.gap
            genenum, retainnum = 0, 0
            for j in range(start, end):
                if((j >= int(len(arr))) or (j < 0)):
                    continue
                else:
                    retainnum += arr[j]
                    genenum += 1
            a.append(float(retainnum/genenum))
        return a

    def run(self):
        alignment = pd.read_csv(self.alignment, header=None, index_col=0)
        alignment.replace('\w+', 1, regex=True, inplace=True)
        alignment.replace('.', 0, inplace=True)
        alignment.fillna(0, inplace=True)
        gff = base.newgff(self.gff)
        lens = base.newlens(self.lens, self.position)
        gff = gff[gff['chr'].isin(lens.index)]
        alignment = alignment.join(gff[['chr', self.position]], how='left')
        alignment.dropna(axis=0, how='any', inplace=True)
        p = self.cal_pindex(alignment)
        print('Polyploidy-index: ', p)
        sys.exit(0)

    def cal_pindex(self, alignment):
        data, df = [], []
        columns = alignment.columns[:-2].tolist()
        for i in range(len(columns)-1):
            for j in range(i+1, len(columns)):
                b = []
                for chr, group in alignment.groupby(['chr']):
                    sub1 = group.loc[:, columns[i]].tolist()
                    sub2 = group.loc[:, columns[j]].tolist()
                    p = self.Pindex(sub1, sub2)
                    b.append(p)
                    df.append([i, j, chr]+p)
                sub_diver = sum([abs(k[0]+k[1]) for k in b])
                if self.remove_delta in ['true', 'TRUE', '1']:
                    sub_total = sum([abs(k[1])+abs(k[0]) for k in b])
                    if sub_total == 0:
                        c = 0
                    else:
                        c = sub_diver/sub_total
                else:
                    sub_total = sum([abs(k[1])+abs(k[0])+abs(k[2]) for k in b])
                    c = sub_diver/sub_total
                data.append(c)
        df = pd.DataFrame(df, columns=[
                          'sub1', 'sub2', 'chr', 'sub1_high', 'sub2_high', 'No_diff', 'Total'])
        df['sub2_high'] = df['sub2_high'].abs()
        self.infomation(df)
        print('\nPolyploidy-index between subgenomes are ', data)
        return sum(data)/len(data)

    def turn_percentage(self, x):
        return '(%.2f%%)' % (x * 100)

    def infomation(self, df):
        data = []
        for names, group in df.groupby(['sub1', 'sub2']):
            newgroup = pd.concat([group.head(1), group],
                                 axis=0, ignore_index=True)
            cols = ['sub1_high', 'sub2_high', 'No_diff', 'Total']
            newgroup.loc[0, cols] = group.loc[:, cols].sum()
            group1 = newgroup.copy()
            group1[cols] = group1[cols].astype(str)
            newgroup['sub1_high'] = (
                newgroup['sub1_high'] / newgroup['Total']).apply(self.turn_percentage)
            newgroup['sub2_high'] = (
                newgroup['sub2_high'] / newgroup['Total']).apply(self.turn_percentage)
            newgroup['No_diff'] = (
                newgroup['No_diff'] / newgroup['Total']).apply(self.turn_percentage)
            newgroup['Total'] = (
                newgroup['Total'] / group['Total'].sum()).apply(self.turn_percentage)
            newgroup[cols] = group1[cols]+newgroup[cols]
            group_list = []
            a = newgroup[['chr']+cols].columns.to_numpy()
            a[0] = 'Chromosome'
            a[1], a[2] = 'Sub_'+str(names[0]+1), 'Sub_'+str(names[1]+1)
            group_list.append(a)
            b = newgroup[['chr']+cols].to_numpy()
            b[0][0] = 'Total'
            for k in b:
                group_list.append(k)
            group_list = np.array(group_list).T
            for k in group_list:
                data.append(k)
        data = pd.DataFrame(data)
        data.to_csv(self.savefile, header=None, index=None)
