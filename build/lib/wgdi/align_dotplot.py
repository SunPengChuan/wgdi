import re
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base


class align_dotplot():
    def __init__(self, options):
        self.position = 'order'
        self.figsize = 'default'
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        if hasattr(self, 'ks_area'):
            self.ks_area = [float(k) for k in self.ks_area.split(',')]
        else:
            self.ks_area = [-1, 3]
        if hasattr(self, 'colors'):
            self.colors = [str(k) for k in self.colors.split(',')]
        else:
            self.colors = ['red', 'blue', 'green', 'black', 'orange']
        if not hasattr(self, 'blockinfo_reverse'):
            self.blockinfo_reverse = 'false'

    def pair_positon(self, alignment, loc1, loc2, colors):
        alignment.index = alignment.index.map(loc1)
        data, i = [], 0
        for k in alignment.columns:
            df = alignment[k].map(loc2)
            df.dropna(axis=0, how='any', inplace=True)
            for index, row in df.iteritems():
                data.append([index, row, colors[i]])
            i += 1
        df = pd.DataFrame(data, columns=['loc1', 'loc2', 'color'])
        return df

    def run(self):
        axis = [0, 1, 1, 0]
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        if re.search('\d', self.figsize):
            self.figsize = [float(k) for k in self.figsize.split(',')]
        else:
            self.figsize = np.array(
                [1, float(lens1.sum())/float(lens2.sum())])*10
        plt.rcParams['ytick.major.pad'] = 0
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.xaxis.set_ticks_position('top')
        step1 = 1 / float(lens1.sum())
        step2 = 1 / float(lens2.sum())
        base.dotplot_frame(fig, ax, lens1, lens2, step1, step2,
                           self.genome1_name, self.genome2_name, [0, 1])
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gff1 = base.gene_location(gff1, lens1, step1, self.position)
        gff2 = base.gene_location(gff2, lens2, step2, self.position)
        bkinfo = pd.read_csv(self.blockinfo, index_col=0)
        if self.blockinfo_reverse == True or self.blockinfo_reverse.upper() == 'TRUE':
            bkinfo[['chr1', 'chr2']] = bkinfo[['chr2', 'chr1']]
            bkinfo[['block1', 'block2']] = bkinfo[['block2', 'block1']]
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo['class'] = bkinfo['class'].astype(str)
        align = self.alignment(gff1, gff2, bkinfo)
        alignment = align[gff1.columns[-int(
            len(bkinfo['class'].drop_duplicates())):]]
        alignment.to_csv(self.savefile, header=None)
        df = self.pair_positon(
            alignment, gff1['loc'], gff2['loc'], self.colors)
        plt.scatter(df['loc2'], df['loc1'], s=float(self.markersize), c=df['color'],
                    alpha=0.5, edgecolors=None, linewidths=0, marker='o')
        ax.axis(axis)
        plt.subplots_adjust(left=0.07, right=0.97, top=0.93, bottom=0.03)
        plt.savefig(self.savefig, dpi=500)
        sys.exit(0)

    def alignment(self, gff1, gff2, bkinfo):
        for cl, group in bkinfo.groupby(['class']):
            name = 'l'+cl
            gff1[name] = np.nan
            group = group.sort_values(by=['length'], ascending=[True])
            for index, row in group.iterrows():
                b1 = row['block1'].split('_')
                b2 = row['block2'].split('_')
                ks = row['ks'].split('_')
                ks = list(map(float, ks))
                block1, block2 = [], []
                for i in range(len(ks)):
                    if self.ks_area[0] <= ks[i] <= self.ks_area[1]:
                        block1.append(int(b1[i]))
                        block2.append(int(b2[i]))
                block1 = list(map(int, block1))
                block2 = list(map(int, block2))
                area = gff1[(gff1['chr'] == row['chr1']) & (
                    gff1['order'] >= min(block1)) & (gff1['order'] <= max(block1))].index
                index1 = gff1[(gff1['chr'] == row['chr1']) &
                              (gff1['order'].isin(block1))].index
                index2 = gff2[(gff2['chr'] == row['chr2']) &
                              (gff2['order'].isin(block2))].index
                if int(b1[0]) - int(b1[1]) > 0:
                    index2 = index2[::-1]
                gff1.loc[index1, name] = index2
                gff1.loc[gff1.index.isin(area) & gff1[name].isna(), name] = '.'
        return gff1
