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
        self.colors = [str(k) for k in self.colors.split(',')]

    def pair_positon(self, alignment, loc1, loc2, colors):
        alignment.index = alignment.index.map(loc1)
        data, i = [], 0
        for k in alignment.columns:
            alignment[k] = alignment[k].map(loc2)
            alignment[k].dropna(axis=0, how='any', inplace=True)
            for index, row in alignment[k].iteritems():
                data.append([index, row, colors[i]])
            i+=1
        df = pd.DataFrame(data, columns=['loc1', 'loc2', 'color'])
        return df

    def run(self):
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
                           self.genome1_name, self.genome2_name)
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gff1 = base.gene_location(gff1, lens1, step1, self.position)
        gff2 = base.gene_location(gff2, lens2, step2, self.position)
        block_list = pd.read_csv(self.block_list, header=None)
        bkinfo = pd.read_csv(self.blockinfo,index_col=0)
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        align = self.alignment(gff1, gff2, block_list, bkinfo)
        alignment = align[gff1.columns[-int(len(block_list[0].drop_duplicates())):]]
        alignment.to_csv(self.savefile, sep='\t', header=None)
        df = self.pair_positon(
            alignment, gff1['loc'], gff2['loc'], self.colors)
        plt.scatter(df['loc2'], df['loc1'], s=float(self.markersize), c=df['color'],
                    alpha=0.5, edgecolors=None, linewidths=0, marker='o')
        plt.subplots_adjust(left=0.07, right=0.97, top=0.93, bottom=0.03)
        plt.savefig(self.savefig, dpi=500)
        sys.exit(0)

    def alignment(self, gff1, gff2, block_list, bkinfo):
        bkl = block_list[0].drop_duplicates().sort_values(ascending=True)
        for col in bkl:
            bl = block_list.loc[block_list[0]==col,1].values
            name = 'l'+str(int(col))
            gff1[name] = np.nan
            for block in [int(k) for k in bl]:
                block1 = bkinfo.loc[block, 'block1'].split(',')
                block2 = bkinfo.loc[block, 'block2'].split(',')
                block1 = [int(k) for k in block1]
                block2 = [int(k) for k in block2]
                # print(block1)
                index = gff1[(gff1['chr'] == bkinfo.loc[block, 'chr1']) & (
                    gff1['order'] >= min(block1)) & (gff1['order'] <= max(block1))].index
                index1 = gff1[(gff1['chr'] == bkinfo.loc[block, 'chr1']) & (
                    gff1['order'].isin(block1))].index
                index2 = gff2[(gff2['chr'] == bkinfo.loc[block, 'chr2']) & (
                    gff2['order'].isin(block2))].index
                if block1[0] -block1[1]>0:
                    index2 = index2[::-1]
                # old_index = gff1[gff1.index.isin(
                #     index) & gff1[name].str.match(r'\w+') == True].index
                # inters = np.intersect1d(index1, old_index)
                # if len(inters)/len(index1) > 0.2:
                #     continue
                # print(len(index1),len(index2),block1)
                gff1.loc[index1, name] = index2
                gff1.loc[gff1.index.isin(index) & gff1[name].isna(), name] = '.'
        return gff1
