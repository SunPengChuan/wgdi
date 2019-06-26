import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base


class block_ks():
    def __init__(self, options):
        self.markersize = 0.8
        self.figsize = 'default'
        self.area = [0, 3]
        self.position = 'order'
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.area = [float(k) for k in self.area.split(',')]

    def block_position(self, colinearity, gff1, gff2, ks):
        pos, pairs = [], []
        for block in colinearity:
            a, b, blk_ks = [], [], []
            if len(block[1]) < int(self.block_length):
                continue
            if block[1][0][0] not in gff1.index or block[1][0][2] not in gff2.index:
                continue
            for k in block[1]:
                if k[0]+","+k[2] in ks.index:
                    pair_ks = ks.at[str(k[0])+","+str(k[2]), 3]
                    loc1, loc2 = gff1.loc[k[0], 'loc'], gff2.loc[k[2], 'loc']
                    a.append(loc1)
                    b.append(loc2)
                    pairs.append([k[0], k[2], loc1, loc2, pair_ks])
                    blk_ks.append(pair_ks)
            if len(a) == 0:
                continue
            x, y, l, h = min(a), min(b), max(a)-min(a), max(b)-min(b)
            pos.append([x, y, l, h, base.get_median(blk_ks)])
        return pos, pairs

    def run(self):
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        if re.search('\d', self.figsize):
            self.figsize = [float(k) for k in self.figsize.split(',')]
        else:
            self.figsize = np.array(
                [1, float(lens1.sum())/float(lens2.sum())])*10
        step1 = 1 / float(lens1.sum())
        step2 = 1 / float(lens2.sum())
        plt.rcParams['ytick.major.pad'] = 0
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.xaxis.set_ticks_position('top')
        base.dotplot_frame(fig, ax, lens1, lens2, step1, step2,
                           self.genome1_name, self.genome2_name)
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gff1 = base.gene_location(gff1, lens1, step1, self.position)
        gff2 = base.gene_location(gff2, lens2, step2, self.position)
        colinearity = base.read_colinearscan(self.colinearity)
        ks = base.read_ks(self.ks)
        pos, pairs = self.block_position(colinearity, gff1, gff2, ks)
        cm = plt.cm.get_cmap('gist_rainbow')
        df = pd.DataFrame(pairs, columns=['id1', 'id2', 'loc1', 'loc2', 'ks'])
        df.drop_duplicates(inplace=True)
        for k in pos:
            x, y = k[0]+0.5*k[2], k[1]+0.5*k[3]
            plt.text(y, x, round(k[4], 2), color='red', fontsize=6)
        sc = plt.scatter(df['loc2'], df['loc1'], s=float(self.markersize), c=df['ks'],
                         alpha=0.5, edgecolors=None, linewidths=0, marker='o', vmin=self.area[0], vmax=self.area[1], cmap=cm)
        cbar = fig.colorbar(sc, shrink=0.5, pad=0.03, fraction=0.1)
        align = dict(family='Arial', style='normal',
                     horizontalalignment="center", verticalalignment="center")
        cbar.set_label('Ks', labelpad=12.5, fontsize=18, **align)
        plt.subplots_adjust(left=0.09, right=0.96, top=0.93, bottom=0.03)
        plt.savefig(self.savefile, dpi=500)
        sys.exit(0)
