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
        self.area = [0,3]
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.area = [float(k) for k in self.area.split(',')]

    def block_position(self, colinearity, loc1, loc2, ks):
        pos,pairs = [],[]
        for block in colinearity:
            a, b, blk_ks = [], [], []
            if len(block[0]) < int(self.block_length):
                continue
            for k in block[0]:
                if (k[0] not in loc1) or (k[2] not in loc2):
                    continue
                if k[0]+","+k[2] in ks.index:
                    pair_ks = ks.at[str(k[0])+","+str(k[2]), 3]
                    blk_ks.append(pair_ks)
                    a.append(loc1[k[0]])
                    b.append(loc2[k[2]])
                    pairs.append([k[0],k[2],loc1[k[0]],loc2[k[2]],pair_ks])
            if len(block[0]) < int(self.block_length) or len(a) == 0:
                continue
            x, y, l, h = min(a), min(b), max(a)-min(a), max(b)-min(b)
            pos.append([x, y, l, h, base.get_median(blk_ks)])
        return pos,pairs

    def frame(self,fig,ax,lens1,lens2,step1,step2):
        for k in lens1.cumsum()[:-1]*step1:
            ax.axhline(y=k, alpha=1, color='black', lw=0.5)
        for k in lens2.cumsum()[:-1]*step2:
            ax.axvline(x=k, alpha=1, color='black', lw=0.5)
        align = dict(family='Times New Roman', style='normal',
                    horizontalalignment="center", verticalalignment="center")
        my_yticks = lens1.cumsum()*step1-0.5*lens1*step1
        plt.yticks(my_yticks.values, lens1.index, fontsize=12, **align)
        my_xticks = lens2.cumsum()*step2-0.5*lens2*step2
        plt.xticks(my_xticks.values, lens2.index, fontsize=12, **align)
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')
        ax.axis([0, 1, 1, 0])
        ax.set_ylabel("Test title",labelpad = 12.5,fontsize=18, **align)
        fig.suptitle('Test title', fontsize=18, **align)
        

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
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.xaxis.set_ticks_position('top')
        self.frame(fig,ax,lens1,lens2,step1,step2)
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gene_loc1 = base.gene_location(gff1, lens1, step1, self.position)
        gene_loc2 = base.gene_location(gff2, lens2, step2, self.position)
        colinearity = base.read_colinearscan(self.colinearity)
        ks = base.read_ks(self.ks)
        pos,pairs = self.block_position(colinearity, gene_loc1, gene_loc2, ks)
        colors = ['red', 'blue', 'grey']
        cm = plt.cm.get_cmap('RdYlBu')
        df = pd.DataFrame(pairs,columns=['id1','id2','loc1','loc2','ks'])
        df.drop_duplicates(inplace=True)
        for k in pos:
            x,y =k[0]+0.5*k[2], k[1]+0.5*k[3]
            plt.text(y, x, round(k[4], 2), color='red', fontsize=6)
        sc=plt.scatter(df['loc2'], df['loc1'], s=float(self.markersize), c=df['ks'],
                    alpha=0.8, edgecolors=None, linewidths=0, marker='o',vmin=self.area[0],vmax=self.area[1],cmap=cm)
        cbar = fig.colorbar(sc, shrink=0.5,pad=0.03,frameon=False)
        align = dict(family='Times New Roman', style='normal',
                     horizontalalignment="center", verticalalignment="center")
        cbar.set_label('Ks',labelpad = 12.5,fontsize=18,**align)
        plt.subplots_adjust(left=0.06, right=1, top=0.95, bottom=0.05)
        plt.savefig(self.savefile, dpi=500)
        plt.show()
        sys.exit(0)
