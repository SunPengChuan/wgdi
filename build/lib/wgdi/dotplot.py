import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base


class dotplot():
    def __init__(self, options):
        self.wgd = 1
        self.score = 200
        self.evalue = 1e-5
        self.repnum = 20
        self.markersize = 0.5
        self.figsize = 'default'
        self.position = 'order'
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def run(self):
        if re.search('\d', self.figsize):
            self.figsize = [float(k) for k in self.figsize.split(',')]
        else:
            self.figsize = np.array(
                [1, float(lens1.sum())/float(lens2.sum())])*10
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.xaxis.set_ticks_position('top')
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        step1 = 1 / float(lens1.sum())
        step2 = 1 / float(lens2.sum())
        base.dotplot_frame(fig, ax, lens1, lens2, step1, step2)
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gene_loc1 = base.gene_location(gff1, lens1, step1, self.position)
        gene_loc2 = base.gene_location(gff2, lens2, step2, self.position)
        blast = base.newblast(self.blast, int(self.score), float(self.evalue), gene_loc1, gene_loc2)
        df = self.pair_positon(blast, gene_loc1, gene_loc2, int(self.wgd), int(self.repnum))
        print(df)
        plt.scatter(df['loc2'], df['loc1'], s=float(self.markersize), c=df['color'],
                    alpha=0.5, edgecolors=None, linewidths=0, marker='o')
        plt.subplots_adjust(left=0.02, right=1, top=0.98, bottom=0)
        plt.savefig(self.savefile, dpi=500)
        sys.exit(0)

    def pair_positon(self, blast, loc1, loc2, rednum, repnum):
        blast = blast.head(1000)
        print(blast.head())
        blast['color'] = ''
        blast['loc1'] = blast[0].map(loc1)
        blast['loc2'] = blast[1].map(loc2)
        bluenum = 4+rednum
        for name, group in blast.groupby([0]):
            blast.loc[group[:rednum].index, 'color'] = 'red'
            blast.loc[group[rednum:bluenum].index, 'color'] = 'blue'
            blast.loc[group[bluenum:repnum].index, 'color'] = 'gray'
        return blast[blast['color'] == 'red']
