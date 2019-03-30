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
        base.dotplot_frame(fig, ax, lens1, lens2, step1, step2,
                           self.genome1_name, self.genome2_name)
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gff1 = base.gene_location(gff1, lens1, step1, self.position)
        gff2 = base.gene_location(gff2, lens2, step2, self.position)
        blast = base.newblast(self.blast, int(self.score),
                              float(self.evalue), gff1, gff2)
        df = self.pair_positon(blast, gff1, gff2,
                               int(self.wgd), int(self.repnum))
        plt.scatter(df['loc2'], df['loc1'], s=float(self.markersize), c=df['color'],
                    alpha=0.5, edgecolors=None, linewidths=0, marker='o')
        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
        plt.savefig(self.savefile, dpi=500)
        sys.exit(0)

    def pair_positon(self, blast, gff1, gff2, rednum, repnum):
        blast['color'] = ''
        blast['loc1'] = blast[0].map(gff1['loc'])
        blast['loc2'] = blast[1].map(gff2['loc'])
        bluenum = 4+rednum
        index = [group[:repnum].index.tolist()
                 for name, group in blast.groupby([0])]
        redindex = np.concatenate(np.array([k[:rednum] for k in index]))
        blueindex = np.concatenate(
            np.array([k[rednum:bluenum] for k in index]))
        grayindex = np.concatenate(
            np.array([k[bluenum:repnum] for k in index]))
        blast.loc[redindex, 'color'] = 'red'
        blast.loc[blueindex, 'color'] = 'blue'
        blast.loc[grayindex, 'color'] = 'gray'
        return blast[blast['color'].str.contains('\w')]
