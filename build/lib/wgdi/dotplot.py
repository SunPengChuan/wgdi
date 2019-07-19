import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base


class dotplot():
    def __init__(self, options):
        self.multiple  = 1
        self.score = 200
        self.evalue = 1e-5
        self.repnum = 20
        self.markersize = 0.5
        self.figsize = 'default'
        self.position = 'order'
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
    
    def pair_positon(self, blast, gff1, gff2, rednum, repnum):
        blast['color'] = ''
        blast['loc1'] = blast[0].map(gff1['loc'])
        blast['loc2'] = blast[1].map(gff2['loc'])
        bluenum = 5+rednum
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
        blast = base.newblast(self.blast, int(self.score),
                              float(self.evalue), gff1, gff2)
        df = self.pair_positon(blast, gff1, gff2,
                               int(self.multiple ), int(self.repnum))
        plt.scatter(df['loc2'], df['loc1'], s=float(self.markersize), c=df['color'],
                    alpha=0.5, edgecolors=None, linewidths=0, marker='o')
        plt.subplots_adjust(left=0.07, right=0.97, top=0.93, bottom=0.03)
        plt.savefig(self.savefile, dpi=500)
        sys.exit(0)
