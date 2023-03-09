import re
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import wgdi.base as base


class dotplot():
    def __init__(self, options):
        self.multiple = 1
        self.score = 100
        self.evalue = 1e-5
        self.repeat_number = 20
        self.markersize = 0.5
        self.figsize = 'default'
        self.position = 'order'
        self.ancestor_top = None
        self.ancestor_left = None
        self.blast_reverse = 'False'
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        if self.ancestor_top == 'none' or self.ancestor_top == '':
            self.ancestor_top = None
        if self.ancestor_left == 'none' or self.ancestor_left == '':
            self.ancestor_left = None

    def pair_positon(self, blast, gff1, gff2, rednum, repeat_number):
        blast['color'] = ''
        blast['loc1'] = blast[0].map(gff1['loc'])
        blast['loc2'] = blast[1].map(gff2['loc'])
        bluenum = 5+rednum
        index = [group.sort_values(by=[11], ascending=[False])[:repeat_number].index.tolist()
                 for name, group in blast.groupby([0])]
        reddata = np.array([k[:rednum] for k in index], dtype=object)
        bluedata = np.array([k[rednum:bluenum] for k in index], dtype=object)
        graydata = np.array([k[bluenum:repeat_number] for k in index], dtype=object)
        if len(reddata):
            redindex = np.concatenate(reddata)
        else:
            redindex = []
        if len(bluedata):
            blueindex = np.concatenate(bluedata)
        else:
            blueindex = []
        if len(graydata):
            grayindex = np.concatenate(graydata)
        else:
            grayindex = []
        blast.loc[redindex, 'color'] = 'red'
        blast.loc[blueindex, 'color'] = 'blue'
        blast.loc[grayindex, 'color'] = 'gray'
        return blast[blast['color'].str.contains('\w')]

    def run(self):
        axis = [0, 1, 1, 0]
        left, right, top, bottom = 0.07, 0.97, 0.93, 0.03
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        step1 = 1 / float(lens1.sum())
        step2 = 1 / float(lens2.sum())
        if self.ancestor_left != None:
            axis[0] = -0.02
            lens_ancestor_left = pd.read_csv(
                self.ancestor_left, sep="\t", header=None)
            lens_ancestor_left[0] = lens_ancestor_left[0].astype(str)
            lens_ancestor_left[3] = lens_ancestor_left[3].astype(str)
            lens_ancestor_left[4] = lens_ancestor_left[4].astype(int)
            lens_ancestor_left[4] = lens_ancestor_left[4] / lens_ancestor_left[4].max()
            lens_ancestor_left = lens_ancestor_left[lens_ancestor_left[0].isin(
                lens1.index)]
        if self.ancestor_top != None:
            axis[3] = -0.02
            lens_ancestor_top = pd.read_csv(
                self.ancestor_top, sep="\t", header=None)
            lens_ancestor_top[0] = lens_ancestor_top[0].astype(str)
            lens_ancestor_top[3] = lens_ancestor_top[3].astype(str)
            lens_ancestor_top[4] = lens_ancestor_top[4].astype(int)
            lens_ancestor_top[4] = lens_ancestor_top[4] / lens_ancestor_top[4].max()
            lens_ancestor_top = lens_ancestor_top[lens_ancestor_top[0].isin(
                lens2.index)]
        if re.search('\d', self.figsize):
            self.figsize = [float(k) for k in self.figsize.split(',')]
        else:
            self.figsize = np.array(
                [1, float(lens1.sum())/float(lens2.sum())])*10
        plt.rcParams['ytick.major.pad'] = 0
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.xaxis.set_ticks_position('top')
        base.dotplot_frame(fig, ax, lens1, lens2, step1, step2,
                           self.genome1_name, self.genome2_name, [axis[0], axis[3]])
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gff1 = base.gene_location(gff1, lens1, step1, self.position)
        gff2 = base.gene_location(gff2, lens2, step2, self.position)
        if self.ancestor_top != None:
            top = top
            self.aree_left = self.ancestor_posion(ax, gff2, lens_ancestor_top, 'top')
        if self.ancestor_left != None:
            left = left
            self.aree_top = self.ancestor_posion(ax, gff1, lens_ancestor_left, 'left')

        blast = base.newblast(self.blast, int(self.score),
                              float(self.evalue), gff1, gff2, self.blast_reverse)
        if len(blast) ==0:
            print('Stoped! \n\nThe gene id in blast file does not correspond to gff1 and gff2.')
            exit(0)
        df = self.pair_positon(blast, gff1, gff2,
                               int(self.multiple), int(self.repeat_number))
        ax.scatter(df['loc2'], df['loc1'], s=float(self.markersize), c=df['color'],
                   alpha=0.5, edgecolors=None, linewidths=0, marker='o')
        ax.axis(axis)
        plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
        plt.savefig(self.savefig, dpi=300)
        plt.show()
        sys.exit(0)

    def ancestor_posion(self, ax, gff, lens, mark):
        data = []
        for index, row in lens.iterrows():
            loc1 = gff[(gff['chr'] == row[0]) & (
                gff['order'] == int(row[1]))].index
            loc2 = gff[(gff['chr'] == row[0]) & (
                gff['order'] == int(row[2]))].index
            loc1, loc2 = gff.loc[[loc1[0], loc2[0]], 'loc']
            if mark == 'top':
                width = abs(loc1-loc2)
                loc = [min(loc1, loc2), 0]
                height = -0.02
                base.Rectangle(ax, loc, height, width, row[3], row[4])
            if mark == 'left':
                height = abs(loc1-loc2)
                loc = [-0.02, min(loc1, loc2), ]
                width = 0.02
                base.Rectangle(ax, loc, height, width, row[3], row[4])
            data.append([loc, height, width, row[3], row[4]])
        return data
