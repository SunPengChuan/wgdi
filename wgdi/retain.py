import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base


class retain():
    def __init__(self, options):
        self.position = 'order'
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        if hasattr(self, 'ylim'):
            self.ylim = [float(k) for k in self.ylim.split(',')]
        else:
            self.ylim = [0, 1]
        self.colors = [str(k) for k in self.colors.split(',')]
        self.figsize = [float(k) for k in self.figsize.split(',')]

    def run(self):
        gff = base.newgff(self.gff)
        lens = base.newlens(self.lens, self.position)
        gff = gff[gff['chr'].isin(lens.index)]
        alignment = pd.read_csv(self.alignment, header=None, index_col=0)
        alignment = alignment.join(gff[['chr', self.position]], how='left')
        self.retain = self.align_chr(alignment)
        self.retain[self.retain.columns[:-2]
                    ].to_csv(self.savefile, sep='\t', header=None)
        fig, axs = plt.subplots(
            len(lens), 1, sharex=True, sharey=True, figsize=tuple(self.figsize))
        fig.add_subplot(111, frameon=False)
        align = dict(family='Arial', verticalalignment="center",
                     horizontalalignment="center")
        plt.ylabel(self.ylabel+'\n\n\n\n', fontsize=20, **align)
        for spine in plt.gca().spines.values():
            spine.set_visible(False)
        plt.tick_params(top=False, bottom=False, left=False,
                        right=False, labelleft=False, labelbottom=False)
        groups = self.retain.groupby(['chr'])
        for i in range(len(lens)):
            group = groups.get_group(lens.index[i])
            for j in self.retain.columns[:-2]:
                axs[i].plot(group['order'].values, group[j].values,
                            linestyle='-', color=self.colors[j-1], linewidth=1)
            axs[i].spines['right'].set_visible(False)
            axs[i].spines['top'].set_visible(False)
            axs[i].set_ylim(self.ylim)
            axs[i].tick_params(labelsize=12)
        align = dict(family='Arial', verticalalignment="center",
                     horizontalalignment="left")
        for i in range(len(lens)):
            x, y = axs[i].get_xlim()[1]*0.90, axs[i].get_ylim()[1]*0.5
            axs[i].text(x, y, self.refgenome +
                        str(lens.index[i]), fontsize=18, **align)
        plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
        sys.exit(0)

    def align_chr(self, alignment):
        for i in alignment.columns[:-2]:
            alignment.loc[alignment[i].str.contains('\w', na=False), i] = 1
            alignment.loc[alignment[i] == '.', i] = 0
            alignment.loc[alignment[i] == ' ', i] = 0
            alignment[i].fillna(0, inplace=True)
            for chr, group in alignment.groupby(['chr']):
                a = self.retain(group[i].values.tolist())
                alignment.loc[group.index, i] = a
        return alignment

    def retain(self, arr):
        a = []
        for i in range(0, len(arr)):
            start, end = i-int(self.step), i+int(self.step)
            if start < 0:
                start = 0
            if end > len(arr):
                end = len(arr)
            ave = sum(arr[start:end])/(end-start)
            a.append(ave)
        return a
