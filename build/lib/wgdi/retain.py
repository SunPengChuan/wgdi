import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class retain():
    def __init__(self, options):
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.colors = [str(k) for k in self.colors.split(',')]
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.main()

    def main(self):
        gff = pd.read_csv(self.gff, sep="\t", header=None, index_col=1)
        alignment = pd.read_csv(self.alignment, sep="\t",
                                header=None, index_col=0)
        gff = gff[[0, 5]]
        gff.columns = ['chr', 'order']
        alignment = alignment.join(gff, how='left')
        self.retain = self.align_chr(alignment)
        self.retain[self.retain.columns[:-2]
                    ].to_csv(self.savefile, sep='\t', header=None)

    def run(self):
        chrnum = self.retain['chr'].drop_duplicates().values
        fig, axs = plt.subplots(
            len(chrnum), 1, sharex=True, sharey=True, figsize=tuple(self.figsize))
        fig.add_subplot(111, frameon=False)
        align = dict(family='Times New Roman', style='italic',
                     verticalalignment="center", horizontalalignment="center")
        plt.ylabel(self.ylabel+'\n\n\n\n', fontsize=16, **align)
        for spine in plt.gca().spines.values():
            spine.set_visible(False)
        plt.tick_params(top=False, bottom=False, left=False,
                        right=False, labelleft=False, labelbottom=False)
        groups = self.retain.groupby(['chr'])
        for i in range(len(chrnum)):
            group = groups.get_group(chrnum[i])
            for j in self.retain.columns[:-2]:
                axs[i].plot(group['order'].values, group[j].values,
                            linestyle='-', color=self.colors[j-1], linewidth=1)
            axs[i].spines['right'].set_visible(False)
            axs[i].spines['top'].set_visible(False)

        for i in range(len(chrnum)):
            x, y = axs[i].get_xlim()[1]*0.95, axs[i].get_ylim()[1]*0.5
            axs[i].text(x, y, self.refgenome+' ' +
                        str(chrnum[i]), fontsize=16, **align)
        plt.savefig(self.figurefile, dpi=500)
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
