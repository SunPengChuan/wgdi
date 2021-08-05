import sys

import matplotlib.pyplot as plt
import pandas as pd

import wgdi.base as base


class karyotype():
    def __init__(self, options):
        self.width = 0.5
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        if hasattr(self, 'figsize'):
            self.figsize = [float(k) for k in self.figsize.split(',')]
        else:
            self.figsize = 10, 6.18
        if hasattr(self, 'width'):
            self.width = float(self.width)
        else:
            self.width = 0.5

    def run(self):
        fig, ax = plt.subplots(figsize=self.figsize)
        ancestor_lens = pd.read_csv(
            self.ancestor, sep="\t", header=None)
        ancestor_lens[0] = ancestor_lens[0].astype(str)
        ancestor_lens[3] = ancestor_lens[3].astype(str)
        ancestor_lens[4] = ancestor_lens[4].astype(int)
        ancestor_lens[4] = ancestor_lens[4] / ancestor_lens[4].max()
        chrs = ancestor_lens[0].drop_duplicates().to_list()
        ax.bar(chrs, 10, color='white', alpha=0)
        for index, row in ancestor_lens.iterrows():
            base.Rectangle(ax, [chrs.index(row[0])-self.width*0.5,
                                row[1]], row[2]-row[1], self.width, row[3], row[4])
        ax.tick_params(labelsize=15)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
