import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


class ksfigure():
    def __init__(self, options):
        self.figsize = 10, 6.18
        self.legendfontsize = 9
        self.labelfontsize = 9
        self.area = 0, 3
        self.mode = 'median'
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        if self.xlabel == 'none' or self.xlabel == '':
            self.xlabel = r'Synonymous nucleotide subsititution (${K_{s}}$)'
        if self.ylabel == 'none' or self.ylabel == '':
            self.ylabel = 'No. of syntenic blocks kernel density'
        if self.title == 'none' or self.title == '':
            self.title = 'ks'
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.area = [float(k) for k in self.area.split(',')]

    def Gaussian_distribution(self, t, k):
        y = np.zeros(len(t))
        for i in range(0, int((len(k) - 1) / 3)+1):
            if np.isnan(k[3 * i + 2]):
                continue
            k[3 * i + 2] = float(k[3 * i + 2])/np.sqrt(2)
            k[3 * i + 0] = float(k[3 * i + 0])*np.sqrt(2*np.pi)*float(k[3 * i + 2])
            y1 = stats.norm.pdf(
                t, float(k[3 * i + 1]), float(k[3 * i + 2])) * float(k[3 * i + 0])
            y = y+y1
        return y

    def run(self):
        plt.rcParams['ytick.major.pad'] = 0
        fig, ax = plt.subplots(figsize=self.figsize)
        ksfit = pd.read_csv(self.ksfit, sep='\t', index_col=0)
        t = np.arange(self.area[0], self.area[1], 0.005)
        col = [k for k in ksfit.columns if re.match('Unnamed:', k)]
        for index, row in ksfit.iterrows():
            ax.plot(t, self.Gaussian_distribution(
                t, row[col].values), linestyle=row['linestyle'], color=row['color'], label=index, linewidth=row['linewidth'])

        align = dict(family='Arial', verticalalignment="center",
                     horizontalalignment="center")
        ax.set_xlabel(self.xlabel, fontsize=self.labelfontsize,
                      labelpad=20, **align)
        ax.set_ylabel(self.ylabel, fontsize=self.labelfontsize,
                      labelpad=20, **align)
        ax.set_title(self.title, weight='bold',
                     fontsize=self.labelfontsize, **align)
        plt.tick_params(labelsize=10)
        plt.legend(loc='upper right', fontsize=self.legendfontsize,
                   prop={'family': 'Arial', 'style': 'italic'}, )
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.savefig(self.savefig, dpi=300)
        plt.show()
        sys.exit(0)

