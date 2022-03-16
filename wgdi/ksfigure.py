import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


class ksfigure():
    def __init__(self, options):
        self.figsize = 10, 6.18
        self.legendfontsize = 30
        self.labelfontsize = 9
        self.area = 0, 3
        self.shadow = True
        self.mode = 'median'
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        if self.xlabel == 'none' or self.xlabel == '':
            self.xlabel = r'Synonymous nucleotide subsititution (${K_{s}}$)'
        if self.ylabel == 'none' or self.ylabel == '':
            self.ylabel = 'kernel density of syntenic blocks'
        if self.title == 'none' or self.title == '':
            self.title = ''
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.area = [float(k) for k in self.area.split(',')]

    def Gaussian_distribution(self, t, k):
        y = np.zeros(len(t))
        for i in range(0, int((len(k) - 1) / 3)+1):
            if np.isnan(k[3 * i + 2]):
                continue
            k[3 * i + 2] = float(k[3 * i + 2])/np.sqrt(2)
            k[3 * i + 0] = float(k[3 * i + 0]) * \
                np.sqrt(2*np.pi)*float(k[3 * i + 2])
            y1 = stats.norm.pdf(
                t, float(k[3 * i + 1]), float(k[3 * i + 2])) * float(k[3 * i + 0])
            y = y+y1
        return y

    def run(self):
        plt.rcParams['ytick.major.pad'] = 0
        fig, ax = plt.subplots(figsize=self.figsize)
        ksfit = pd.read_csv(self.ksfit, index_col=0)
        t = np.arange(self.area[0], self.area[1], 0.005)
        col = [k for k in ksfit.columns if re.match('Unnamed:', k)]
        for index, row in ksfit.iterrows():
            ax.plot(t, self.Gaussian_distribution(
                t, row[col].values), linestyle=row['linestyle'], color=row['color'],alpha=0.8, label=index, linewidth=row['linewidth'])
            if self.shadow == True or self.shadow.upper() == 'TRUE':
                ax.fill_between(t, 0, self.Gaussian_distribution(t, row[col].values),  color=row['color'], alpha=0.15, interpolate=True, edgecolor=None, label=index,)
        align = dict(family='Arial', verticalalignment="center",
                     horizontalalignment="center")
        ax.set_xlabel(self.xlabel, fontsize=self.labelfontsize,
                      labelpad=20, **align)
        ax.set_ylabel(self.ylabel, fontsize=self.labelfontsize,
                      labelpad=20, **align)
        ax.set_title(self.title, weight='bold',
                     fontsize=self.labelfontsize, **align)
        plt.tick_params(labelsize=10)
        handles,labels = ax.get_legend_handles_labels()
        if self.shadow == True or self.shadow.upper() == 'TRUE':
            plt.legend(handles=handles[:int(len(labels)/2)],labels=labels[:int(len(labels)/2)],loc='upper right', prop={
                   'family': 'Arial', 'style': 'italic', 'size': self.legendfontsize})
        else:
            plt.legend(handles=handles[:int(len(labels))],labels=labels[:int(len(labels))],loc='upper right', prop={
                   'family': 'Arial', 'style': 'italic', 'size': self.legendfontsize})
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
        sys.exit(0)
