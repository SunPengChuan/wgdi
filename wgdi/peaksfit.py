import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde, linregress

import wgdi.base as base


class peaksfit():
    def __init__(self, options):
        self.figsize = 10, 6.18
        self.fontsize = 9
        self.area = 0, 3
        self.mode = 'median'
        self.histogram_only = 'false'
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.area = [float(k) for k in self.area.split(',')]
        self.bins_number = int(self.bins_number)
        self.peaks = 1

    def ks_values(self, df):
        df.loc[df['ks'].str.startswith('_'),'ks']= df.loc[df['ks'].str.startswith('_'),'ks'].str[1:]
        ks = df['ks'].str.split('_')
        ks_total = []
        ks_average = []
        for v in ks.values:
            ks_total.extend([float(k) for k in v])
        ks_average = df['ks_average'].values
        ks_median = df['ks_median'].values
        return [ks_median, ks_average, ks_total]

    def gaussian_fuc(self, x, *params):
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            amp = float(params[i])
            ctr = float(params[i+1])
            wid = float(params[i+2])
            y = y + amp * np.exp(-((x - ctr)/wid)**2)
        return y

    def kde_fit(self, data, x):
        kde = gaussian_kde(data)
        kde.set_bandwidth(bw_method=kde.factor/3.)
        p = kde(x)
        guess = [1,1, 1]*self.peaks
        popt, pcov = curve_fit(self.gaussian_fuc, x, p, guess, maxfev = 80000)
        popt = [abs(k) for k in popt]
        data = []
        y = self.gaussian_fuc(x, *popt)
        for i in range(0, len(popt), 3):
            array = [popt[i], popt[i+1], popt[i+2]]
            data.append(self.gaussian_fuc(x, *array))
        slope, intercept, r_value, p_value, std_err = linregress(p, y)
        print("\nR-square: "+str(r_value**2))
        print("The gaussian fitting curve parameters are :")
        print('  |  '.join([str(k) for k in popt]))
        return y, data

    def run(self):
        plt.rcParams['ytick.major.pad'] = 0
        fig, ax = plt.subplots(figsize=self.figsize)
        bkinfo = pd.read_csv(self.blockinfo)
        ks_median, ks_average, ks_total = self.ks_values(bkinfo)
        data = eval('ks_'+self.mode)
        data = [k for k in data if self.area[0] <= k <= self.area[1]]
        x = np.linspace(self.area[0], self.area[1], self.bins_number)
        n, bins, patches = ax.hist(data, int(
            self.bins_number), density=1, facecolor='blue', alpha=0.3, label='histogram')
        if self.histogram_only == True or self.histogram_only.upper() == 'TRUE':
            pass
        else:
            y, fit = self.kde_fit(data, x)
            ax.plot(x, y, color='black', linestyle='-', label='Gaussian fit')
        ax.grid()
        align = dict(family='Arial', verticalalignment="center",
                     horizontalalignment="center")
        ax.set_xlabel(r'${K_{s}}$', fontsize=20)
        ax.set_ylabel('Frequency', fontsize=20)
        ax.tick_params(labelsize=18)
        ax.legend(fontsize=20)
        ax.set_xlim(self.area)
        plt.subplots_adjust(left=0.09, right=0.96, top=0.93, bottom=0.12)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
        sys.exit(0)
