import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats.kde import gaussian_kde

import wgdi.base as base


class kspeaks():
    def __init__(self, options):
        self.tandem_length = 200
        self.figsize = 10, 6.18
        self.fontsize = 9
        self.block_length = 3
        self.area = 0, 3
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.homo = [float(k) for k in self.homo.split(',')]
        self.ks_area = [float(k) for k in self.ks_area.split(',')]
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.area = [float(k) for k in self.area.split(',')]

    def remove_tandem(self, bkinfo):
        group = bkinfo[bkinfo['chr1'] == bkinfo['chr2']].copy()
        group.loc[:, 'start'] = group.loc[:, 'start1']-group.loc[:, 'start2']
        group.loc[:, 'end'] = group.loc[:, 'end1']-group.loc[:, 'end2']
        index = group[(group['start'].abs() <= int(self.tandem_length)) | (
            group['end'].abs() <= int(self.tandem_length))].index
        bkinfo = bkinfo.drop(index)
        return bkinfo

    def ks_kde(self, df):
        ks = df['ks'].str.split('_')
        arr = []
        ks_ave = []
        for v in ks.values:
            v = [float(k) for k in v if float(k) >= 0]
            if len(v) == 0:
                continue
            arr.extend(v)
            ks_ave.append(sum([float(k) for k in v])/len(v))
        kdemedian = gaussian_kde(df['ks_median'].values)
        kdemedian.set_bandwidth(bw_method=kdemedian.factor / 3.)
        kdeaverage = gaussian_kde(ks_ave)
        kdeaverage.set_bandwidth(bw_method=kdeaverage.factor / 3.)
        kdetotal = gaussian_kde(arr)
        kdetotal.set_bandwidth(bw_method=kdetotal.factor / 3.)
        return [kdemedian, kdeaverage, kdetotal]

    def run(self):
        plt.rcParams['ytick.major.pad'] = 0
        fig, ax = plt.subplots(figsize=self.figsize)
        bkinfo = pd.read_csv(self.blockinfo)
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo = bkinfo[(bkinfo['length'] > int(self.block_length))
                        & (bkinfo['pvalue'] < float(self.pvalue))]
        if self.tandem == True or self.tandem == 'true' or self.tandem == 1:
            bkinfo = self.remove_tandem(bkinfo)
        bkinfo = bkinfo[bkinfo['length'] >= int(self.block_length)]
        bkinfo = bkinfo[bkinfo['homo' +
                               str(self.multiple)] >= float(self.homo[0])]
        bkinfo = bkinfo[bkinfo['homo' +
                               str(self.multiple)] <= float(self.homo[1])]
        bkinfo = bkinfo[bkinfo['ks_median'] >= float(self.ks_area[0])]
        bkinfo = bkinfo[bkinfo['ks_median'] <= float(self.ks_area[1])]
        kdemedian, kdeaverage, kdetotal = self.ks_kde(bkinfo)
        dist_space = np.linspace(self.area[0], self.area[1], 500)
        ax.plot(dist_space, kdemedian(dist_space), color='red',
                label='block median')
        ax.plot(dist_space, kdeaverage(dist_space),
                color='black', label='block average')
        ax.plot(dist_space, kdetotal(dist_space),
                color='blue', label='all pairs')
        ax.grid()
        align = dict(family='Arial', verticalalignment="center",
                     horizontalalignment="center")
        ax.set_xlabel(r'${K_{s}}$', fontsize=20)
        ax.set_ylabel('Frequency', fontsize=20)
        ax.tick_params(labelsize=18)
        ax.legend(fontsize=20)
        plt.subplots_adjust(left=0.09, right=0.96, top=0.93, bottom=0.12)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
        bkinfo.to_csv(self.savefile, index=False)
        sys.exit(0)
