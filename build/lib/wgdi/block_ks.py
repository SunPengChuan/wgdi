import re
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base


class block_ks():
    def __init__(self, options):
        self.markersize = 0.8
        self.figsize = 'default'
        self.tandem_length = 200
        self.tandem = False
        self.area = [0, 3]
        self.position = 'order'
        self.ks_col = 'ks_NG86'
        self.pvalue = 0.01
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.area = [float(k) for k in self.area.split(',')]
        if not hasattr(self, 'blockinfo_reverse'):
            self.blockinfo_reverse = 'false'

    def block_position(self, bkinfo, lens1, lens2, step1, step2):
        pos, pairs = [], []
        dict_y_chr = dict(zip(lens1.index, np.append(
            np.array([0]), lens1.cumsum()[:-1].values)))
        dict_x_chr = dict(zip(lens2.index, np.append(
            np.array([0]), lens2.cumsum()[:-1].values)))
        for index, row in bkinfo.iterrows():
            block1 = row['block1'].split('_')
            block2 = row['block2'].split('_')
            ks = row['ks'].split('_')
            locy_median = (dict_y_chr[row['chr1']] +
                           0.5*(row['end1']+row['start1']))*step1
            locx_median = (dict_x_chr[row['chr2']] +
                           0.5*(row['end2']+row['start2']))*step2
            pos.append([locx_median, locy_median, row['ks_median']])
            if len(block1) != len(ks):
                print(row['ks'])
                print(row['id'], len(block1), len(ks))
            for i in range(len(block1)):
                locy = (dict_y_chr[row['chr1']]+float(block1[i]))*step1
                locx = (dict_x_chr[row['chr2']]+float(block2[i]))*step2
                pairs.append([locx, locy, float(ks[i])])
        return pos, pairs

    def remove_tandem(self, bkinfo):
        group = bkinfo[bkinfo['chr1'] == bkinfo['chr2']].copy()
        group.loc[:, 'start'] = group.loc[:, 'start1']-group.loc[:, 'start2']
        group.loc[:, 'end'] = group.loc[:, 'end1']-group.loc[:, 'end2']
        index = group[(group['start'].abs() <= int(self.tandem_length)) | (
            group['end'].abs() <= int(self.tandem_length))].index
        bkinfo = bkinfo.drop(index)
        return bkinfo

    def run(self):
        axis = [0, 1, 1, 0]
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        if re.search('\d', self.figsize):
            self.figsize = [float(k) for k in self.figsize.split(',')]
        else:
            self.figsize = np.array(
                [1, float(lens1.sum())/float(lens2.sum())])*10
        step1 = 1 / float(lens1.sum())
        step2 = 1 / float(lens2.sum())
        fig, ax = plt.subplots(figsize=self.figsize)
        plt.rcParams['ytick.major.pad'] = 0
        ax.xaxis.set_ticks_position('top')
        base.dotplot_frame(fig, ax, lens1, lens2, step1, step2,
                           self.genome1_name, self.genome2_name, [1, 1])
        bkinfo = pd.read_csv(self.blockinfo)
        if self.blockinfo_reverse == True or self.blockinfo_reverse.upper() == 'TRUE':
            bkinfo[['chr1', 'chr2']] = bkinfo[['chr2', 'chr1']]
            bkinfo[['block1', 'block2']] = bkinfo[['block2', 'block1']]
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo = bkinfo[(bkinfo['length'] > int(self.block_length)) & (bkinfo['chr1'].isin(
            lens1.index)) & (bkinfo['chr2'].isin(lens2.index)) & (bkinfo['pvalue'] < float(self.pvalue))]
        if self.tandem == True or self.tandem == 'true' or self.tandem == 1:
            bkinfo = self.remove_tandem(bkinfo)
        pos, pairs = self.block_position(bkinfo, lens1, lens2, step1, step2)
        cm = plt.cm.get_cmap('gist_rainbow')  # gist_rainbow,hsv
        df = pd.DataFrame(pairs, columns=['loc1', 'loc2', 'ks'])
        df = df[(df['ks'] >= self.area[0]) & (df['ks'] <= self.area[1])]
        df.drop_duplicates(inplace=True)
        sc = plt.scatter(df['loc1'], df['loc2'], s=float(self.markersize), c=df['ks'],
                         alpha=0.9, edgecolors=None, linewidths=0, marker='o', vmin=self.area[0], vmax=self.area[1], cmap=cm)
        # sc = plt.scatter(df['loc1'], df['loc2'], s=float(self.markersize), c='dimgray',
        #                  alpha=0.7, edgecolors=None, linewidths=0, marker='o', vmin=self.area[0], vmax=self.area[1])
        cbar = fig.colorbar(sc, shrink=0.5, pad=0.03, fraction=0.1)
        align = dict(family='Arial', style='normal',
                     horizontalalignment="center", verticalalignment="center")
        cbar.set_label('Ks', labelpad=12.5, fontsize=18, **align)
        ax.axis(axis)
        plt.subplots_adjust(left=0.09, right=0.96, top=0.93, bottom=0.03)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
        sys.exit(0)
