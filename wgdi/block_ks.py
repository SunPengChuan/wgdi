import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base


class block_ks:
    def __init__(self, options):
        # Default parameters
        self.markersize = 0.8
        self.figsize = 'default'
        self.tandem_length = 200
        self.blockinfo_reverse = False
        self.tandem = False
        self.area = [0, 3]
        self.position = 'order'
        self.ks_col = 'ks_NG86'
        self.pvalue = 0.01
        
        # Overriding default parameters with options
        for k, v in options:
            setattr(self, str(k), v)
            print(f"{k} = {v}")
        
        # Parsing area as a float list
        self.area = [float(k) for k in str(self.area).split(',')]
        self.markersize =  float(self.markersize)
        self.tandem_length =  int(self.tandem_length)
        
        self.blockinfo_reverse =  base.str_to_bool(self.blockinfo_reverse)
        self.remove_tandem =  base.str_to_bool(self.remove_tandem)

    def block_position(self, bkinfo, lens1, lens2, step1, step2):
        pos, pairs = [], []
        
        # Create mappings for chromosome positions
        dict_y_chr = dict(zip(lens1.index, np.append([0], lens1.cumsum()[:-1].values)))
        dict_x_chr = dict(zip(lens2.index, np.append([0], lens2.cumsum()[:-1].values)))
        
        # Iterate through block information
        for _, row in bkinfo.iterrows():
            block1 = row['block1'].split('_')
            block2 = row['block2'].split('_')
            ks = row['ks'].split('_')
            
            locy_median = (dict_y_chr[row['chr1']] + 0.5 * (row['end1'] + row['start1'])) * step1
            locx_median = (dict_x_chr[row['chr2']] + 0.5 * (row['end2'] + row['start2'])) * step2
            pos.append([locx_median, locy_median, row['ks_median']])
            
            # Ensure ks length matches block length
            if len(block1) != len(ks):
                ks = ks[1:]
                
            for i in range(len(block1)):
                locy = (dict_y_chr[row['chr1']] + float(block1[i])) * step1
                locx = (dict_x_chr[row['chr2']] + float(block2[i])) * step2
                pairs.append([locx, locy, float(ks[i])])
        
        return pos, pairs

    def remove_tandem(self, bkinfo):
        # Filter for same-chromosome blocks
        group = bkinfo[bkinfo['chr1'] == bkinfo['chr2']].copy()
        
        # Calculate block start and end differences
        group['start'] = group['start1'] - group['start2']
        group['end'] = group['end1'] - group['end2']
        
        # Remove tandems based on threshold
        index = group[(group['start'].abs() <= self.tandem_length) |
                      (group['end'].abs() <= self.tandem_length)].index
        return bkinfo.drop(index)

    def run(self):
        # Initialize axis and chromosome lens
        axis = [0, 1, 1, 0]
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        
        # Parse figsize
        if re.search(r'\d', self.figsize):
            self.figsize = [float(k) for k in self.figsize.split(',')]
        else:
            self.figsize = np.array([1, float(lens1.sum()) / float(lens2.sum())]) * 10
        
        # Calculate step sizes
        step1 = 1 / float(lens1.sum())
        step2 = 1 / float(lens2.sum())
        
        # Create figure and axes
        fig, ax = plt.subplots(figsize=self.figsize)
        plt.rcParams['ytick.major.pad'] = 0
        ax.xaxis.set_ticks_position('top')
        
        # Plot dotplot frame
        base.dotplot_frame(fig, ax, lens1, lens2, step1, step2,
                           self.genome1_name, self.genome2_name, [0, 1])
        
        # Load block information
        bkinfo = pd.read_csv(self.blockinfo)
        
        # Handle reverse block information
        if self.blockinfo_reverse == True:
            bkinfo[['chr1', 'chr2']] = bkinfo[['chr2', 'chr1']]
            bkinfo[['block1', 'block2']] = bkinfo[['block2', 'block1']]
        
        # Filter block information
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo = bkinfo[(bkinfo['length'] >= int(self.block_length)) & 
                        (bkinfo['chr1'].isin(lens1.index)) & 
                        (bkinfo['chr2'].isin(lens2.index)) & 
                        (bkinfo['pvalue'] < float(self.pvalue))]
        
        # Remove tandem duplicates if required
        if self.tandem == False:
            bkinfo = self.remove_tandem(bkinfo)
        
        # Calculate positions and pairs
        pos, pairs = self.block_position(bkinfo, lens1, lens2, step1, step2)
        
        # Filter pairs by ks value
        df = pd.DataFrame(pairs, columns=['loc1', 'loc2', 'ks'])
        df = df[(df['ks'] >= self.area[0]) & (df['ks'] <= self.area[1])]
        df.drop_duplicates(inplace=True)
        
        # Plot scatter
        cm = plt.cm.get_cmap('gist_rainbow')
        sc = plt.scatter(df['loc1'], df['loc2'], s=self.markersize, c=df['ks'],
                         alpha=0.9, edgecolors=None, linewidths=0, marker='o', 
                         vmin=self.area[0], vmax=self.area[1], cmap=cm)
        
        # Add colorbar
        cbar = fig.colorbar(sc, shrink=0.5, pad=0.03, fraction=0.1)
        align = dict(family='DejaVu Sans', style='normal',
                     horizontalalignment="center", verticalalignment="center")
        cbar.set_label('Ks', labelpad=12.5, fontsize=16, **align)
        
        # Set axis and save figure
        ax.axis(axis)
        plt.subplots_adjust(left=0.09, right=0.96, top=0.93, bottom=0.03)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
