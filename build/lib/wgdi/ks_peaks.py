import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats.kde import gaussian_kde

import wgdi.base as base

class kspeaks:
    def __init__(self, options):
        # Default values
        self.tandem_length = 200
        self.figsize = 10, 6.18
        self.fontsize = 9
        self.block_length = 3
        self.area = 0, 3
        self.tandem =  True

        # Set options passed in
        for k, v in options:
            setattr(self, str(k), v)
            print(f'{str(k)} = {v}')

        # Convert string values to lists of floats
        self.homo = [float(k) for k in self.homo.split(',')]
        self.ks_area = [float(k) for k in self.ks_area.split(',')]
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.area = [float(k) for k in self.area.split(',')]
        self.pvalue = float(self.pvalue)
        self.block_length = int(self.block_length)
        self.tandem = base.str_to_bool(self.tandem)

    def remove_tandem(self, bkinfo):
        """
        Remove tandem duplications based on start and end position differences.
        """
        group = bkinfo[bkinfo['chr1'] == bkinfo['chr2']].copy()
        group.loc[:, 'start'] = group.loc[:, 'start1'] - group.loc[:, 'start2']
        group.loc[:, 'end'] = group.loc[:, 'end1'] - group.loc[:, 'end2']
        
        # Drop rows where start or end difference is within tandem length
        index = group[(group['start'].abs() <= self.tandem_length) | 
                      (group['end'].abs() <= self.tandem_length)].index
        bkinfo = bkinfo.drop(index)
        return bkinfo

    def ks_kde(self, df):
        """
        Perform kernel density estimation (KDE) on Ks data.
        """
        # Clean up 'ks' column by removing leading underscores
        df.loc[df['ks'].str.startswith('_'), 'ks'] = df.loc[df['ks'].str.startswith('_'), 'ks'].str[1:]
        
        ks = df['ks'].str.split('_')
        arr = []
        ks_ave = []
        
        # Collect individual Ks values and calculate average Ks per row
        for v in ks.values:
            v = [float(k) for k in v if float(k) >= 0]
            if len(v) == 0:
                continue
            arr.extend(v)
            ks_ave.append(sum(v) / len(v))  # Mean of each row's Ks values
        
        # KDE for three distributions: median, average, total
        kdemedian = gaussian_kde(df['ks_median'].values)
        kdemedian.set_bandwidth(bw_method=kdemedian.factor / 3.)
        
        kdeaverage = gaussian_kde(ks_ave)
        kdeaverage.set_bandwidth(bw_method=kdeaverage.factor / 3.)
        
        kdetotal = gaussian_kde(arr)
        kdetotal.set_bandwidth(bw_method=kdetotal.factor / 3.)

        return [kdemedian, kdeaverage, kdetotal]

    def run(self):
        """
        Main method to process the data, perform KDE, and generate the plot.
        """
        plt.rcParams['ytick.major.pad'] = 0
        fig, ax = plt.subplots(figsize=self.figsize)

        # Read the block info file
        bkinfo = pd.read_csv(self.blockinfo)
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo['length'] = bkinfo['length'].astype(int)

        # Filter based on block length and p-value
        bkinfo = bkinfo[(bkinfo['length'] > self.block_length) &
                        (bkinfo['pvalue'] < self.pvalue)]

        # Remove tandem duplications if needed
        if self.tandem == False:
            bkinfo = self.remove_tandem(bkinfo)

        # Further filtering based on homozygous range and Ks area
        bkinfo = bkinfo[bkinfo[f'homo{self.multiple}'] >= self.homo[0]]
        bkinfo = bkinfo[bkinfo[f'homo{self.multiple}'] <= self.homo[1]]
        bkinfo = bkinfo[bkinfo['ks_median'] >= self.ks_area[0]]
        bkinfo = bkinfo[bkinfo['ks_median'] <= self.ks_area[1]]

        # Perform KDE on the Ks data
        kdemedian, kdeaverage, kdetotal = self.ks_kde(bkinfo)

        # Define the range for the x-axis (Ks values)
        dist_space = np.linspace(self.area[0], self.area[1], 500)

        # Plot the KDE results
        ax.plot(dist_space, kdemedian(dist_space), color='red', label='block median')
        ax.plot(dist_space, kdeaverage(dist_space), color='black', label='block average')
        ax.plot(dist_space, kdetotal(dist_space), color='blue', label='all pairs')

        # Set plot labels, grid, and limits
        ax.grid()
        ax.set_xlabel(r'${K_{s}}$', fontsize=20)
        ax.set_ylabel('Frequency', fontsize=20)
        ax.tick_params(labelsize=18)
        ax.set_xlim(self.area)
        ax.legend(fontsize=20)

        # Adjust layout for better display
        plt.subplots_adjust(left=0.09, right=0.96, top=0.93, bottom=0.12)

        # Save the figure
        plt.savefig(self.savefig, dpi=500)
        plt.show()

        # Save the filtered data to CSV
        bkinfo.to_csv(self.savefile, index=False)