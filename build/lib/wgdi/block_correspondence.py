import re
import numpy as np
import pandas as pd
import wgdi.base as base

class block_correspondence():
    def __init__(self, options):
        # Default values
        self.tandem = True
        self.pvalue = 0.2
        self.position = 'order'
        self.block_length = 5
        self.tandem_length = 200
        self.tandem_ratio = 1
        self.ks_hit = 0.5

        # Set user-defined options
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

        # Parse ks_area and homo if present
        self.ks_area = [float(k) for k in getattr(self, 'ks_area', '-1,3').split(',')]
        self.homo = [float(k) for k in self.homo.split(',')]
        self.tandem_ratio = float(self.tandem_ratio)

    def run(self):
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        
        # Load block information from CSV
        bkinfo = pd.read_csv(self.blockinfo)
        bkinfo = self.preprocess_blockinfo(bkinfo, lens1, lens2)
        
        # Initialize correspondence DataFrame
        cor = self.initialize_correspondence(lens1, lens2)
        
        # If no tandem allowed, remove tandem regions
        if not self.tandem:
            bkinfo = self.remove_tandem(bkinfo)
        
        # Remove low KS hits
        bkinfo = self.remove_ks_hit(bkinfo)

        # Find collinearity regions and save results
        collinear_indices = self.collinearity_region(cor, bkinfo, lens1)
        bkinfo.loc[bkinfo.index.isin(collinear_indices), :].to_csv(self.savefile, index=False)

    def preprocess_blockinfo(self, bkinfo, lens1, lens2):
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        
        # Filter by length, chromosome indices, and p-value
        bkinfo = bkinfo[(bkinfo['length'] >= int(self.block_length)) & 
                        (bkinfo['chr1'].isin(lens1.index)) & 
                        (bkinfo['chr2'].isin(lens2.index)) & 
                        (bkinfo['pvalue'] <= float(self.pvalue))]
        
        # Filter by tandem ratio if the column exists
        if 'tandem_ratio' in bkinfo.columns:
            bkinfo = bkinfo[bkinfo['tandem_ratio'] <= self.tandem_ratio]
        
        return bkinfo

    def initialize_correspondence(self, lens1, lens2):
        # Create correspondence DataFrame with initial values
        cor = [[k, i, 0, lens1[i], j, 0, lens2[j], float(self.homo[0]), float(self.homo[1])] 
               for k in range(1, int(self.multiple) + 1) 
               for i in lens1.index 
               for j in lens2.index]
        
        cor = pd.DataFrame(cor, columns=['sub', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'homo1', 'homo2'])
        cor['chr1'] = cor['chr1'].astype(str)
        cor['chr2'] = cor['chr2'].astype(str)
        
        return cor

    def remove_tandem(self, bkinfo):
        # Remove tandem regions from the DataFrame
        group = bkinfo[bkinfo['chr1'] == bkinfo['chr2']].copy()
        group['start'] = group['start1'] - group['start2']
        group['end'] = group['end1'] - group['end2']
        tandem_condition = (group['start'].abs() <= int(self.tandem_length)) | (group['end'].abs() <= int(self.tandem_length))
        index_to_remove = group[tandem_condition].index
        return bkinfo.drop(index_to_remove)

    def remove_ks_hit(self, bkinfo):
        # Remove records with insufficient KS hits
        for index, row in bkinfo.iterrows():
            ks = self.get_ks_value(row['ks'])
            ks_ratio = len([k for k in ks if self.ks_area[0] <= k <= self.ks_area[1]]) / len(ks)
            if ks_ratio < self.ks_hit:
                bkinfo.drop(index, inplace=True)
        return bkinfo

    def get_ks_value(self, ks_str):
        # Extract and return KS values as floats
        ks = ks_str.split('_')
        ks = list(map(float, ks[1:])) if ks[0] == '' else list(map(float, ks))
        return ks

    def collinearity_region(self, cor, bkinfo, lens):
        collinear_indices = []
        for (chr1, chr2), group in bkinfo.groupby(['chr1', 'chr2']):
            group = group.sort_values(by=['length'], ascending=False)
            df = pd.Series(0, index=range(1, int(lens[str(chr1)]) + 1))
            for index, row in group.iterrows():
                # Check homology conditions
                if not self.is_valid_homo(row):
                    continue
                # Update the block series and compute ratio
                b1 = [int(k) for k in row['block1'].split('_')]
                df1 = df.copy()
                df1[b1] += 1
                ratio = (len(df1[df1 > 0]) - len(df[df > 0])) / len(b1)
                if ratio < 0.5:
                    continue
                df[b1] += 1
                collinear_indices.append(index)
        
        return collinear_indices

    def is_valid_homo(self, row):
        # Check if the homology values are within the specified range
        return self.homo[0] <= row['homo' + self.multiple] <= self.homo[1]
