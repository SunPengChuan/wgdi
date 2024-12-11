import pandas as pd
import wgdi.base as base

class shared_fusion:
    def __init__(self, options):
        for k, v in options:
            setattr(self, str(k), v)
            print(f"{k} = {v}")
        
        # Handle classid and limit_length options
        self.classid = [str(k) for k in self.classid.split(',')] if hasattr(self, 'classid') else ['class1', 'class2']
        self.limit_length = int(self.limit_length) if hasattr(self, 'limit_length') else 20
        
        # Clean and split lens files
        self.lens1 = self.lens1.replace(' ', '').split(',')
        self.lens2 = self.lens2.replace(' ', '').split(',')

    def run(self):
        # Read classification files and block information
        ancestor_left = base.read_classification(self.ancestor_left)
        ancestor_top = base.read_classification(self.ancestor_top)
        bkinfo = pd.read_csv(self.blockinfo)

        # Preprocess blockinfo columns
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo['start1'] = bkinfo['start1'].astype(int)
        bkinfo['end1'] = bkinfo['end1'].astype(int)
        bkinfo['start2'] = bkinfo['start2'].astype(int)
        bkinfo['end2'] = bkinfo['end2'].astype(int)
        
        # Filter based on ancestor chromosomes
        bkinfo = bkinfo[(bkinfo['chr1'].isin(ancestor_left[0].values)) & 
                        (bkinfo['chr2'].isin(ancestor_top[0].values))]

        # Read lens files
        lens1 = pd.read_csv(self.lens1[0], sep='\t', header=None)
        lens2 = pd.read_csv(self.lens2[0], sep='\t', header=None)
        lens1[0] = lens1[0].astype(str)
        lens2[0] = lens2[0].astype(str)

        # Perform block fusion analysis
        blockinfoout = self.block_fusions(bkinfo, ancestor_left, ancestor_top)

        # Apply filters based on breakpoints and length
        blockinfoout = blockinfoout[(blockinfoout['breakpoints1'] == 1) & 
                                     (blockinfoout['breakpoints2'] == 1)]
        blockinfoout = blockinfoout[(blockinfoout['break_length1'] >= self.limit_length) & 
                                     (blockinfoout['break_length2'] >= self.limit_length)]

        # Save the filtered block info
        blockinfoout.to_csv(self.filtered_blockinfo, index=False)

        # Filter lens data based on the blockinfoout
        lens1 = lens1[lens1[0].isin(blockinfoout['chr1'].values)]
        lens2 = lens2[lens2[0].isin(blockinfoout['chr2'].values)]

        # Save filtered lens data
        lens1.to_csv(self.lens1[1], sep='\t', index=False, header=False)
        lens2.to_csv(self.lens2[1], sep='\t', index=False, header=False)

    def block_fusions(self, bkinfo, ancestor_left, ancestor_top):
        # Initialize new columns in the bkinfo dataframe
        bkinfo['breakpoints1'] = 0
        bkinfo['breakpoints2'] = 0
        bkinfo['break_length1'] = 0
        bkinfo['break_length2'] = 0

        for index, row in bkinfo.iterrows():
            # Process species 1 (chr1)
            a, b = sorted([row['start1'], row['end1']])
            d1 = ancestor_left[(ancestor_left[0] == row['chr1']) & 
                               (ancestor_left[2] >= a) & (ancestor_left[1] <= b)]
            if len(d1) > 1:
                bkinfo.loc[index, 'breakpoints1'] = 1
                breaklength_max = 0
                for _, row2 in d1.iterrows():
                    length_in = len([k for k in range(a, b) if k in range(row2[1], row2[2])])
                    length_out = (b - a) - length_in
                    breaklength_max = max(breaklength_max, min(length_in, length_out) + 1)
                bkinfo.loc[index, 'break_length1'] = breaklength_max

            # Process species 2 (chr2)
            c, d = sorted([row['start2'], row['end2']])
            d2 = ancestor_top[(ancestor_top[0] == row['chr2']) & 
                              (ancestor_top[2] >= c) & (ancestor_top[1] <= d)]
            if len(d2) > 1:
                bkinfo.loc[index, 'breakpoints2'] = 1
                breaklength_max = 0
                for _, row2 in d2.iterrows():
                    length_in = len([k for k in range(c, d) if k in range(row2[1], row2[2])])
                    length_out = (d - c) - length_in
                    breaklength_max = max(breaklength_max, min(length_in, length_out) + 1)
                bkinfo.loc[index, 'break_length2'] = breaklength_max

        return bkinfo
