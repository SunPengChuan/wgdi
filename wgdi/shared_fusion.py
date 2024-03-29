import pandas as pd

import wgdi.base as base


class shared_fusion():
    def __init__(self, options):
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        if hasattr(self, 'classid'):
            self.classid = [str(k) for k in self.classid.split(',')]
        else:
            self.classid = ['class1', 'class2']
        if hasattr(self, 'limit_length'):
            self.limit_length = int(self.limit_length)
        else:
            self.limit_length = 20
        self.lens1 = self.lens1.replace(' ', '').split(',')
        self.lens2 = self.lens2.replace(' ', '').split(',')

    def run(self):
        ancestor_left = base.read_calassfication(self.ancestor_left)
        ancestor_top = base.read_calassfication(self.ancestor_top)
        bkinfo = pd.read_csv(self.blockinfo)
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo['start1'] = bkinfo['start1'].astype(int)
        bkinfo['end1'] = bkinfo['end1'].astype(int)
        bkinfo['start2'] = bkinfo['start2'].astype(int)
        bkinfo['end2'] = bkinfo['end2'].astype(int)
        bkinfo = bkinfo[(bkinfo['chr1'].isin(ancestor_left[0].values)) & (
            bkinfo['chr2'].isin(ancestor_top[0].values))]
        lens1 = pd.read_csv(self.lens1[0], sep='\t', header=None)
        lens2 = pd.read_csv(self.lens2[0], sep='\t', header=None)
        lens1[0] = lens1[0].astype(str)
        lens2[0] = lens2[0].astype(str)
        blockinfoout = self.block_fusions(bkinfo, ancestor_left, ancestor_top)
        blockinfoout = blockinfoout[(blockinfoout['breakpoints1'] == 1) & (
            blockinfoout['breakpoints2'] == 1)]
        blockinfoout = blockinfoout[(blockinfoout['break_length1'] >= self.limit_length) & (
            blockinfoout['break_length2'] >= self.limit_length)]
        blockinfoout.to_csv(self.filtered_blockinfo, index=False)
        lens1 = lens1[lens1[0].isin(blockinfoout['chr1'].values)]
        lens2 = lens2[lens2[0].isin(blockinfoout['chr2'].values)]
        lens1.to_csv(self.lens1[1], sep='\t', index=False, header=False)
        lens2.to_csv(self.lens2[1], sep='\t', index=False, header=False)

    def block_fusions(self, bkinfo, ancestor_left, ancestor_top):
        bkinfo['breakpoints1'] = 0
        bkinfo['breakpoints2'] = 0
        bkinfo['break_length1'] = 0
        bkinfo['break_length2'] = 0
        for index, row in bkinfo.iterrows():
            # for specie1
            a, b = sorted([row['start1'], row['end1']])
            d1 = ancestor_left[(ancestor_left[0] == row['chr1']) & (
                ancestor_left[2] >= a) & (ancestor_left[1] <= b)]
            if len(d1) > 1:
                bkinfo.loc[index, 'breakpoints1'] = 1
                breaklength_max = 0
                for index2, row2 in d1.iterrows():
                    length_in = len(
                        [k for k in range(a, b) if k in range(row2[1], row2[2])])
                    length_out = (b-a)-length_in
                    breaklength_max =breaklength_max  if breaklength_max > min(length_in, length_out)+1 else min(length_in, length_out)+1
                bkinfo.loc[index, 'break_length1'] = breaklength_max
            #for species2
            c, d = sorted([row['start2'], row['end2']])
            d2 = ancestor_top[(ancestor_top[0] == row['chr2']) & (
                ancestor_top[2] >= c) & (ancestor_top[1] <= d)]
            if len(d2) > 1:
                bkinfo.loc[index, 'breakpoints2'] = 1
                breaklength_max = 0
                for index2, row2 in d2.iterrows():
                    length_in = len(
                        [k for k in range(c, d) if k in range(row2[1], row2[2])])
                    length_out = (d-c)-length_in
                    breaklength_max =breaklength_max  if breaklength_max > min(length_in, length_out)+1 else min(length_in, length_out)+1
                bkinfo.loc[index, 'break_length2'] = breaklength_max
        return bkinfo
