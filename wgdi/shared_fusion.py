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
        bkinfo = bkinfo[(bkinfo['chr1'].isin(ancestor_left[0].values)) & (
            bkinfo['chr2'].isin(ancestor_top[0].values))]
        lens1 = pd.read_csv(self.lens1[0], sep='\t', header=None)
        lens2 = pd.read_csv(self.lens2[0], sep='\t', header=None)
        lens1[0] = lens1[0].astype(str)
        lens2[0] = lens2[0].astype(str)
        data = []
        for name, group in bkinfo.groupby('chr1'):
            d1 = ancestor_left[ancestor_left[0] == name]
            for index1, row1 in group.iterrows():
                a, b = sorted([row1['start1'], row1['end1']])
                a, b = int(a), int(b)
                for index2, row2 in d1.iterrows():
                    c, d = sorted([row2[1], row2[2]])
                    length_in = len(
                        [k for k in range(a, b) if k in range(c, d)])
                    length_out = (b-a)-length_in
                    if length_in > self.limit_length and length_out > self.limit_length:
                        data.append(
                            [row1['id'], row2[3], row2[4], length_in, length_out])

        for name, group in bkinfo.groupby('chr2'):
            d2 = ancestor_top[ancestor_top[0] == name]
            for index1, row1 in group.iterrows():
                a, b = sorted([row1['start2'], row1['end2']])
                a, b = int(a), int(b)
                for index2, row2 in d2.iterrows():
                    c, d = sorted([row2[1], row2[2]])
                    length_in = len(
                        [k for k in range(a, b) if k in range(c, d)])
                    length_out = (b-a)-length_in
                    if length_in > self.limit_length and length_out > self.limit_length:
                        data.append(
                            [row1['id'], row2[3], row2[4], length_in, length_out])

        df = pd.DataFrame(data, columns=['id', 'color', 'class', 'in', 'out'])
        df.to_csv(self.savefile, index=False)
        df.drop_duplicates(subset=['id'], keep='first', inplace=True)
        blockinfoout = bkinfo[bkinfo['id'].isin(df['id'].values)]
        blockinfoout.to_csv(self.filtered_blockinfo, index=False)
        lens1 = lens1[lens1[0].isin(blockinfoout['chr1'].values)]
        lens2 = lens2[lens2[0].isin(blockinfoout['chr2'].values)]
        lens1.to_csv(self.lens1[1], sep='\t', index=False, header=False)
        lens2.to_csv(self.lens2[1], sep='\t', index=False, header=False)
