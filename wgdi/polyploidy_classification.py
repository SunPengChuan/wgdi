import re
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base


class polyploidy_classification():
    def __init__(self, options):
        self.diff = 0.05
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        if hasattr(self, 'classid'):
            self.classid = [str(k) for k in self.classid.split(',')]
        else:
            self.classid = ['class1', 'class2']

    def run(self):
        ancestor_left = base.read_calassfication(self.ancestor_left)
        ancestor_top = base.read_calassfication(self.ancestor_top)
        bkinfo = pd.read_csv(self.blockinfo)
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo = bkinfo[(bkinfo['chr1'].isin(ancestor_left[0].values)) & (
            bkinfo['chr2'].isin(ancestor_top[0].values))]
        bkinfo[self.classid[0]] = 0
        bkinfo[self.classid[1]] = 0
        for name, group in bkinfo.groupby('chr1'):
            d1 = ancestor_left[ancestor_left[0] == name]
            for index1, row1 in group.iterrows():
                a, b = sorted([row1['start1'], row1['end1']])
                a, b = int(a), int(b)
                for index2, row2 in d1.iterrows():
                    c, d = sorted([row2[1], row2[2]])
                    h = len([k for k in range(a, b) if k in range(c, d)])/(b-a)
                    if h >= float(self.diff):
                        bkinfo.loc[index1, self.classid[0]] = row2[4]
        for name, group in bkinfo.groupby('chr2'):
            d2 = ancestor_top[ancestor_top[0] == name]
            for index1, row1 in group.iterrows():
                a, b = sorted([row1['start2'], row1['end2']])
                a, b = int(a), int(b)
                for index2, row2 in d2.iterrows():
                    c, d = sorted([row2[1], row2[2]])
                    h = len([k for k in range(a, b) if k in range(c, d)])/(b-a)
                    if h >= float(self.diff):
                        bkinfo.loc[index1, self.classid[1]] = row2[4]
        bkinfo.to_csv(self.savefile, index=None)
        sys.exit(0)
