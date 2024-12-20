import pandas as pd
import wgdi.base as base


class polyploidy_classification:
    def __init__(self, options):
        self.same_protochromosome = False
        self.same_subgenome = False
        for k, v in options:
            setattr(self, str(k), v)
            print(f"{k} = {v}")
        
        self.same_protochromosome = base.str_to_bool(self.same_protochromosome)
        self.same_subgenome = base.str_to_bool(self.same_subgenome)
        
        # Initialize classid with a default value if not provided
        self.classid = [str(k) for k in getattr(self, 'classid', 'class1,class2').split(',')]

    def run(self):
        # Read input files
        ancestor_left = base.read_classification(self.ancestor_left)
        ancestor_top = base.read_classification(self.ancestor_top)
        bkinfo = pd.read_csv(self.blockinfo)

        # Ensure chr1 and chr2 are treated as strings
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)

        # Filter rows where chr1 and chr2 match ancestor values
        bkinfo = bkinfo[bkinfo['chr1'].isin(ancestor_left[0].values) & bkinfo['chr2'].isin(ancestor_top[0].values)]

        # Initialize additional columns
        bkinfo[self.classid[0]] = 0
        bkinfo[self.classid[1]] = 0
        bkinfo[self.classid[0] + '_color'] = ''
        bkinfo[self.classid[1] + '_color'] = ''
        bkinfo['diff'] = 0.0

        # Processing the first classification (ancestor_left vs chr1)
        for name, group in bkinfo.groupby('chr1'):
            d1 = ancestor_left[ancestor_left[0] == name]
            for index1, row1 in group.iterrows():
                a, b = sorted([row1['start1'], row1['end1']])
                a, b = int(a), int(b)
                for index2, row2 in d1.iterrows():
                    c, d = sorted([row2[1], row2[2]])
                    h = len([k for k in range(a, b) if k in range(c, d)]) / (b - a)
                    if h > bkinfo.loc[index1, 'diff']:
                        bkinfo.loc[index1, 'diff'] = float(h)
                        bkinfo.loc[index1, self.classid[0]] = row2[4]
                        bkinfo.loc[index1, self.classid[0] + '_color'] = row2[3]

        # Reset 'diff' and process the second classification (ancestor_top vs chr2)
        bkinfo['diff'] = 0.0
        for name, group in bkinfo.groupby('chr2'):
            d2 = ancestor_top[ancestor_top[0] == name]
            for index1, row1 in group.iterrows():
                a, b = sorted([row1['start2'], row1['end2']])
                a, b = int(a), int(b)
                for index2, row2 in d2.iterrows():
                    c, d = sorted([row2[1], row2[2]])
                    h = len([k for k in range(a, b) if k in range(c, d)]) / (b - a)
                    if h > bkinfo.loc[index1, 'diff']:
                        bkinfo.loc[index1, 'diff'] = float(h)
                        bkinfo.loc[index1, self.classid[1]] = row2[4]
                        bkinfo.loc[index1, self.classid[1] + '_color'] = row2[3]

        # Uncomment if you want to filter rows where both colors match
        if self.same_protochromosome == True:
            bkinfo = bkinfo[bkinfo[self.classid[1] + '_color'] == bkinfo[self.classid[0] + '_color']]
        if self.same_subgenome == True:
            bkinfo = bkinfo[bkinfo[self.classid[1]] == bkinfo[self.classid[0]]]  

        # Save the result to a CSV file
        bkinfo.to_csv(self.savefile, index=False)
