import numpy as np
import pandas as pd

import wgdi.base as base


class karyotype_mapping:
    def __init__(self, options):
        # Initialize default attributes
        self.blast_reverse = False
        self.blockinfo_reverse = False
        self.position = 'order'
        self.block_length = 5
        self.limit_length = 5
        self.repeat_number = 10
        self.score = 100
        self.evalue = 1e-5

        # Update attributes with provided keyword arguments and print them
        for k, v in options:
            setattr(self, k, v)
            print(f"{k} = {v}")
        
        self.blast_reverse = base.str_to_bool(self.blast_reverse)
        self.blockinfo_reverse = base.str_to_bool(self.blockinfo_reverse)
        self.limit_length = int(self.limit_length)

    def karyotype_left(self, pairs, ancestor, gff1, gff2):
        # Loop through each row in ancestor to set color and classification in gff1
        for _, row in ancestor.iterrows():
            loc_min, loc_max = sorted([row[1], row[2]])
            index1 = gff1[(gff1['chr'] == row[0]) &
                          (gff1['order'] >= loc_min) &
                          (gff1['order'] <= loc_max)].index
            gff1.loc[index1, ['color', 'classification']] = row[3], row[4]

        # Merge pairs with gff1 and update gff2 with color and classification
        data = pd.merge(pairs, gff1, left_on=0, right_index=True, how='left')
        data.drop_duplicates(subset=[1], inplace=True)
        data.set_index(1, inplace=True)
        gff2.loc[data.index, ['color', 'classification']] = data[['color', 'classification']]
        return gff2

    def karyotype_top(self, pairs, ancestor, gff1, gff2):
        # Loop through each row in ancestor to set color and classification in gff2
        for _, row in ancestor.iterrows():
            loc_min, loc_max = sorted([row[1], row[2]])
            index1 = gff2[(gff2['chr'] == row[0]) &
                          (gff2['order'] >= loc_min) &
                          (gff2['order'] <= loc_max)].index
            gff2.loc[index1, ['color', 'classification']] = row[3], row[4]

        # Merge pairs with gff2 and update gff1 with color and classification
        data = pd.merge(pairs, gff2, left_on=1, right_index=True, how='left')
        data.drop_duplicates(subset=[0], inplace=True)
        data.set_index(0, inplace=True)
        gff1.loc[data.index, ['color', 'classification']] = data[['color', 'classification']]
        return gff1

    def karyotype_map(self, gff, lens):
        # Filter gff based on lens index and non-null color
        gff = gff[gff['chr'].isin(lens.index) & gff['color'].notnull()]
        ancestor = []
        # Group by chromosome and process each group to create ancestor records
        for chr, group in gff.groupby('chr'):
            color, class_id, arr = '', 1, []
            for _, row in group.iterrows():
                if color ==  row['color'] and class_id == row['classification']:
                    arr.append(row['order'])
                else:
                    if len(arr) >= self.limit_length:
                        ancestor.append([chr, min(arr), max(arr), color, class_id, len(arr)])
                    color, class_id = row['color'], row['classification']
                    arr = []
                    if len(ancestor) >= 1 and color == ancestor[-1][3] and class_id == ancestor[-1][4] and chr == ancestor[-1][0]:
                        arr.append(ancestor[-1][1])
                        arr += np.random.randint(ancestor[-1][1], ancestor[-1][2], size=ancestor[-1][5]-1).tolist()
                        ancestor.pop()
                    arr.append(row['order'])
            if len(arr) >= self.limit_length:
                ancestor.append([chr, min(arr), max(arr), color, class_id, len(arr)])

        ancestor = pd.DataFrame(ancestor)
        # Adjust min and max positions for each chromosome group
        for chr, group in ancestor.groupby(0):
            ancestor.loc[group.index[0], 1] = 1
            ancestor.loc[group.index[-1], 2] = lens[chr]
        ancestor[4] = ancestor[4].astype(int)
        return ancestor[[0, 1, 2, 3, 4, 5]]

    def colinear_gene_pairs(self, bkinfo, gff1, gff2):
        gff1 = gff1.reset_index()
        gff2 = gff2.reset_index()
        
        gff1_indexed = gff1.set_index(['chr', 'order'])
        gff2_indexed = gff2.set_index(['chr', 'order'])
        
        data = []
        for _, row in bkinfo.iterrows():
            b1 = list(map(int, row['block1'].split('_')))
            b2 = list(map(int, row['block2'].split('_')))

            for order1, order2 in zip(b1, b2):
                a = gff1_indexed.loc[(row['chr1'], order1), 1]
                b = gff2_indexed.loc[(row['chr2'], order2), 1]
                data.append([a, b])
        return pd.DataFrame(data)
    
    def new_ancestor(self, ancestor, gff1, gff2, blast):
        # Iterate through ancestor rows to adjust positions based on neighboring rows
        for i in range(1, len(ancestor)):
            if ancestor.iloc[i, 0] == ancestor.iloc[i-1, 0]:
                area = ancestor.iloc[i, 1] - ancestor.iloc[i-1, 2]
                if area <= 5 and ancestor.iloc[i, 5] > 200:
                    ancestor.iloc[i-1, 2] = ancestor.iloc[i, 1] - 1
                else:
                    if ancestor.iloc[i, 5] > 200:
                        index1 = gff1[(gff1['chr'] == ancestor.iloc[i, 0]) &
                                    (gff1['order'] >= ancestor.iloc[i-1, 2]+1) &
                                    (gff1['order'] <= ancestor.iloc[i, 1]-1)].index
                    else:
                        index1 = gff1[(gff1['chr'] == ancestor.iloc[i, 0]) &
                                    (gff1['order'] >= ancestor.iloc[i-1, 2]+1) &
                                    (gff1['order'] <= ancestor.iloc[i, 2])].index
                    index2 = gff2[gff2['color'] == ancestor.iloc[i-1, 3]].index
                    index3 = gff2[gff2['color'] == ancestor.iloc[i, 3]].index

                    newblast1 = blast[(blast[0].isin(index1)) & (blast[1].isin(index2))]
                    newblast2 = blast[(blast[0].isin(index1)) & (blast[1].isin(index3))]

                    if len(newblast1) >= len(newblast2) and ancestor.iloc[i, 5] > 200:
                        ancestor.iloc[i-1, 2] = ancestor.iloc[i, 1] - 1
                    elif ancestor.iloc[i, 5] < 200:
                        ancestor.iloc[i, 3] = ancestor.iloc[i-1, 3]
                    else:
                        ancestor.iloc[i, 1] = ancestor.iloc[i-1, 2] + 1
        
        ancestor['group'] = (ancestor[0].shift(1) != ancestor[0]) | (ancestor[3].shift(1) != ancestor[3]) | (ancestor[4].shift(1) != ancestor[4])
        ancestor['group'] = ancestor['group'].cumsum()

        # 使用groupby进行分组，并应用聚合函数
        result = ancestor.groupby('group').agg({
            0: 'first',
            1: 'min',
            2: 'max',
            3: 'first',
            4: 'first',
        }).reset_index(drop=True)

        return result

    def run(self):
        # Read and process block information
        bkinfo = pd.read_csv(self.blockinfo, index_col='id')
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        if self.blockinfo_reverse == True:
            bkinfo[['chr1', 'chr2']] =  bkinfo[['chr2', 'chr1']]
            bkinfo[['block1', 'block2']] =  bkinfo[['block2', 'block1']]
        bkinfo = bkinfo[bkinfo['length'] > int(self.block_length)]

        # Read GFF and lens data
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        lens = base.newlens(self.the_other_lens, self.position)
        blast = base.newblast(self.blast, int(self.score), float(self.evalue), gff1, gff2, self.blast_reverse)
        # blast.drop_duplicates(subset=[0], keep='first', inplace=True)

        # Find colinear gene pairs
        pairs = self.colinear_gene_pairs(bkinfo, gff1, gff2)

        # Depending on available attributes, call either karyotype_top or karyotype_left
        if hasattr(self, 'ancestor_top'):
            ancestor = base.read_classification(self.ancestor_top)
            data = self.karyotype_top(pairs, ancestor, gff1, gff2)
        elif hasattr(self, 'ancestor_left'):
            ancestor = base.read_classification(self.ancestor_left)
            data = self.karyotype_left(pairs, ancestor, gff1, gff2)
        else:
            print('Missing ancestor file.')
            exit(0)

        # Map the data and create the final ancestor file
        the_other_ancestor_file = self.karyotype_map(data, lens)
        # print(the_other_ancestor_file)
        the_other_ancestor_file = self.new_ancestor(the_other_ancestor_file, gff1, gff2, blast)
        the_other_ancestor_file.to_csv(self.the_other_ancestor_file, sep='\t', header=False, index=False)