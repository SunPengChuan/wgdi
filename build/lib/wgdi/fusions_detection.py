import pandas as pd


class fusions_detection:
    def __init__(self, options):
        self.min_genes_per_side = 5
        for k, v in options:
            setattr(self, k, v)
            print(f"{k} = {v}")
        self.min_genes_per_side = int(self.min_genes_per_side)

    def run(self):
        # Load the ancestor file and process the positions
        ancestor = pd.read_csv(self.ancestor, sep='\t', header=None)
        position = ancestor.groupby(0)[2].unique().apply(pd.Series)
        bkinfo = pd.read_csv(self.blockinfo)
        bkinfo['break_length1'] = 0
        bkinfo['break_length2'] = 0
        newbkinfo = bkinfo.head(0)
        
        # Iterate over each row in the position dataframe
        for index, row in position.iterrows():
            group = bkinfo[bkinfo['chr2'] == index].copy()
            group.loc[:, 'break_length1'] = group['start2'] - row[0]
            group.loc[:, 'break_length2'] = group['end2'] - row[0]
            group = group[
                ((group['break_length1'] > self.min_genes_per_side) & (group['break_length2'] < -self.min_genes_per_side)) |
                ((group['break_length1'] < -self.min_genes_per_side) & (group['break_length2'] > self.min_genes_per_side))
            ]
            df = group['block2'].str.split('_', expand=True).stack().astype(int)
            group['greater'] = (df > row[0]).groupby(level=0).sum()
            group['less'] = (df < row[0]).groupby(level=0).sum()
            group = group[(group['greater'] >= 5) & (group['less'] >= 5)]
            newbkinfo = pd.concat([newbkinfo, group])

        # Get and print the shared fusion positions
        shared_fusion_positions = newbkinfo['chr2'].unique()
        print("\nThe following are the shared fusion breakpoints:")
        print(", ".join(shared_fusion_positions))
        newbkinfo.to_csv(self.filtered_blockinfo, header=True, index=False)
