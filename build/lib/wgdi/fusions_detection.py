import pandas as pd
from tabulate import tabulate

class fusions_detection:
    def __init__(self, options):
        self.min_genes_per_side = 5
        self.density = 0.3
        for k, v in options:
            setattr(self, k, v)
            print(f"{k} = {v}")
        self.min_genes_per_side = int(self.min_genes_per_side)
        self.density = float(self.density)

    def run(self):
        # Load the ancestor file and process the positions
        ancestor = pd.read_csv(self.ancestor, sep='\t', header=None)
        position = ancestor.groupby(0)[2].unique().apply(pd.Series)
        bkinfo = pd.read_csv(self.blockinfo)
        newbkinfo = bkinfo.head(0)
        
        # Iterate over each row in the position dataframe
        for index, row in position.iterrows():
            # Filter the bkinfo dataframe based on chr2 and density
            filtered_group = bkinfo[(bkinfo['chr2'] == index) & (bkinfo['density2'] >= self.density)].copy()
            # Split the block2 column and stack the resulting series
            df = filtered_group['block2'].str.split('_', expand=True).stack().astype(int)
            # Count the number of genes greater and less than the current position
            filtered_group['greater'] = (df > row[0]).groupby(level=0).sum()
            filtered_group['less'] = (df < row[0]).groupby(level=0).sum()
            # Filter the group based on the minimum number of genes per side
            filtered_group = filtered_group[(filtered_group['greater'] >= self.min_genes_per_side) & (filtered_group['less'] >= self.min_genes_per_side)]
            # Concatenate the filtered group with the newbkinfo dataframe
            newbkinfo = pd.concat([newbkinfo, filtered_group])

        # Get and print the shared fusion positions
        newbkinfo.to_csv(self.filtered_blockinfo, header=True, index=False)

        non_overlap_counts = newbkinfo.groupby('chr2').apply(self.count_non_overlapping)
        data = [(chr2, count) for chr2, count in non_overlap_counts.items()]
        print("\nThe following are the shared fusion breakpoints and counts:")
        print(tabulate(data, headers=["Fusion Breakpoint", "Count"], tablefmt="github"))

    def count_non_overlapping(self, group):
        if len(group) == 1:
            return 1
        grouped = group.groupby('chr1')
        total_count = 0
        for chr1, chr_group in grouped:
            chr_group = chr_group.sort_values(by='start1').reset_index(drop=True)
            count = 0
            current_end = -1 
            for _, row in chr_group.iterrows():
                start1, end1 = row['start1'], row['end1']
                if start1 > current_end:
                    count += 1
                    current_end = end1 
            total_count += count
        return total_count