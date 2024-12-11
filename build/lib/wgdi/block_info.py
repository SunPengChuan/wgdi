import numpy as np
import pandas as pd
import wgdi.base as base


class block_info:
    def __init__(self, options):
        self.repeat_number = 20
        self.ks_col = 'ks_NG86'
        self.blast_reverse = False
        for k, v in options:
            setattr(self, str(k), v)
            print(f"{k} = {v}")
        
        self.repeat_number = int(self.repeat_number)
        self.blast_reverse = base.str_to_bool(self.blast_reverse)

    def block_position(self, collinearity, blast, gff1, gff2, ks):
        data = []
        for block in collinearity:
            blk_homo, blk_ks = [], []

            # Skip blocks with missing gene coordinates in GFF files
            if block[1][0][0] not in gff1.index or block[1][0][2] not in gff2.index:
                continue
            
            # Extract chromosome info
            chr1, chr2 = gff1.at[block[1][0][0], 'chr'], gff2.at[block[1][0][2], 'chr']
            
            # Extract start and end positions
            array1, array2 = [float(i[1]) for i in block[1]], [float(i[3]) for i in block[1]]
            start1, end1 = array1[0], array1[-1]
            start2, end2 = array2[0], array2[-1]
            
            block1, block2 = [], []
            for k in block[1]:
                block1.append(int(float(k[1])))
                block2.append(int(float(k[3])))
                
                # Check for KS values
                pair_ks = self.get_ks_value(ks, k)
                blk_ks.append(pair_ks)

                # Retrieve blast homo data
                if k[0]+","+k[2] in blast.index:
                    blk_homo.append(blast.loc[k[0]+","+k[2], [f'homo{i}' for i in range(1, 6)]].values.tolist())
            
            ks_median, ks_average = self.calculate_ks_statistics(blk_ks)
            homo = self.calculate_homo_statistics(blk_homo)

            blkks = '_'.join([str(k) for k in blk_ks])
            block1 = '_'.join([str(k) for k in block1])
            block2 = '_'.join([str(k) for k in block2])
            
            # Calculate tandem ratio
            tandem_ratio = self.tandem_ratio(blast, gff2, block[1])
            
            # Store the results
            data.append([
                block[0], chr1, chr2, start1, end1, start2, end2, block[2], len(block[1]), 
                ks_median, ks_average, *homo, block1, block2, blkks, tandem_ratio
            ])
        
        # Create a DataFrame with the results
        data_df = pd.DataFrame(data, columns=[
            'id', 'chr1', 'chr2', 'start1', 'end1', 'start2', 'end2', 'pvalue', 'length', 
            'ks_median', 'ks_average', 'homo1', 'homo2', 'homo3', 'homo4', 'homo5', 
            'block1', 'block2', 'ks', 'tandem_ratio'
        ])

        # Calculate density
        data_df['density1'] = data_df['length'] / ((data_df['end1'] - data_df['start1']).abs() + 1)
        data_df['density2'] = data_df['length'] / ((data_df['end2'] - data_df['start2']).abs() + 1)

        return data_df

    def get_ks_value(self, ks, k):
        """Return KS value for the given pair of genes."""
        pair = f"{k[0]},{k[2]}"
        if pair in ks.index:
            return ks[pair]
        pair_rev = f"{k[2]},{k[0]}"
        if pair_rev in ks.index:
            return ks[pair_rev]
        return -1

    def calculate_ks_statistics(self, blk_ks):
        """Calculate KS statistics: median and average."""
        ks_arr = [k for k in blk_ks if k >= 0]
        if len(ks_arr) == 0:
            return -1, -1
        ks_median = base.get_median(ks_arr)
        ks_average = sum(ks_arr) / len(ks_arr)
        return ks_median, ks_average

    def calculate_homo_statistics(self, blk_homo):
        """Calculate homo statistics by averaging across all blocks."""
        df = pd.DataFrame(blk_homo)
        homo = df.mean().values if len(df) > 0 else [-1, -1, -1, -1, -1]
        return homo

    def blast_homo(self, blast, gff1, gff2, repeat_number):
        """Assign homo values based on blast data."""
        index = [group.sort_values(by=11, ascending=False)[:repeat_number].index.tolist() for name, group in blast.groupby([0])]
        blast = blast.loc[np.concatenate([k[:repeat_number] for k in index], dtype=object), [0, 1]]
        blast = blast.assign(homo1=np.nan, homo2=np.nan, homo3=np.nan, homo4=np.nan, homo5=np.nan)

        # Assign homo values
        for i in range(1, 6):
            bluenum = i + 5
            redindex = np.concatenate([k[:i] for k in index], dtype=object)
            blueindex = np.concatenate([k[i:bluenum] for k in index], dtype=object)
            grayindex = np.concatenate([k[bluenum:repeat_number] for k in index], dtype=object)
            blast.loc[redindex, f'homo{i}'] = 1
            blast.loc[blueindex, f'homo{i}'] = 0
            blast.loc[grayindex, f'homo{i}'] = -1
        
        blast['chr1_order'] = blast[0].map(gff1['order'])
        blast['chr2_order'] = blast[1].map(gff2['order'])
        return blast

    def tandem_ratio(self, blast, gff2, block):
        """Calculate tandem ratio for a block."""
        block = pd.DataFrame(block)[[0, 2]].rename(columns={0: 'id1', 2: 'id2'})
        block['order2'] = block['id2'].map(gff2['order'])

        # Filter block_blast data
        block_blast = blast[(blast[0].isin(block['id1'].values)) & (blast[1].isin(block['id2'].values))].copy()
        block_blast = pd.merge(block_blast, block, left_on=0, right_on='id1', how='left')
        block_blast['difference'] = (block_blast['chr2_order'] - block_blast['order2']).abs()

        # Filter based on difference and calculate ratio
        block_blast = block_blast[(block_blast['difference'] <= self.repeat_number) & (block_blast['difference'] > 0)]
        return len(block_blast[0].unique()) / len(block) * len(block_blast) / (len(block) + len(block_blast))

    def run(self):
        """Main function to run the analysis."""
        # Initialize required datasets
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)

        # Filter GFF files based on chromosome indices
        gff1 = gff1[gff1['chr'].isin(lens1.index)]
        gff2 = gff2[gff2['chr'].isin(lens2.index)]

        # Load blast data
        blast = base.newblast(self.blast, int(self.score), float(self.evalue), gff1, gff2, self.blast_reverse)
        blast = self.blast_homo(blast, gff1, gff2, self.repeat_number)
        blast.index = blast[0] + ',' + blast[1]

        # Get collinearity data
        collinearity = self.auto_file(gff1, gff2)

        # Load ks data if necessary
        ks = pd.Series([]) if self.ks == 'none' or self.ks == '' or not hasattr(self, 'ks') else base.read_ks(self.ks, self.ks_col)

        # Get the block position data
        data = self.block_position(collinearity, blast, gff1, gff2, ks)
        data['class1'] = 0
        data['class2'] = 0

        # Save results
        data.to_csv(self.savefile, index=None)

    def auto_file(self, gff1, gff2):
        """Auto-detect and read collinearity file."""
        with open(self.collinearity) as f:
            p = ' '.join(f.readlines()[0:30])
        
        # Handle different file formats
        if 'path length' in p or 'MAXIMUM GAP' in p:
            return base.read_colinearscan(self.collinearity)
        elif 'MATCH_SIZE' in p or '## Alignment' in p:
            return self.process_mcscanx(gff1, gff2)
        elif '# Alignment' in p:
            return base.read_collinearity(self.collinearity)
        elif '###' in p:
            return self.process_jcvi(gff1, gff2)

    def process_mcscanx(self, gff1, gff2):
        """Process MCScanX format collinearity data."""
        col = base.read_mcscanx(self.collinearity)
        collinearity = []
        for block in col:
            newblock = [k for k in block[1] if k[0] in gff1.index and k[2] in gff2.index]
            if newblock:
                for k in newblock:
                    k[1], k[3] = gff1.at[k[0], 'order'], gff2.at[k[2], 'order']
                collinearity.append([block[0], newblock, block[2]])
        return collinearity

    def process_jcvi(self, gff1, gff2):
        """Process JCVI format collinearity data."""
        col = base.read_jcvi(self.collinearity)
        collinearity = []
        for block in col:
            newblock = [k for k in block[1] if k[0] in gff1.index and k[2] in gff2.index]
            if newblock:
                for k in newblock:
                    k[1], k[3] = gff1.at[k[0], 'order'], gff2.at[k[2], 'order']
                collinearity.append([block[0], newblock, block[2]])
        return collinearity
