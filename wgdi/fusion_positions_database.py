import pandas as pd
import os
from Bio import SeqIO

class fusion_positions_database:
    def __init__(self, options):
        self.NumGeneSetsPerSide = 100
        for k, v in options:
            setattr(self, k, v)
            print(f'{k} = {v}')
        self.NumGeneSetsPerSide = int(self.NumGeneSetsPerSide)

    def run(self):
        # Load and remove duplicates from data
        gff = pd.read_csv(self.gff, sep="\t", header=None, dtype={0: str, 5: int}).drop_duplicates()
        pep = SeqIO.to_dict(SeqIO.parse(self.pep, "fasta"))
        df = pd.read_csv(self.fusion_positions, sep="\t", header=None, dtype={0: str, 1: int}).drop_duplicates()
        
        # Load ancestral sequence file if it exists
        seqs = SeqIO.to_dict(SeqIO.parse(self.ancestor_pep, "fasta")) if os.path.exists(self.ancestor_pep) else {}

        sf_gff, sf_lens = [], []

        # Process fusion positions
        for _, row in df.iterrows():
            newchr = row[2]
            newgff = gff[(gff[0] == row[0]) & 
                         (gff[5] >= row[1] - self.NumGeneSetsPerSide) & 
                         (gff[5] < row[1] + self.NumGeneSetsPerSide)].copy()
            newgff['id'] = [f"{newchr}s{str(row[0]).zfill(2)}g{str(i).zfill(3)}" for i in range(1, len(newgff) + 1)]
            sf_position = row[1] - newgff.iloc[0, 5]
            sf_lens.append([newchr, sf_position, len(newgff)])
            
            # For each gene in the filtered GFF region
            for _, gff_row in newgff.iterrows():
                if gff_row[1] in pep and gff_row['id'] not in seqs:
                    gene = pep[gff_row[1]][:]
                    gene.id, gene.description = gff_row['id'], ''
                    seqs[gff_row['id']] = gene
                    # Collect data for the final GFF output
                    sf_gff.append([gff_row['id'], newchr, sf_position, gff_row[2], gff_row[3], gff_row[4], gff_row[1]])

        # Write sequences to FASTA file
        SeqIO.write(seqs.values(), self.ancestor_pep, 'fasta')

        # Save filtered GFF data
        if sf_gff:
            print('run')
            sf_gff = pd.DataFrame(sf_gff)
            sf_gff.rename(columns={3: 'start', 4: 'end', 5: 'strand'}, inplace=True)
            sf_gff['order'] = sf_gff[0].str[-3:].astype(int)
            sf_gff[[1, 0, 'start', 'end', 'strand', 'order', 6]].to_csv(self.ancestor_gff, sep="\t", mode='a', index=False, header=None)
            sf_lens = pd.DataFrame(sf_lens).drop_duplicates()
            sf_lens.to_csv(self.ancestor_lens, sep="\t", mode='a', index=False, header=None)

            # Generate ancestral sequence data
            ancestor = []
            for _, row in sf_lens.iterrows():
                ancestor.append([row[0], 1, row[1], 'red', 1])
                ancestor.append([row[0], row[1] + 1, row[2], 'blue', 1])
            pd.DataFrame(ancestor).to_csv(self.ancestor_file, sep="\t", mode='a', index=False, header=None)

        # Remove duplicates from the output files
        for file in [self.ancestor_gff, self.ancestor_lens, self.ancestor_file]:
            df = pd.read_csv(file, header=None).drop_duplicates().to_csv(file, index=False, header=None)
