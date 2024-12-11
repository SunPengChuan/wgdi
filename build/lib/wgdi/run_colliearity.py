import gc
import re
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

import wgdi.base as base
import wgdi.collinearity as improvedcollinearity


class mycollinearity():
    def __init__(self, options):
        # Initialize parameters with default values
        self.repeat_number = 10
        self.multiple = 1
        self.score = 100
        self.evalue = 1e-5
        self.blast_reverse = False
        self.over_gap  = 5
        self.comparison = 'genomes'
        self.options = options

        for k, v in options:
            setattr(self, str(k), v)
            print(f"{str(k)} = {v}")
        self.position = 'order'
        # Parse grading values
        if hasattr(self, 'grading'):
            self.grading = [int(k) for k in self.grading.split(',')]
        else:
            self.grading = [50, 40, 25]
        # Ensure process is an integer
        if hasattr(self, 'process'):
            self.process = int(self.process)
        else:
            self.process = 4
        self.over_gap  = int(self.over_gap )
        base.str_to_bool(self.blast_reverse)

    def deal_blast_for_chromosomes(self, blast, rednum, repeat_number):
        bluenum = rednum
        blast = blast.sort_values(by=[0, 11], ascending=[True, False])
        def assign_grading(group):
            group['cumcount'] = group.groupby(1).cumcount()
            group = group[group['cumcount'] <= repeat_number]
            group['grading'] = pd.cut(
                group['cumcount'],
                bins=[-1, 0, bluenum, repeat_number],
                labels=self.grading,
                right=True
            )
            return group
        newblast = blast.groupby(['chr1', 'chr2']).apply(assign_grading).reset_index(drop=True)
        newblast['grading'] = newblast['grading'].astype(int)
        return newblast[newblast['grading'] > 0]
    
    def deal_blast_for_genomes(self, blast, rednum, repeat_number):
        # Initialize the grading column
        blast['grading'] = 0
        
        # Define the blue number as the sum of rednum and the predefined constant
        bluenum = 4 + rednum
        
        # Get the indices for each group by sorting the 11th column in descending order
        index = [group.sort_values(by=[11], ascending=[False])[:repeat_number].index.tolist()
                for name, group in blast.groupby([0])]
        
        # Split the indices into red, blue, and gray groups
        reddata = np.array([k[:rednum] for k in index], dtype=object)
        bluedata = np.array([k[rednum:bluenum] for k in index], dtype=object)
        graydata = np.array([k[bluenum:repeat_number] for k in index], dtype=object)
        
        # Concatenate the results into flat lists
        redindex = np.concatenate(reddata) if reddata.size else []
        blueindex = np.concatenate(bluedata) if bluedata.size else []
        grayindex = np.concatenate(graydata) if graydata.size else []

        # Update the grading column based on the group indices
        blast.loc[redindex, 'grading'] = self.grading[0]
        blast.loc[blueindex, 'grading'] = self.grading[1]
        blast.loc[grayindex, 'grading'] = self.grading[2]

        # Return only the rows with non-zero grading
        return blast[blast['grading'] > 0]

    def run(self):
        # Read and process lens files
        lens1 = base.newlens(self.lens1, 'order')
        lens2 = base.newlens(self.lens2, 'order')
        # Read and process gff files
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        # Filter gff data based on lens indices
        gff1 = gff1[gff1['chr'].isin(lens1.index)]
        gff2 = gff2[gff2['chr'].isin(lens2.index)]
        # Process blast data

        blast = base.newblast(self.blast, int(self.score), float(self.evalue),gff1, gff2, self.blast_reverse)

        # Map positions and chromosome information
        blast['loc1'] = blast[0].map(gff1[self.position])
        blast['loc2'] = blast[1].map(gff2[self.position])
        blast['chr1'] = blast[0].map(gff1['chr'])
        blast['chr2'] = blast[1].map(gff2['chr'])
        # Apply blast filtering and grading
        if self.comparison.lower() == 'genomes':
            blast = self.deal_blast_for_genomes(blast, int(self.multiple), int(self.repeat_number))
        if self.comparison.lower() == 'chromosomes':
            blast = self.deal_blast_for_chromosomes(blast, int(self.multiple), int(self.repeat_number))
        print(f"The filtered homologous gene pairs are {len(blast)}.\n")
        if len(blast) < 1:
            print("Stopped!\n\nIt may be that the id1 and id2 in the BLAST file do not match with (gff1, lens1) and (gff2, lens2).")
            sys.exit(1)
        # Group blast data by 'chr1' and 'chr2'
        total = []
        for (chr1, chr2), group in blast.groupby(['chr1', 'chr2']):
            total.append([chr1, chr2, group])
        del blast, group
        gc.collect()
        # Determine chunk size for multiprocessing
        n = int(np.ceil(len(total) / float(self.process)))
        result, data = '', []
        try:
            # Initialize multiprocessing Pool
            pool = Pool(self.process)
            for i in range(0, len(total), n):
                # Apply single_pool function asynchronously
                data.append(pool.apply_async(
                    self.single_pool, args=(total[i:i + n], gff1, gff2, lens1, lens2)
                ))
            pool.close()
            pool.join()
        except:
            pool.terminate()
        for k in data:
            # Collect results from async tasks
            text = k.get()
            if text:
                result += text
        # Write final output to file
        result = re.split('\n', result)
        fout = open(self.savefile, 'w')
        num = 1
        for line in result:
            if re.match(r"# Alignment", line):
                # Replace alignment number
                s = f'# Alignment {num}:'
                fout.write(s + line.split(':')[1] + '\n')
                num += 1
                continue
            if len(line) > 0:
                fout.write(line + '\n')
        fout.close()
        sys.exit(0)

    def single_pool(self, group, gff1, gff2, lens1, lens2):
        text = ''
        for bk in group:
            chr1, chr2 = str(bk[0]), str(bk[1])
            print(f'Running {chr1} vs {chr2}')
            # Extract and sort points
            points = bk[2][['loc1', 'loc2', 'grading']].sort_values(
                by=['loc1', 'loc2'], ascending=[True, True]
            )
            # Initialize collinearity analysis
            collinearity = improvedcollinearity.collinearity(
                self.options, points)
            data = collinearity.run()
            if not data:
                continue
            # Extract gene information
            gf1 = gff1[gff1['chr'] == chr1].reset_index().set_index('order')[[1, 'strand']]
            gf2 = gff2[gff2['chr'] == chr2].reset_index().set_index('order')[[1, 'strand']]
            n = 1
            for block, evalue, score in data:
                if len(block) < self.over_gap:
                    continue
                # Map gene names and strands
                block['name1'] = block['loc1'].map(gf1[1])
                block['name2'] = block['loc2'].map(gf2[1])
                block['strand1'] = block['loc1'].map(gf1['strand'])
                block['strand2'] = block['loc2'].map(gf2['strand'])
                block['strand'] = np.where(
                    block['strand1'] == block['strand2'], '1', '-1'
                )
                # Prepare text output
                block['text'] = block.apply(
                    lambda x: f"{x['name1']} {x['loc1']} {x['name2']} {x['loc2']} {x['strand']}\n",
                    axis=1
                )
                # Determine alignment mark
                a, b = block['loc2'].head(2).values
                mark = 'plus' if a < b else 'minus'
                # Append alignment information
                text += f'# Alignment {n}: score={score} pvalue={evalue} N={len(block)} {chr1}&{chr2} {mark}\n'
                text += ''.join(block['text'].values)
                n += 1
        return text