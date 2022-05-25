import re
import sys

import numpy as np
import pandas as pd
from Bio import SeqIO

import wgdi.base as base


class ancestral_karyotype_repertoire():
    def __init__(self, options):
        self.gap = 5
        self.direction = 0.01
        self.mark = 'aak1s'
        self.blockinfo_reverse = 'false'
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)

    def run(self):
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        bkinfo = pd.read_csv(self.blockinfo, index_col='id')
        if self.blockinfo_reverse == True or self.blockinfo_reverse.upper() == 'TRUE':
            bkinfo[['chr1', 'chr2']] = bkinfo[['chr2', 'chr1']]
            bkinfo[['block1', 'block2']] = bkinfo[['block2', 'block1']]
        bkinfo = bkinfo.head(5)
        for index, row in bkinfo.iterrows():
            block1, block2 = row['block1'].split('_'), row['block2'].split('_')
            block1, block2 = [int(k) for k in block1], [int(k) for k in block2]
            if int(block1[1])-int(block1[0]) < 0:
                self.direction = -0.01
            for i in range(1, len(block2)):
                if abs(block1[i]-block1[i-1]) == 1 and abs(block2[i]-block2[i-1]) < int(self.gap):
                    gff1_id = gff1[(gff1['chr'] == str(row['chr1'])) & (
                        gff1['order'] == block1[i])].index[0]
                    order = gff1.loc[gff1_id, 'order']
                    gff1_row = gff1.loc[gff1_id, :].copy()
                    for num in range(block2[i-1], block2[i]):
                        order = order + self.direction
                        id = gff2[(gff2['chr'] == str(row['chr2']))
                                  & (gff2['order'] == num)].index[0]
                        gff1_row['order'] = order
                        gff1.loc[id, :] = gff1_row
        gff1.to_csv('out.csv')
        df = gff1.copy()
        df = df.sort_values(by=['chr', 'order'])
        for name, group in df.groupby(['chr']):
            df.loc[group.index, 'order'] = list(range(1, len(group)+1))
            df.loc[group.index, 'newname'] = list(
                [str(self.mark)+str(name)+'g'+str(i).zfill(5) for i in range(1, len(group)+1)])
        df['order'] = df['order'].astype('int')
        df['oldname'] = df.index
        columns = ['chr', 'newname', 'start',
                   'end', 'stand', 'order', 'oldname']
        df[columns].to_csv(self.ancestor_gff, sep="\t",
                           index=False, header=None)

        lens = df.groupby('chr').max()[['end', 'order']]
        lens.to_csv(self.ancestor_lens, sep="\t", header=None)
        ancestor = base.read_calassfication(self.ancestor)
        for index, row in ancestor.iterrows():
            id1 = gff1[(gff1['chr'] == str(row[0])) & (
                gff1['order'] == row[1])].index[0]
            id2 = gff1[(gff1['chr'] == str(row[0])) & (
                gff1['order'] == row[2])].index[0]
            ancestor.loc[index, 1] = df.loc[id1, 'order']
            ancestor.loc[index, 2] = df.loc[id2, 'order']
        ancestor.to_csv(self.ancestor_new, sep="\t", index=False, header=None)
        id_dict = df['newname'].to_dict()
        seqs = []
        for seq_record in SeqIO.parse(self.ancestor_pep, "fasta"):
            if seq_record.id in id_dict:
                seq_record.id = id_dict[seq_record.id]
            else:
                continue
            seq_record.description = ''
            seqs.append(seq_record)
        SeqIO.write(seqs, self.ancestor_pep_new, "fasta")
