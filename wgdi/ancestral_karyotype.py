import re
import sys

import numpy as np
import pandas as pd
from Bio import SeqIO

import wgdi.base as base


class ancestral_karyotype():
    def __init__(self, options):
        self.mark = 'aak'
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)

    def run(self):
        gff = base.newgff(self.gff)
        ancestor = base.read_calassfication(self.ancestor)
        gff = gff[gff['chr'].isin(ancestor[0].values.tolist())]
        newgff = gff.copy()
        data = []
        chr_arr = ancestor[3].drop_duplicates().to_list()
        chr_dict = dict(zip(chr_arr, range(1, len(chr_arr)+1)))
        dict1, dict2 = {}, {}
        for (cla, color), group in ancestor.groupby([4, 3], sort=[False, False]):
            num = chr_dict[color] + len(chr_arr)*(int(cla)-1)
            for index, row in group.iterrows():
                index1 = gff[(gff['chr'] == row[0]) & (
                    gff['order'] >= row[1]) & (gff['order'] <= row[2])].index
                newgff.loc[index1, 'chr'] = str(num)
                for k in index1:
                    data.append(newgff.loc[k, :].values.tolist()+[k])
            dict1[str(num)] = cla
            dict2[str(num)] = color
        df = pd.DataFrame(data)
        pep = SeqIO.to_dict(SeqIO.parse(self.pep_file, "fasta"))
        df = df[df[6].isin(pep.keys())]
        for name, group in df.groupby([0]):
            df.loc[group.index, 'order'] = list(range(1, len(group)+1))
            df.loc[group.index, 'newname'] = list(
                [str(self.mark)+str(name)+'g'+str(i).zfill(5) for i in range(1, len(group)+1)])
        df['order'] = df['order'].astype('int')
        df = df[[0, 'newname', 1, 2, 3, 'order', 6]]
        df = df.sort_values(by=[0, 'order'])
        df.to_csv(self.ancestor_gff, sep="\t", index=False, header=None)
        lens = df.groupby(0).max()[[2, 'order']]
        lens.to_csv(self.ancestor_lens, sep="\t", header=None)
        lens[1] = 1
        lens['color'] = lens.index.map(dict2)
        lens['class'] = lens.index.map(dict1)
        lens[[1, 'order', 'color', 'class']].to_csv(
            self.ancestor_file, sep="\t", header=None)
        id_dict = df.set_index(6).to_dict()['newname']
        seqs = []
        for seq_record in SeqIO.parse(self.pep_file, "fasta"):
            if seq_record.id in id_dict:
                seq_record.id = id_dict[seq_record.id]
            else:
                continue
            seqs.append(seq_record)
        SeqIO.write(seqs, self.ancestor_pep, "fasta")
