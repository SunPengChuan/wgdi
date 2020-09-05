import fileinput
import os
import re
import sys
from io import StringIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MafftCommandline, MuscleCommandline


class trees():
    def __init__(self, options):
        base_conf = base.config()
        self.position = 'order'
        self.alignfile = ''
        self.cdsfile = ''
        for k, v in base_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)

    def grouping(self, alignment):
        data = []
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)
        cds = SeqIO.to_dict(SeqIO.parse(self.cds_file, "fasta"))
        for index, row in alignment.iterrows():
            ids = []
            for i in range(len(row[:-2])):
                if row[i] is np.nan:
                    continue
                gene_cds = cds[row[i]]
                gene_cds.id = str(int(i)+1)
                ids.append(gene_cds)
            file = str(row['chr'])+'g'+str(row[self.position])
            self.cdsfile = os.path.join(self.dir, file+'.fasta')
            self.alignfile = os.path.join(self.dir, file+'.aln')
            print(self.alignfile)
            if os.path.isfile(self.alignfile):
                continue
            SeqIO.write(ids, self.cdsfile, "fasta")
            self.align()
            self.buildtrees()
            data.append(self.alignfile)
        return data

    def align(self):
        if self.align_software == 'mafft':
            mafft_cline = MafftCommandline(
                cmd=self.mafft_path, input=self.cdsfile, auto=True)
            stdout, stderr = mafft_cline()
            align = AlignIO.read(StringIO(stdout), "fasta")
            AlignIO.write(align, self.alignfile, "fasta")
        if self.align_software == 'muscle':
            muscle_cline = MuscleCommandline(
                cmd=self.muscle_path, input=self.cdsfile, out=self.alignfile, seqtype="nucleo", clwstrict=True)
            stdout, stderr = muscle_cline()

    def buildtrees(self):
        args = [self.iqtree_path, '-s', self.alignfile, '-m', self.model]
        command = ' '.join(args)
        try:
            os.system(command)
        except:
            return False
        return True

    def run(self):
        alignment = pd.read_csv(self.alignment, header=None)
        alignment.replace('.', np.nan, inplace=True)
        gff = base.newgff(self.gff)
        lens = base.newlens(self.lens, self.position)
        gff = gff[gff['chr'].isin(lens.index)]
        print(alignment.head())
        alignment.dropna(thresh=3, inplace=True)
        alignment = pd.merge(
            alignment, gff[['chr', self.position]], left_on=0, right_on=gff.index, how='left')
        data = self.grouping(alignment)
        data = [k+'.treefile' for k in data]
        with open(self.trees_file, 'w') as fout, fileinput.input(data) as fin:
            for line in fin:
                fout.write(line)
