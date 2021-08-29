import fileinput
import os
import re
import sys
from io import StringIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MafftCommandline, MuscleCommandline

import wgdi.base as base


class trees():
    def __init__(self, options):
        base_conf = base.config()
        self.position = 'order'
        self.alignfile = ''
        self.cdsfile = ''
        self.trimming = 'trimal'
        self.minimum = 3
        self.delete_detail = 'true'
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
            file = base.gen_md5_id(str(row.values))
            self.cdsfile = os.path.join(self.dir, file+'.fasta')
            self.alignfile = os.path.join(self.dir, file+'.aln')
            self.treefile = os.path.join(self.dir, file+'.aln.treefile')
            if os.path.isfile(self.treefile):
                data.append(self.treefile)
                continue
            ids = []
            for i in range(len(row)):
                if type(row[i]) == float and np.isnan(row[i]):
                    continue
                print(row[i])
                gene_cds = cds[row[i]]
                gene_cds.id = str(int(i)+1)
                ids.append(gene_cds)

            SeqIO.write(ids, self.cdsfile, "fasta")
            self.align()
            if self.trimming.upper() == 'TRIMAL':
                self.trimal()
            self.buildtrees()
            data.append(self.treefile)
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

    def trimal(self):
        args = [self.trimal_path, '-in', self.alignfile,
                '-out', self.alignfile, '-automated1']
        command = ' '.join(args)
        try:
            os.system(command)
        except:
            return False
        return True

    def buildtrees(self):
        args = [self.iqtree_path, '-s', self.alignfile, '-m', self.model]
        command = ' '.join(args)
        try:
            os.system(command)
        except:
            return False
        if self.delete_detail.upper() == 'TRUE':
            for file in (self.cdsfile, self.alignfile, self.alignfile+'.bionj', self.alignfile+'.iqtree',
                         self.alignfile+'.log', self.alignfile+'.mldist', self.alignfile+'.model.gz'):
                try:
                    os.remove(file)
                except OSError:
                    pass
        return True

    def run(self):
        alignment = pd.read_csv(self.alignment, header=None)
        alignment.replace('.', np.nan, inplace=True)
        alignment.dropna(thresh=int(self.minimum), inplace=True)
        if hasattr(self, 'gff') and hasattr(self, 'lens'):
            gff = base.newgff(self.gff)
            lens = base.newlens(self.lens, self.position)
            alignment = pd.merge(
                alignment, gff[['chr', self.position]], left_on=0, right_on=gff.index, how='left')
            alignment.dropna(subset=['chr', 'order'], inplace=True)
            alignment['order'] = alignment['order'].astype(int)
            alignment = alignment[alignment['chr'].isin(lens.index)]
            alignment.drop(alignment.columns[-2:], axis=1, inplace=True)
        data = self.grouping(alignment)
        fout = open(self.trees_file, 'w')
        fout.close()
        for i in range(0, len(data), 100):
            trees = ' '.join([str(k) for k in data[i:i+100]])
            args = ['cat', trees, '>>', self.trees_file]
            command = ' '.join([str(k) for k in args])
            os.system(command)
        df = pd.read_csv(self.trees_file, header=None, sep='\t')
        df[1] = data
        df[0].to_csv(self.trees_file, index=None, sep='\t', header=False)
