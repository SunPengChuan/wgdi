import os
import re
import sys
from io import StringIO

import numpy as np
import pandas as pd
import wgdi.base as base
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import MafftCommandline, MuscleCommandline
from Bio.Phylo.PAML import yn00


class ks():
    def __init__(self, options):
        base_conf = base.config()
        self.pair_pep_file = 'pair.pep'
        self.pair_cds_file = 'pair.cds'
        self.prot_align_file = 'prot.aln'
        self.mrtrans = 'pair.mrtrans'
        self.pair_yn = 'pair.yn'
        self.cds_file = 'cds'
        self.pep_file = 'pep'
        for k, v in base_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)

    def auto_file(self):
        pairs = []
        p = pd.read_csv(self.pairs_file, sep='\n', header=None, nrows=30)
        p = '\n'.join(p[0])
        if 'path length' in p or 'MAXIMUM GAP' in p:
            collinearity = base.read_colinearscan(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif 'MATCH_SIZE' in p or '## Alignment' in p:
            collinearity = base.read_mcscanx(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif '# Alignment' in p:
            collinearity = base.read_coliearity(self.pairs_file)
            pairs = [[v[0], v[2]] for k in collinearity for v in k[1]]
        elif ',' in p:
            collinearity = pd.read_csv(self.pairs_file, header=None)
            pairs = collinearity.values.tolist()
        else:
            collinearity = pd.read_csv(self.pairs_file, header=None, sep='\t')
            pairs = collinearity.values.tolist()
        df = pd.DataFrame(pairs)
        df = df.drop_duplicates()
        df[0] = df[0].astype(str)
        df[1] = df[1].astype(str)
        df.index = df[0]+','+df[1]
        return df

    def run(self):
        path = os.getcwd()
        if not os.path.exists(self.pep_file):
            base.cds_to_pep(os.path.join(path, self.cds_file),
                            os.path.join(path, self.pep_file))
        cds = SeqIO.to_dict(SeqIO.parse(self.cds_file, "fasta"))
        pep = SeqIO.to_dict(SeqIO.parse(self.pep_file, "fasta"))
        df_pairs = self.auto_file()
        if os.path.exists(self.ks_file):
            ks = pd.read_csv(self.ks_file, sep='\t')
            ks = ks.drop_duplicates()
            ks['id'] = ks['id1']+','+ks['id2']
            df_pairs.drop(np.intersect1d(df_pairs.index,
                                         ks['id'].to_numpy()), inplace=True)
            ks_file = open(self.ks_file, 'a+')
        else:
            ks_file = open(self.ks_file, 'w')
            ks_file.write(
                '\t'.join(['id1', 'id2', 'ka_NG86', 'ks_NG86', 'ka_YN00', 'ks_YN00'])+'\n')
        df_pairs = df_pairs[(df_pairs[0].isin(cds.keys())) & (df_pairs[1].isin(
            cds.keys())) & (df_pairs[0].isin(pep.keys())) & (df_pairs[1].isin(pep.keys()))]
        pairs = df_pairs[[0, 1]].to_numpy()
        for k in pairs:
            SeqIO.write([cds[k[0]], cds[k[1]]], self.pair_cds_file, "fasta")
            SeqIO.write([pep[k[0]], pep[k[1]]], self.pair_pep_file, "fasta")
            kaks = self.pair_kaks(k)
            print(k, kaks)
            if kaks == None:
                continue
            ks_file.write('\t'.join([str(i) for i in list(k)+list(kaks)])+'\n')
        ks_file.close()
        for file in (self.pair_pep_file, self.pair_cds_file, self.mrtrans, self.pair_yn, self.prot_align_file, '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'):
            try:
                os.remove(file)
            except OSError:
                pass

    def pair_kaks(self, k):
        self.align()
        pal = self.pal2nal()
        if not pal:
            return []
        kaks = self.run_yn00()
        if kaks == None:
            return []
        kaks_new = [kaks[k[0]][k[1]]['NG86']['dN'], kaks[k[0]][k[1]]['NG86']
                    ['dS'], kaks[k[0]][k[1]]['YN00']['dN'], kaks[k[0]][k[1]]['YN00']['dS']]
        return kaks_new

    def align(self):
        if self.align_software == 'mafft':
            mafft_cline = MafftCommandline(
                cmd=self.mafft_path, input=self.pair_pep_file, auto=True)
            stdout, stderr = mafft_cline()
            align = AlignIO.read(StringIO(stdout), "fasta")
            AlignIO.write(align, self.prot_align_file, "fasta")
        if self.align_software == 'muscle':
            muscle_cline = MuscleCommandline(
                cmd=self.muscle_path, input=self.pair_pep_file, out=self.prot_align_file, seqtype="protein", clwstrict=True)
            stdout, stderr = muscle_cline()

    def pal2nal(self):
        args = ['perl', self.pal2nal_path, self.prot_align_file,
                self.pair_cds_file, '-output paml -nogap', '>'+self.mrtrans]
        command = ' '.join(args)
        try:
            os.system(command)
        except:
            return False
        return True

    def run_yn00(self):
        yn = yn00.Yn00()
        yn.alignment = self.mrtrans
        yn.out_file = self.pair_yn
        yn.set_options(icode=0, commonf3x4=0, weighting=0, verbose=1)
        try:
            run_result = yn.run(command=self.yn00_path)
        except:
            run_result = None
        return run_result
