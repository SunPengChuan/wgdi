import os
import re
import sys
from io import StringIO

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
            colinearity = base.read_colinearscan(self.pairs_file)
            pairs = [[v[0], v[2]] for k in colinearity for v in k[0]]
        elif 'MATCH_SIZE' in p or '## Alignment' in p:
            colinearity = base.read_mcscanx(self.pairs_file)
            pairs = [k[1:] for k in colinearity]
        elif ',' in p:
            colinearity = pd.read_csv(self.pairs_file, header=None)
            pairs = colinearity.values.tolist()
        else:
            colinearity = pd.read_csv(self.pairs_file, header=None, sep='\t')
            pairs = colinearity.values.tolist()
        df = pd.DataFrame(pairs)
        df = df.drop_duplicates()
        return df.values.tolist()

    def run(self):
        path = os.getcwd()
        if self.pep_file == 'pep':
            base.cds_to_pep(os.path.join(path, self.cds_file),
                            os.path.join(path, self.pep_file))
        cds = SeqIO.to_dict(SeqIO.parse(self.cds_file, "fasta"))
        pep = SeqIO.to_dict(SeqIO.parse(self.pep_file, "fasta"))
        pairs = self.auto_file()
        ks_file = open(self.ks_file, 'w')
        ks_file.write(
            '\t'.join(['id1', 'id2', 'ka_NG86', 'ks_NG86', 'ka_YN00', 'ks_YN00'])+'\n')
        for k in pairs[0:20]:
            self.pair = str(k[0]+','+str(k[1]))
            if k[0] in cds.keys() and k[1] in cds.keys() and k[0] in pep.keys() and k[1] in pep.keys():
                cds[k[0]].id = cds[k[0]].id.replace('.', '_')
                cds[k[1]].id = cds[k[1]].id.replace('.', '_')
                pep[k[0]].id = pep[k[0]].id.replace('.', '_')
                pep[k[1]].id = pep[k[1]].id.replace('.', '_')
                SeqIO.write([cds[k[0]], cds[k[1]]],
                            self.pair_cds_file, "fasta")
                SeqIO.write([pep[k[0]], pep[k[1]]],
                            self.pair_pep_file, "fasta")
            else:
                continue
            print(k)
            kaks = self.pair_kaks(k)
            print(kaks)
            if kaks == None:
                continue
            ks_file.write('\t'.join([str(i) for i in k+kaks])+'\n')
        for file in (self.pair_pep_file, self.pair_cds_file, self.pep_file, self.mrtrans, self.pair_yn, 
            self.prot_align_file, '2YN.dN', '2YN.dS', '2YN.t', 'rst', 'rst1', 'yn00.ctl', 'rub'):
            os.remove(file)

    def pair_kaks(self, k):
        k[0], k[1] = k[0].replace('.', '_'), k[1].replace('.', '_')
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
                mafft_path=self.mafft_path, input=self.pair_pep_file, auto=True)
            stdout, stderr = mafft_cline()
            self.prot_align_file = AlignIO.read(StringIO(stdout), "fasta")
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
