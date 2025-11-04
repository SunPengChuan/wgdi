import os
import shutil
from io import StringIO

import numpy as np
import pandas as pd
from Bio import AlignIO, Seq, SeqIO, SeqRecord
import subprocess

import wgdi.base as base


class trees():
    def __init__(self, options):
        base_conf = base.config()
        self.position = 'order'
        self.alignfile = ''
        self.align_trimming = ''
        self.trimming = 'trimal'
        self.threads = '1'
        self.minimum = 4
        self.tree_software = 'iqtree'
        self.delete_detail = True
        for k, v in base_conf:
            setattr(self, str(k), v)
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        if hasattr(self, 'codon_position'):
            self.codon_position = [
                int(k)-1 for k in self.codon_position.split(',')]
        else:
            self.codon_position = [0, 1, 2]
        self.delete_detail = base.str_to_bool(self.delete_detail)

    def grouping(self, alignment):
        data = []
        indexs = []
        if not os.path.exists(self.dir):
            os.makedirs(self.dir)
        sequence = SeqIO.to_dict(SeqIO.parse(self.sequence_file, "fasta"))
        if hasattr(self, 'cds_file'):
            seq_cds = SeqIO.to_dict(SeqIO.parse(self.cds_file, "fasta"))
        for index, row in alignment.iterrows():
            file = base.gen_md5_id(str(row.values))
            self.sequencefile = os.path.join(self.dir, file+'.fasta')
            self.alignfile = os.path.join(self.dir, file+'.aln')
            self.align_trimming = self.alignfile+'.trimming'
            self.treefile = os.path.join(self.dir, file+'.aln.treefile')
            if os.path.isfile(self.treefile) and os.path.isfile(self.alignfile):
                data.append(self.treefile)
                indexs.append(index)
                continue
            ids = []
            ids_cds = []
            for i in range(len(row)):
                if type(row[i]) == float and np.isnan(row[i]):
                    continue
                gene_sequence = sequence[row[i]]
                gene_sequence.id = str(int(i)+1)
                gene_sequence.description = ''
                ids.append(gene_sequence)
            SeqIO.write(ids, self.sequencefile, "fasta")
            self.align()
            if hasattr(self, 'cds_file'):
                self.seqcdsfile = os.path.join(self.dir, file+'.cds.fasta')
                for i in range(len(row)):
                    if type(row[i]) == float and np.isnan(row[i]):
                        continue
                    gene_cds = seq_cds[row[i]]
                    gene_cds.id = str(int(i)+1)
                    ids_cds.append(gene_cds)
                SeqIO.write(ids_cds, self.seqcdsfile, "fasta")
                self.pal2nal()
                self.codon()
            if self.trimming.upper() == 'TRIMAL':
                self.trimal()
            if self.trimming.upper() == 'DIVVIER':
                self.divvier()
            self.buildtrees()
            if os.path.isfile(self.treefile):
                data.append(self.treefile)
        return data

    def codon(self):
        if self.codon_position == [0, 1, 2]:
            shutil.move(self.alignfile+'.mrtrans', self.alignfile)
            return True
        records = list(SeqIO.parse(self.alignfile+'.mrtrans', 'fasta'))
        if len(records) == 0:
            return False
        newrecords = []
        def final_list(test_list, x, y): return [
            test_list[i+j] for i in range(0, len(test_list), x) for j in y]
        for k in records:
            if len(k.seq) % 3 > 0:
                return False
            seq = final_list(k.seq, 3, self.codon_position)
            k.seq = ''.join(seq)
            newrecords.append(SeqRecord.SeqRecord(
                Seq.Seq(k.seq), id=k.id, description=''))
        SeqIO.write(newrecords, self.alignfile, 'fasta')
        return True

    def pal2nal(self):
        args = ['perl', self.pal2nal_path, self.alignfile,
                self.seqcdsfile, '-output fasta', '>'+self.alignfile+'.mrtrans']
        command = ' '.join(args)
        try:
            os.system(command)
        except:
            return False
        return True

    def align(self):
        if self.align_software == 'mafft':
            try:
                command = [self.mafft_path,'--quiet', self.sequencefile, '>', self.alignfile]
                subprocess.run(" ".join(command), shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while running MAFFT: {e}")

        if self.align_software == 'muscle':
            try:
                command = [self.muscle_path,'-align', self.sequencefile, '-output', self.alignfile, '-quiet']
                subprocess.run(" ".join(command), shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while running Muscle: {e}")

    def trimal(self):
        args = [self.trimal_path, '-in', self.alignfile,
                '-out', self.align_trimming, '-automated1']
        command = ' '.join(args)
        try:
            os.system(command)
        except:
            return False
        return True

    def divvier(self):
        args = [self.divvier_path, '-mincol', '4', '-divvygap', self.alignfile]
        command = ' '.join(args)
        try:
            os.system(command)
            os.rename(self.alignfile+'.divvy.fas', self.align_trimming)
        except:
            return False
        return True

    def buildtrees(self):
        try:
            if self.tree_software.upper() == 'IQTREE':
                args = [self.iqtree_path, '-s', self.align_trimming,
                        '-m', self.model, '-T', self.threads, '--quiet']
                command = ' '.join(args)
                os.system(command)
                os.rename(self.align_trimming+'.treefile', self.treefile)
            elif self.tree_software.upper() == 'FASTTREE':
                args = [self.fasttree_path,
                        self.align_trimming, '>', self.treefile]
                command = ' '.join(args)
                os.system(command)
        except:
            return False
        if self.delete_detail == True:
            for file in (self.sequencefile, self.align_trimming+'.bionj', self.align_trimming+'.iqtree', self.align_trimming+'.ckp.gz',
                         self.align_trimming+'.log', self.align_trimming+'.mldist', self.align_trimming+'.model.gz'):
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
        df[0].to_csv(self.trees_file, index=None, sep='\t', header=False)
        print("done")