import os
import re
import shutil
import sys

import numpy as np
import pandas as pd


class colinearscan():
    def __init__(self, options):
        self.partial = False
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        if hasattr(self, 'mg'):
            self.mg = self.mg.split(',')
        else:
            self.mg = [50, 50]
        if os.path.exists(self.dir):
            shutil.rmtree(self.dir)
        os.makedirs(self.dir)
        os.makedirs(self.dir+'/pair/')
        os.makedirs(self.dir+'/block/')

    def run(self):
        gff_1 = pd.read_csv(self.gff1, sep="\t", header=None)
        gff_2 = pd.read_csv(self.gff2, sep="\t", header=None)
        gff_1[0] = gff_1[0].astype('str')
        gff_2[0] = gff_2[0].astype('str')
        lens_1 = pd.read_csv(self.lens1, sep="\t", header=None, index_col=0)
        lens_2 = pd.read_csv(self.lens2, sep="\t", header=None, index_col=0)
        lens_1.index = lens_1.index.astype('str')
        lens_2.index = lens_2.index.astype('str')
        blast = pd.read_csv(self.blast, sep="\t", header=None)
        df = self.deal_blast(blast, gff_1, gff_2)
        dict={}
        for (chr1, chr2), group in df.groupby(['chr1', 'chr2']):
            if self.partial == True:
                if chr1+'\t'+chr2 in dict:
                    continue
            if str(chr1) not in lens_1.index or str(chr2) not in lens_2.index:
                continue
            dict[chr1+'\t'+chr2], dict[chr2+'\t'+chr1] = 1,1
            group = group.drop_duplicates()
            group = group.sort_values(by=['loc1', 'loc2'])
            dir1 = './'+self.dir+'/pair/'+str(chr1)+'.vs.'+str(chr2)+'.pair'
            dir2 = './'+self.dir+'/block/'+str(chr1)+'.vs.'+str(chr2)+'.blk'
            group[['id1', 'stand1', 'loc1', 'id2', 'stand2', 'loc2']].to_csv(
                dir1, sep=' ', index=None, header=None)
            args = ['blockscan', '-chr1len', lens_1.at[str(chr1), 1], '-chr2len', lens_2.at[str(
                chr2), 1], '-mg1', self.mg[0], '-mg2', self.mg[1], dir1, '>'+dir2]
            command = ' '.join([str(k) for k in args])
            os.system(command)
        args = ['cat', self.dir+'/block/*.blk', '>', self.dir+'.block.txt']
        command = ' '.join([str(k) for k in args])
        os.system(command)

    def deal_blast(self, blast, gff1, gff2):
        blast = blast[(blast[11] >= float(self.score)) & (
            blast[10] <= float(self.evalue)) & (blast[1] != blast[0])]
        blast.drop_duplicates(subset=[0,1],keep='first',inplace=True)
        n = 5
        array = []
        for name, group in blast.groupby([0]):
            index = group.index[n:]
            array += index.tolist()
        index = list(set(array))
        blast.drop(index, inplace=True)
        blast = blast[[0, 1]]
        blast.columns = ['id1', 'id2']
        gff1 = gff1[[1, 0, 4, 5]]
        gff1.columns = ['id1', 'chr1', 'stand1', 'loc1']
        gff2 = gff2[[1, 0, 4, 5]]
        gff2.columns = ['id2', 'chr2', 'stand2', 'loc2']
        df = pd.merge(blast, gff1, on='id1')
        df = pd.merge(df, gff2, on='id2')
        df.replace('+', '1', inplace=True)
        df.replace('-', '-1', inplace=True)
        return df
