import os
import re
import shutil
import sys

import numpy as np
import pandas as pd
import wgdi.base as base


class colinearscan():
    def __init__(self, options):
        self.repeat_number = 20
        self.score = 100
        self.evalue = 1e-5
        self.position = 'order'
        self.blast_reverse = 'False'
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

    def deal_blast(self, blast, gff1, gff2, repeat_number):
        index1 = [group.sort_values(by=11, ascending=False)[:repeat_number].index.tolist()
                  for name, group in blast.groupby([0])]
        index1 = np.concatenate(np.array(index1))
        index2 = [group.sort_values(by=11, ascending=False)[:repeat_number].index.tolist()
                  for name, group in blast.groupby([1])]
        index2 = np.concatenate(np.array(index2))
        index = np.intersect1d(index1, index2)
        blast = blast.loc[index, [0, 1]]
        gff1 = gff1[['chr', 'stand', 'order']]
        gff2 = gff2[['chr', 'stand', 'order']]
        gff1.columns = ['chr1', 'stand1', 'loc1']
        gff2.columns = ['chr2', 'stand2', 'loc2']
        blast = pd.merge(blast, gff1, left_on=0, right_on=gff1.index)
        blast = pd.merge(blast, gff2, left_on=1, right_on=gff2.index)
        blast.replace({'+': '1', '-': '-1'}, inplace=True)
        return blast

    def rewriteblock(self, blast, file, fout):
        num = 0
        fout = open(fout, 'w')
        collinearity = base.read_colinearscan(file)
        for block in collinearity:
            num += 1
            if block[1][-1][0]+','+block[1][-1][2] not in blast.index:
                block[1] = block[1][:-1]
            fout.write('the '+str(num)+'th path length ' +
                       str(len(block[1]))+'\n')
            for k in block[1]:
                s = ' '.join([str(i) for i in k])
                fout.write(s+'\n')
            fout.write('>LOCALE p-value :'+block[2]+'\n\n')

    def run(self):
        lens1 = base.newlens(self.lens1, 'order')
        lens2 = base.newlens(self.lens2, 'order')
        lens1 = lens1[lens1 > 4]
        lens2 = lens2[lens2 > 4]
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gff1 = gff1[gff1['chr'].isin(lens1.index)]
        gff2 = gff2[gff2['chr'].isin(lens2.index)]
        blast = base.newblast(self.blast, int(self.score), float(
            self.evalue), gff1, gff2, self.blast_reverse)
        df = self.deal_blast(blast, gff1, gff2, int(self.repeat_number))
        for (chr1, chr2), group in df.groupby(['chr1', 'chr2']):
            group = group.sort_values(by=['loc1', 'loc2'])
            dir1 = './'+self.dir+'/pair/'+str(chr1)+'.vs.'+str(chr2)+'.pair'
            dir2 = './'+self.dir+'/block/'+str(chr1)+'.vs.'+str(chr2)+'.blk'
            group[[0, 'stand1', 'loc1', 1, 'stand2', 'loc2']].to_csv(
                dir1, sep=' ', index=None, header=None)
            args = ['blockscan', '-chr1len', lens1[str(chr1)], '-chr2len', lens2[str(
                chr2)], '-mg1', self.mg[0], '-mg2', self.mg[1], dir1, '>'+dir2]
            command = ' '.join([str(k) for k in args])
            os.system(command)
        args = ['cat', self.dir+'/block/*.blk', '>', self.dir+'.block.old.txt']
        command = ' '.join([str(k) for k in args])
        os.system(command)
        self.rewriteblock(blast, self.dir+'.block.old.txt',
                          self.dir+'.block.txt')
        shutil.rmtree(self.dir)
        sys.exit(0)
