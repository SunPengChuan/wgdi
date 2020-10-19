import os
import re
import shutil
import sys

import numpy as np
import pandas as pd

import wgdi.base as base
import wgdi.collinearity as improvedcollinearity


class mycollinearity():
    def __init__(self, options):
        self.repeat_number = 20
        self.multiple = 1
        self.score = 100
        self.evalue = 1e-5
        self.blast_reverse = 'False'
        self.over_windows = 5
        self.options = options
        for k, v in options:
            setattr(self, str(k), v)
            # print(str(k), ' = ', v)
        self.position = 'order'
        if hasattr(self, 'grading'):
            self.grading = [int(k) for k in self.grading.split(',')]
        else:
            self.grading = [50, 40, 25]
        if os.path.exists(self.dir):
            shutil.rmtree(self.dir)
        os.makedirs(self.dir)
        self.number = 1

    def deal_blast(self, blast, gff1, gff2, rednum, repeat_number):
        blast['grading'] = 0
        bluenum = 2+rednum
        index = [group.sort_values(by=[11], ascending=[False])[:repeat_number].index.tolist()
                 for name, group in blast.groupby([0])]
        reddata = np.array([k[:rednum] for k in index], dtype=object)
        bluedata = np.array([k[rednum:bluenum] for k in index], dtype=object)
        graydata = np.array([k[bluenum:repeat_number]
                             for k in index], dtype=object)
        if len(reddata):
            redindex = np.concatenate(reddata)
        else:
            redindex = []
        if len(bluedata):
            blueindex = np.concatenate(bluedata)
        else:
            blueindex = []
        if len(graydata):
            grayindex = np.concatenate(graydata)
        else:
            grayindex = []
        blast.loc[redindex, 'grading'] = self.grading[0]
        blast.loc[blueindex, 'grading'] = self.grading[1]
        blast.loc[grayindex, 'grading'] = self.grading[2]
        return blast[blast['grading'] > 0]

    def write_block(self, fout, data, chr1, chr2, gff1, gff2):
        gff1, gff2 = gff1[gff1['chr'] == chr1], gff2[gff2['chr'] == chr2]
        fout = open(fout, 'w')
        blocks, evalues, socres = data
        for i in range(len(blocks)):
            if len(blocks[i]) < self.over_windows:
                continue
            if blocks[i][1][0]-blocks[i][0][0] > 0:
                mark = 'plus'
            else:
                mark = 'minus'
            fout.write('# Alignment '+str(self.number)+': score='+str(socres[i])+' e_value=' + str(
                evalues[i])+' N=' + str(len(blocks[i])) + ' '+str(chr1)+'&'+str(chr2) + ' ' + mark+'\n')
            self.number += 1
            for k in blocks[i]:
                name1 = gff1[gff1['order'] == k[0]].index[0]
                name2 = gff2[gff2['order'] == k[1]].index[0]
                if gff1.loc[name1, 'stand'] == gff2.loc[name2, 'stand']:
                    order = '1'
                else:
                    order = '-1'
                s = ' '.join([name1, str(k[0]), name2, str(k[1]), order])
                fout.write(s+'\n')

    def run(self):
        lens1 = base.newlens(self.lens1, 'order')
        lens2 = base.newlens(self.lens2, 'order')
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gff1 = gff1[gff1['chr'].isin(lens1.index)]
        gff2 = gff2[gff2['chr'].isin(lens2.index)]
        blast = base.newblast(self.blast, int(self.score), float(
            self.evalue), gff1, gff2, self.blast_reverse)
        blast = self.deal_blast(blast, gff1, gff2, int(
            self.multiple), int(self.repeat_number))
        blast['loc1'] = blast[0].map(gff1.loc[:, self.position])
        blast['loc2'] = blast[1].map(gff2.loc[:, self.position])
        blast['chr1'] = blast[0].map(gff1.loc[:, 'chr'])
        blast['chr2'] = blast[1].map(gff2.loc[:, 'chr'])
        total = []
        for (chr1, chr2), group in blast.groupby(['chr1', 'chr2']):
            df = pd.DataFrame(np.zeros((lens1[chr1], lens2[chr2])))
            for index, row in group.iterrows():
                df.loc[row['loc1'], row['loc2']] = row['grading']
            df = df.loc[:, df.sum(axis=0) != 0]
            df = df.loc[df.sum(axis=1) != 0, :]
            collinearity = improvedcollinearity.collinearity(self.options, df)
            data = collinearity.run()
            fp = self.dir+'/'+str(chr1)+'.vs.'+str(chr2)+'.blk'
            gf1, gf2 = gff1[gff1['chr'] == chr1], gff2[gff2['chr'] == chr1]
            self.write_block(fp, data, chr1, chr2, gff1, gff2)
        args = ['cat', self.dir+'/*.blk', '>', self.savefile]
        command = ' '.join([str(k) for k in args])
        os.system(command)
        shutil.rmtree(self.dir)
        sys.exit(0)
