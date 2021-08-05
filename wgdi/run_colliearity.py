import gc
import os
import re
import shutil
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

import wgdi.base as base
import wgdi.collinearity as improvedcollinearity


class mycollinearity():
    def __init__(self, options):
        self.repeat_number = 10
        self.multiple = 1
        self.score = 100
        self.evalue = 1e-5
        self.blast_reverse = 'False'
        self.over_windows = 5
        self.options = options
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.position = 'order'
        if hasattr(self, 'grading'):
            self.grading = [int(k) for k in self.grading.split(',')]
        else:
            self.grading = [50, 40, 25]
        if hasattr(self, 'process'):
            self.process = int(self.process)
        else:
            self.process = 4

    def deal_blast(self, blast, gff1, gff2, rednum, repeat_number):
        blast['grading'] = 0
        bluenum = 4+rednum
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
            total.append([chr1, chr2, group])
        del blast, group
        gc.collect()
        n = int(np.ceil(len(total) / float(self.process)))
        result, data = '', []
        try:
            pool = Pool(self.process)
            for i in range(0, len(total), n):
                data.append(pool.apply_async(self.single_pool, args=(
                    total[i:i + n], gff1, gff2, lens1, lens2)))
            pool.close()
            pool.join()
        except:
            pool.terminate()
        for k in data:
            result += k.get()
        result = re.split('\n', result)
        fout = open(self.savefile, 'w')
        num = 1
        for line in result:
            if re.match(r"# Alignment", line):
                s = '# Alignment ' + str(num) + ':'
                fout.writelines(s+line.split(':')[1]+'\n')
                num += 1
                continue
            if len(line) > 0:
                fout.writelines(line+'\n')
        fout.close()
        sys.exit(0)

    def single_pool(self, group, gff1, gff2, lens1, lens2):
        text = ''
        for bk in group:
            chr1, chr2 = str(bk[0]), str(bk[1])
            print('runing ', chr1, 'vs', chr2)
            df = pd.DataFrame(np.zeros((lens1[chr1], lens2[chr2])),columns=range(1,lens2[chr2]+1),index=range(1,lens1[chr1]+1))
            for index, row in bk[2].iterrows():
                df.loc[row['loc1'], row['loc2']] = row['grading']
            df = df.loc[:, df.sum(axis=0) != 0]
            df = df.loc[df.sum(axis=1) != 0, :]
            collinearity = improvedcollinearity.collinearity(
                self.options, df)
            data = collinearity.run()
            gf1, gf2 = gff1[gff1['chr'] == chr1], gff2[gff2['chr'] == chr2]
            blocks, evalues, socres = data
            n = 1
            for i in range(len(blocks)):
                if len(blocks[i]) < self.over_windows:
                    continue
                if blocks[i][1][0]-blocks[i][0][0] > 0:
                    mark = 'plus'
                else:
                    mark = 'minus'
                text += '# Alignment '+str(n)+': score='+str(socres[i])+' pvalue=' + str(
                    evalues[i])+' N=' + str(len(blocks[i])) + ' '+str(chr1)+'&'+str(chr2) + ' ' + mark+'\n'
                n += 1
                for k in blocks[i]:
                    name1 = gf1[gf1['order'] == k[0]].index[0]
                    name2 = gf2[gf2['order'] == k[1]].index[0]
                    if gf1.loc[name1, 'stand'] == gf2.loc[name2, 'stand']:
                        order = '1'
                    else:
                        order = '-1'
                    s = ' '.join([name1, str(k[0]), name2, str(k[1]), order])
                    text += s+'\n'
        return text
