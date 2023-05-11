import gc
import re
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
        self.over_windows = int(self.over_windows)

    def deal_blast(self, blast, rednum, repeat_number):
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
        blast = self.deal_blast(blast, int(
            self.multiple), int(self.repeat_number))
        blast['loc1'] = blast[0].map(gff1[self.position])
        blast['loc2'] = blast[1].map(gff2[self.position])
        blast['chr1'] = blast[0].map(gff1['chr'])
        blast['chr2'] = blast[1].map(gff2['chr'])
        print('The filtered homologous gene pairs are '+str(len(blast))+'.\n')
        if len(len(blast))<1:
            print('Stoped! \n\nIt may be that the id1 and id2 in the BLAST file do not match with (gff1,lens1) and (gff2, lens2).')
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
            points = bk[2][['loc1', 'loc2', 'grading']].sort_values(
                by=['loc1', 'loc2'], ascending=[True, True])
            collinearity = improvedcollinearity.collinearity(
                self.options, points)
            data = collinearity.run()
            if len(data)==0:
                continue
            gf1, gf2 = gff1[gff1['chr'] == chr1], gff2[gff2['chr'] == chr2]
            gf1, gf2 = gf1.reset_index().set_index(
                'order')[[1, 'stand']], gf2.reset_index().set_index('order')[[1, 'stand']]
            n = 1
            for (block, evalue, socre) in data:
                if len(block) < self.over_windows:
                    continue
                block['name1'] = block['loc1'].map(gf1[1])
                block['name2'] = block['loc2'].map(gf2[1])
                block['stand1'] = block['loc1'].map(gf1['stand'])
                block['stand2'] = block['loc2'].map(gf2['stand'])
                block['stand'] = np.where(
                    block['stand1'] == block['stand2'], '1', '-1')
                block['text'] = block.apply(lambda x: str(x['name1'])+' '+str(
                    x['loc1'])+' '+str(x['name2'])+' '+str(x['loc2'])+' '+x['stand']+'\n', axis=1)
                a, b = block['loc2'].head(2).values
                if a < b:
                    mark = 'plus'
                else:
                    mark = 'minus'
                text += '# Alignment '+str(n)+': score='+str(socre)+' pvalue=' + str(
                    evalue)+' N=' + str(len(block)) + ' '+str(chr1)+'&'+str(chr2) + ' ' + mark+'\n'
                text += ''.join(block['text'].values)
                n += 1
        return text
