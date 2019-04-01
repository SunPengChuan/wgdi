import os
import re
import shutil
import sys

import numpy as np
import pandas as pd
import wgdi.base as base


class colinearscan():
    def __init__(self, options):
        self.repnum = 20
        self.score = 200
        self.evalue = 1e-5
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
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gff1 = gff1[gff1['chr'].isin(lens1.index)]
        gff2 = gff2[gff2['chr'].isin(lens2.index)]
        blast = base.newblast(self.blast, int(self.score),
                              float(self.evalue), gff1, gff2)
        df = self.deal_blast(blast, gff1, gff2, int(self.repnum))
        for (chr1, chr2), group in df.groupby(['chr1', 'chr2']):
            group = group.sort_values(by=['loc1', 'loc2'])
            # print(group.head())
            dir1 = './'+self.dir+'/pair/'+str(chr1)+'.vs.'+str(chr2)+'.pair'
            dir2 = './'+self.dir+'/block/'+str(chr1)+'.vs.'+str(chr2)+'.blk'
            group[[0, 'stand1', 'loc1', 1, 'stand2', 'loc2']].to_csv(
                dir1, sep=' ', index=None, header=None)
            args = ['blockscan', '-chr1len', lens1[str(chr1)], '-chr2len', lens2[str(
                chr2)], '-mg1', self.mg[0], '-mg2', self.mg[1], dir1, '>'+dir2]
            command = ' '.join([str(k) for k in args])
            # os.system(command)
        args = ['cat', self.dir+'/block/*.blk', '>', self.dir+'.block.txt']
        command = ' '.join([str(k) for k in args])
        # os.system(command)
        self.rewriteblock(sys.argv[1],sys.argv[2])

    def deal_blast(self, blast, gff1, gff2, repnum):
        index = [group[:repnum].index.tolist()
                 for name, group in blast.groupby([0])]
        index = np.concatenate(np.array(index))
        blast = blast.loc[index, [0, 1]]
        gff1 = gff1[['chr','stand','order']]
        gff2 = gff2[['chr','stand','order']]
        gff1.columns = ['chr1', 'stand1', 'loc1']
        gff2.columns = ['chr2', 'stand2', 'loc2']
        blast = pd.merge(blast, gff1, left_on=0,right_on=gff1.index)
        blast = pd.merge(blast, gff2, left_on=1,right_on=gff2.index)
        blast.replace({'+': '1', '-': '-1'}, inplace=True)
        return blast

    def rewriteblock(self,file,fout):
        num = 1
        fout = open(fout,'w')
        with open(file) as f:
            for line in f.readlines():
                if re.match(r"the", line):
                    line = line.replace(re.search('\d+',line).group(),str(num),1)
                    num+=1
                fout.write(line)
        f.close()



# for block in colinearity:
#     if len(block[0]) ==0:
#         continue
#     for k in block[0]:
#         ges_1.append(k[0])
# blast = blast[(blast[0].isin(set(ges_1)))]
# homopairs = {}
# for name, group in blast.groupby([0])[1]:
#     newgroup = group.values.tolist()[:repnum]
#     for i, el in enumerate(newgroup, start=1):
#         if i <= dupnum:
#             homopairs[name+","+el] = 1
#         elif i <= hitnum:
#             homopairs[name+","+el] = 0
#         else:
#             homopairs[name+","+el] = -1
# pos = []
# blocks=[]
# for block in colinearity:
#     a, b, blk_ks, homo = [], [], [], 0
#     if len(block[0]) == 0:
#         continue
#     for k in block[0]:
#         if (k[0] not in gff_1.index) or (k[2] not in gff_2.index):
#             continue
#         if k[0]+","+k[2] in ks.index:
#             blk_ks.append(ks.at[k[0]+","+k[2], 3])
#         if k[0]+","+k[2] in homopairs.keys():
#             homo += homopairs[k[0]+","+k[2]]
#     print(homo)
#     homo = homo/len(block[0])
#     print(homo, len(block[0]))
#     strand = sum([int(k[4]) for k in block[0]])
#     array1 = [float(k[1]) for k in block[0]]
#     array2 = [float(k[3]) for k in block[0]]
#     chr1, chr2 = gff_1.at[block[0][0][0], 0], gff_2.at[block[0][0][2], 0]
#     y1, x1, y2, x2 = min(array1), min(array2), max(array1), max(array2)
#     pos.append([chr1, chr2, x1, x2, y1, y2, len(block[0]),
#                 strand, base.get_median(blk_ks), homo])
#     blocks.append(block)
# df = pd.DataFrame(pos, columns=['chr1', 'chr2', 'x1', 'x2',
#                                 'y1', 'y2', 'length', 'strand', 'ks', 'homo'])
# df.insert(0, 'id', range(1, len(df)+1))
# df.to_csv('block_info.new.csv')