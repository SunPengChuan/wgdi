import numpy as np
import pandas as pd

import wgdi.base as base


class block_info():
    def __init__(self, options):
        self.repeat_number = 20
        self.ks_col = 'ks_NG86'
        self.blast_reverse = 'False'
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)

    def block_position(self, collinearity, blast, gff1, gff2, ks):
        data = []
        for block in collinearity:
            blk_homo, blk_ks = [],  []
            if block[1][0][0] not in gff1.index or block[1][0][2] not in gff2.index:
                continue
            chr1, chr2 = gff1.loc[block[1][0][0],
                                  'chr'], gff2.loc[block[1][0][2], 'chr']
            array1, array2 = [float(i[1]) for i in block[1]], [
                float(i[3]) for i in block[1]]
            start1, end1 = array1[0], array1[-1]
            start2, end2 = array2[0], array2[-1]
            block1, block2 = [], []
            for k in block[1]:
                block1.append(int(float(k[1])))
                block2.append(int(float(k[3])))
                if k[0]+","+k[2] in ks.index:
                    pair_ks = ks[str(k[0])+","+str(k[2])]
                    blk_ks.append(pair_ks)
                elif k[2]+","+k[0] in ks.index:
                    pair_ks = ks[str(k[2])+","+str(k[0])]
                    blk_ks.append(pair_ks)
                else:
                    blk_ks.append(-1)
                if k[0]+","+k[2] not in blast.index:
                    continue
                blk_homo.append(
                    blast.loc[k[0]+","+k[2], ['homo'+str(i) for i in range(1, 6)]].values.tolist())
            ks_arr = [k for k in blk_ks if k >= 0]
            if len(ks_arr) == 0:
                ks_median = -1
                ks_average = -1
            else:
                arr_ks = [k for k in blk_ks if k >= 0]
                ks_median = base.get_median(arr_ks)
                ks_average = sum(arr_ks)/len(arr_ks)
            df = pd.DataFrame(blk_homo)
            homo = df.mean().values
            if len(homo) == 0:
                homo = [-1, -1, -1, -1, -1]
            blkks = '_'.join([str(k) for k in blk_ks])
            block1 = '_'.join([str(k) for k in block1])
            block2 = '_'.join([str(k) for k in block2])
            blkks = '_' + blkks
            data.append([block[0], chr1, chr2, start1, end1, start2, end2, block[2], len(
                block[1]), ks_median, ks_average, homo[0], homo[1], homo[2], homo[3], homo[4], block1, block2, blkks])
        data = pd.DataFrame(data, columns=['id', 'chr1', 'chr2', 'start1', 'end1', 'start2', 'end2',
                                           'pvalue', 'length', 'ks_median', 'ks_average', 'homo1', 'homo2', 'homo3',
                                           'homo4', 'homo5', 'block1', 'block2', 'ks'])
        data['density1'] = data['length'] / \
            ((data['end1']-data['start1']).abs()+1)
        data['density2'] = data['length'] / \
            ((data['end2']-data['start2']).abs()+1)
        return data

    def blast_homo(self, blast, gff1, gff2, repeat_number):
        index = [group.sort_values(by=11, ascending=False)[:repeat_number].index.tolist()
                 for name, group in blast.groupby([0])]
        blast = blast.loc[np.concatenate(
            np.array([k[:repeat_number] for k in index], dtype=object)), [0, 1]]
        blast = blast.assign(homo1=np.nan, homo2=np.nan,
                             homo3=np.nan, homo4=np.nan, homo5=np.nan)
        for i in range(1, 6):
            bluenum = i+5
            redindex = np.concatenate(
                np.array([k[:i] for k in index], dtype=object))
            blueindex = np.concatenate(
                np.array([k[i:bluenum] for k in index], dtype=object))
            grayindex = np.concatenate(
                np.array([k[bluenum:repeat_number] for k in index], dtype=object))
            blast.loc[redindex, 'homo'+str(i)] = 1
            blast.loc[blueindex, 'homo'+str(i)] = 0
            blast.loc[grayindex, 'homo'+str(i)] = -1
        return blast

    def run(self):
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        gff1 = gff1[gff1['chr'].isin(lens1.index)]
        gff2 = gff2[gff2['chr'].isin(lens2.index)]
        blast = base.newblast(self.blast, int(self.score), float(
            self.evalue), gff1, gff2, self.blast_reverse)
        blast = self.blast_homo(blast, gff1, gff2, int(self.repeat_number))
        blast.index = blast[0]+','+blast[1]
        collinearity = self.auto_file(gff1, gff2)
        ks = base.read_ks(self.ks, self.ks_col)
        data = self.block_position(collinearity, blast, gff1, gff2, ks)
        data['class1'] = 0
        data['class2'] = 0
        data.to_csv(self.savefile, index=None)

    def auto_file(self, gff1, gff2):
        p = pd.read_csv(self.collinearity, sep='\n', header=None, nrows=30)
        p = '\n'.join(p[0])
        if 'path length' in p or 'MAXIMUM GAP' in p:
            collinearity = base.read_colinearscan(self.collinearity)
        elif 'MATCH_SIZE' in p or '## Alignment' in p:
            col = base.read_mcscanx(self.collinearity)
            collinearity = []
            for block in col:
                newblock = []
                for k in block[1]:
                    if k[0] not in gff1.index or k[2] not in gff2.index:
                        continue
                    k[1], k[3] = gff1.loc[k[0], 'order'], gff2.loc[k[2], 'order']
                    newblock.append(k)
                if len(newblock) == 0:
                    continue
                collinearity.append([block[0], newblock, block[2]])
        elif '# Alignment' in p:
            collinearity = base.read_coliearity(self.collinearity)
        elif '###' in p:
            col = base.read_jcvi(self.collinearity)
            collinearity = []
            for block in col:
                newblock = []
                for k in block[1]:
                    if k[0] not in gff1.index or k[2] not in gff2.index:
                        continue
                    k[1], k[3] = gff1.loc[k[0], 'order'], gff2.loc[k[2], 'order']
                    newblock.append(k)
                if len(newblock) == 0:
                    continue
                collinearity.append([block[0], newblock, block[2]])
        return collinearity
