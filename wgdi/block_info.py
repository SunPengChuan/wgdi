import numpy as np
import pandas as pd
import wgdi.base as base


class block_info():
    def __init__(self, options):
        self.repnum = 20
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)

    def block_position(self, colinearity, blast, gff1, gff2, ks):
        data = []
        for block in colinearity:
            blk_homo, blk_ks = [],  []
            if block[1][0][0] not in gff1.index or block[1][0][2] not in gff2.index:
                continue
            chr1, chr2 = gff1.loc[block[1][0][0],
                                  'chr'], gff2.loc[block[1][0][2], 'chr']
            array1, array2 = [float(i[1]) for i in block[1]], [
                float(i[3]) for i in block[1]]
            start1, end1 = array1[0], array1[-1]
            start2, end2 = array2[0], array2[-1]
            for k in block[1]:
                if k[0]+","+k[2] not in blast.index:
                    continue
                blk_homo.append(
                    blast.loc[k[0]+","+k[2], ['homo'+str(i) for i in range(1, 6)]].values.tolist())
                if k[0]+","+k[2] in ks.index:
                    pair_ks = ks.at[str(k[0])+","+str(k[2]), 3]
                    blk_ks.append(pair_ks)
                else:
                    blk_ks.append(0)
            blkks = ','.join([str(k) for k in blk_ks])
            df = pd.DataFrame(blk_homo)
            homo = df.mean().values
            if len(homo) == 0:
                continue
            data.append([block[0], chr1, chr2, start1, end1, start2, end2, block[2], len(
                block[1]), base.get_median(blk_ks), homo[0], homo[1], homo[2], homo[3], homo[4], blkks])
        data = pd.DataFrame(data, columns=['id', 'chr1', 'chr2', 'start1', 'end1', 'start2', 'end2',
                                           'pvalue', 'length', 'ks_median', 'homo1', 'homo2', 'homo3', 'homo4', 'homo5', 'ks'])
        data.to_csv(self.savefile, index=None)
        return data

    def blast_homo(self, blast, gff1, gff2, repnum):
        index = [group[:repnum].index.tolist()
                 for name, group in blast.groupby([0])]
        blast = blast.loc[np.concatenate(
            np.array([k[:repnum] for k in index])), [0, 1]]
        blast = blast.assign(homo1=np.nan, homo2=np.nan,
                             homo3=np.nan, homo4=np.nan, homo5=np.nan)
        for i in range(1, 6):
            bluenum = i+5
            redindex = np.concatenate(np.array([k[:i] for k in index]))
            blueindex = np.concatenate(np.array([k[i:bluenum] for k in index]))
            grayindex = np.concatenate(
                np.array([k[bluenum:repnum] for k in index]))
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
        blast = base.newblast(self.blast, int(self.score),
                              float(self.evalue), gff1, gff2)
        blast = self.blast_homo(blast, gff1, gff2, int(self.repnum))
        blast.index = blast[0]+','+blast[1]
        colinearity = base.read_colinearscan(self.colinearity)
        ks = base.read_ks(self.ks)
        data = self.block_position(colinearity, blast, gff1, gff2, ks)
