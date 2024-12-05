import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base

class align_dotplot:
    def __init__(self, options):
        # Default values
        self.position = 'order'
        self.figsize = 'default'
        self.classid = 'class1'

        # Initialize from options
        for k, v in options:
            setattr(self, str(k), v)
            print(f'{k} = {v}')
        
        self.ks_area = [float(k) for k in getattr(self, 'ks_area', '-1,3').split(',')]
        self.colors = [str(k) for k in getattr(self, 'colors', 'red,blue,green,black,orange').split(',')]
        self.ancestor_top = None if getattr(self, 'ancestor_top', 'none') == 'none' else self.ancestor_top
        self.ancestor_left = None if getattr(self, 'ancestor_left', 'none') == 'none' else self.ancestor_left

        self.blockinfo_reverse = base.str_to_bool(self.blockinfo_reverse)

    def pair_position(self, alignment, loc1, loc2, colors):
        alignment.index = alignment.index.map(loc1)
        data = []
        for i, k in enumerate(alignment.columns):
            df = alignment[k].map(loc2).dropna()
            for idx, row in df.items():
                data.append([idx, row, colors[i]])
        return pd.DataFrame(data, columns=['loc1', 'loc2', 'color'])

    def run(self):
        axis = [0, 1, 1, 0]

        # Lens generation and figure size
        lens1 = base.newlens(self.lens1, self.position)
        lens2 = base.newlens(self.lens2, self.position)
        
        if re.search(r'\d', self.figsize):
            self.figsize = [float(k) for k in self.figsize.split(',')]
        else:
            self.figsize = np.array([1, float(lens1.sum()) / float(lens2.sum())]) * 10
            
        plt.rcParams['ytick.major.pad'] = 0

        # Create plot
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.xaxis.set_ticks_position('top')
        step1, step2 = 1 / float(lens1.sum()), 1 / float(lens2.sum())

        # Process Ancestor Data
        if self.ancestor_left:
            axis[0] = -0.02
            lens_ancestor_left = self.process_ancestor(self.ancestor_left, lens1.index)

        if self.ancestor_top:
            axis[3] = -0.02
            lens_ancestor_top = self.process_ancestor(self.ancestor_top, lens2.index)

        base.dotplot_frame(fig, ax, lens1, lens2, step1, step2, 
                           self.genome1_name, self.genome2_name, [0, 1])

        # Process GFF files
        gff1, gff2 = base.newgff(self.gff1), base.newgff(self.gff2)
        gff1 = base.gene_location(gff1, lens1, step1, self.position)
        gff2 = base.gene_location(gff2, lens2, step2, self.position)

        if self.ancestor_top:
            self.ancestor_position(ax, gff2, lens_ancestor_top, 'top')

        if self.ancestor_left:
            self.ancestor_position(ax, gff1, lens_ancestor_left, 'left')

        # Process block info and alignment
        bkinfo = self.process_blockinfo(lens1,lens2)
        align = self.alignment(gff1, gff2, bkinfo)
        alignment = align[gff1.columns[-len(bkinfo[self.classid].drop_duplicates()):]]
        alignment.to_csv(self.savefile, header=False)

        # Create scatter plot
        df = self.pair_position(alignment, gff1['loc'], gff2['loc'], self.colors)
        plt.scatter(df['loc2'], df['loc1'], s=float(self.markersize), c=df['color'], 
                    alpha=0.5, edgecolors=None, linewidths=0, marker='o')

        ax.axis(axis)
        plt.subplots_adjust(left=0.07, right=0.97, top=0.93, bottom=0.03)
        plt.savefig(self.savefig, dpi=500)
        plt.show()

    def process_ancestor(self, ancestor_file, lens_index):
        df = pd.read_csv(ancestor_file, sep="\t", header=None)
        df[0] = df[0].astype(str)
        df[3] = df[3].astype(str)
        df[4] = df[4].astype(int)
        df[4] = df[4] / df[4].max()
        return df[df[0].isin(lens_index)]

    def process_blockinfo(self, lens1, lens2):
        bkinfo = pd.read_csv(self.blockinfo, index_col='id')
        if self.blockinfo_reverse ==  True:
            bkinfo[['chr1', 'chr2']] = bkinfo[['chr2', 'chr1']]
            bkinfo[['block1', 'block2']] = bkinfo[['block2', 'block1']]
        bkinfo['chr1'] = bkinfo['chr1'].astype(str)
        bkinfo['chr2'] = bkinfo['chr2'].astype(str)
        bkinfo[self.classid] = bkinfo[self.classid].astype(str)
        return bkinfo[bkinfo['chr1'].isin(lens1.index) & (bkinfo['chr2'].isin(lens2.index))]

    def alignment(self, gff1, gff2, bkinfo):
        gff1['uid'] = gff1['chr'] + 'g' + gff1['order'].astype(str)
        gff2['uid'] = gff2['chr'] + 'g' + gff2['order'].astype(str)
        gff1['id'] = gff1.index
        gff2['id'] = gff2.index
        
        for cl, group in bkinfo.groupby(self.classid):
            name = f'l{cl}'
            gff1[name] = ''
            group = group.sort_values(by=['length'], ascending=True)

            for _, row in group.iterrows():
                block = self.create_block_dataframe(row)
                if block.empty:
                    continue
                block1_min, block1_max = block['block1'].agg(['min', 'max'])
                area = gff1[(gff1['chr'] == row['chr1']) & 
                            (gff1['order'] >= block1_min) & 
                            (gff1['order'] <= block1_max)].index
                
                block['id1'] = (row['chr1'] + 'g' + block['block1'].astype(str)).map(
                    dict(zip(gff1['uid'], gff1.index)))
                block['id2'] = (row['chr2'] + 'g' + block['block2'].astype(str)).map(
                    dict(zip(gff2['uid'], gff2.index)))

                gff1.loc[block['id1'].values, name] = block['id2'].values
                gff1.loc[gff1.index.isin(area) & gff1[name].eq(''), name] = '.'
        return gff1

    def create_block_dataframe(self, row):
        b1, b2, ks = row['block1'].split('_'), row['block2'].split('_'), row['ks'].split('_')
        ks = list(map(float, ks[1:])) if ks[0] == '' else list(map(float, ks))
        block = pd.DataFrame(np.array([b1, b2, ks]).T, columns=['block1', 'block2', 'ks'])
        block['block1'] = block['block1'].astype(int)
        block['block2'] = block['block2'].astype(int)
        block['ks'] = block['ks'].astype(float)
        return block[(block['ks'] <= self.ks_area[1]) & 
                     (block['ks'] >= self.ks_area[0])].drop_duplicates(subset=['block1'], keep='first')

    def ancestor_position(self, ax, gff, lens, mark):
        for _, row in lens.iterrows():
            loc1 = gff[(gff['chr'] == row[0]) & (gff['order'] == int(row[1]))].index
            loc2 = gff[(gff['chr'] == row[0]) & (gff['order'] == int(row[2]))].index
            loc1, loc2 = gff.loc[[loc1[0], loc2[0]], 'loc']
            if mark == 'top':
                width = abs(loc1-loc2)
                loc = [min(loc1, loc2), 0]
                height = -0.02
            if mark == 'left':
                height = abs(loc1-loc2)
                loc = [-0.02, min(loc1, loc2), ]
                width = 0.02
            base.Rectangle(ax, loc, height, width, row[3], row[4])