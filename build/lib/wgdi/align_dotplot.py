import re
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class align_dotplot():
    def __init__(self, options):
        self.position = 'order'
        self.figsize = 'default'
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.colors = [str(k) for k in self.colors.split(',')]

    def gene_location(self, gff, lens, step):
        loc_gene, dict_chr, n = {}, {}, 0
        for i in lens.index:
            dict_chr[str(i)] = n
            n += float(lens[i])
        for k in gff.index:
            if gff.loc[k, 'chr'] not in dict_chr:
                continue
            loc = (float(dict_chr[gff.loc[k, 'chr']]) +
                   float(gff.loc[k, self.position])) * step
            loc_gene[gff.loc[k, 'id']] = loc
        return loc_gene

    def pair_positon(self, alignment, loc1, loc2):
        pos1, pos2 = [], []
        gl_start1, gl_start2 = 0.95, 0.05
        alignment = alignment[alignment.str.match(r'\w+') == True]
        for index in alignment.index:
            if (index not in loc1) or (alignment[index] not in loc2):
                continue
            pos1.append(gl_start1 - loc1[index])
            pos2.append(gl_start2 + loc2[alignment[index]])
        return pos1, pos2

    def plot_chr1(self, lens, gl, gl2, step, mark, name):
        gl_start, n, start_x = 0.95, 0, 0.05
        mark_y = 0.04
        align = dict(family='Times New Roman', style='normal',
                     horizontalalignment="center", verticalalignment="center")
        for k in lens.index:
            n += float(lens[k])
            mark_new = str(mark) + str(k)
            x = gl_start - float(n) * step
            mark_x = x + 0.5 * float(lens[k]) * step
            plt.plot([start_x, start_x + gl2], [x, x],
                     linestyle='-', color='black', linewidth=0.5)
            plt.text(mark_y, mark_x, mark_new, color='black',
                     fontsize=12, rotation=90, weight='semibold', **align)
        plt.plot([start_x, start_x + gl2], [gl_start, gl_start],
                 linestyle='-', color='black', linewidth=1)
        plt.text(mark_y - 0.02, 0.5 * (2 * gl_start - gl), name, color='black', fontsize=18, rotation=90,
                 weight='semibold', **align)

    def plot_chr2(self, lens, gl, gl2, step, mark, name):
        gl_start, n, start_x = 0.05, 0, 0.95
        mark_y = 0.96
        align = dict(family='Times New Roman', style='normal',
                     horizontalalignment="center", verticalalignment="center")
        for k in lens.index:
            n += float(lens[k])
            mark_new = str(mark) + str(k)
            x = gl_start + float(n) * step
            mark_x = x - 0.5 * float(lens[k]) * step
            plt.plot([x, x], [start_x, start_x - gl2],
                     linestyle='-', color='black', linewidth=0.5)
            plt.text(mark_x, mark_y, mark_new, color='black',
                     fontsize=12, rotation=0, weight='semibold', **align)
        plt.plot([gl_start, gl_start], [start_x, start_x - gl2],
                 linestyle='-', color='black', linewidth=1)
        plt.text(0.5 * (2 * gl_start + gl), mark_y + 0.02, name, color='black', fontsize=18, rotation=0,
                 weight='semibold', **align)

    def run(self):
        gff_1 = pd.read_csv(self.gff1, sep="\t", header=None)
        gff_2 = pd.read_csv(self.gff2, sep="\t", header=None)
        gff_1.rename(columns={0: 'chr', 1: 'id', 2: 'start',
                              3: 'end', 5: 'order'}, inplace=True)
        gff_2.rename(columns={0: 'chr', 1: 'id', 2: 'start',
                              3: 'end', 5: 'order'}, inplace=True)
        gff_1['chr'] = gff_1['chr'].astype('str')
        gff_2['chr'] = gff_2['chr'].astype('str')
        lens_1 = pd.read_csv(self.lens1, sep="\t", header=None, index_col=0)
        lens_2 = pd.read_csv(self.lens2, sep="\t", header=None, index_col=0)
        gl1, gl2 = 0.92, 0.92
        if self.position == 'order':
            lens_1 = lens_1[2]
            lens_2 = lens_2[2]

        else:
            lens_1 = lens_1[1]
            lens_2 = lens_2[1]
        if re.search('\d', self.figsize):
            self.figsize = [float(k) for k in self.figsize.split(',')]
        else:
            self.figsize = np.array(
                [1, float(lens_1.sum())/float(lens_2.sum())])*10
        step1 = gl1 / float(lens_1.sum())
        step2 = gl2 / float(lens_2.sum())
        fig = plt.figure(figsize=tuple(self.figsize))
        plt.axis('off')
        self.plot_chr1(lens_1, gl1, gl2, step1, '', self.genome1_name)
        self.plot_chr2(lens_2, gl1, gl2, step2, '', self.genome2_name)
        alignment = pd.read_csv(self.alignment, sep='\t',
                                header=None, index_col=0)
        alignment.replace('\s+', '', inplace=True)
        gene_loc_1 = self.gene_location(gff_1, lens_1, step1)
        gene_loc_2 = self.gene_location(gff_2, lens_2, step2)
        for k in alignment.columns:
            y, x = self.pair_positon(alignment[k], gene_loc_1, gene_loc_2)
            cols = [self.colors[k-1]]*len(x)
            plt.scatter(x, y, s=float(self.markersize), c=cols, alpha=0.5,
                        edgecolors=None, linewidths=0, marker='o')
        plt.subplots_adjust(left=0.02, right=1, top=0.98, bottom=0)
        plt.savefig(self.savefile, dpi=500)
        sys.exit(0)
