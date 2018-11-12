import re
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class align_dotplot():
    def __init__(self, options):
        for k, v in options:
            setattr(self, str(k), v)
            print(str(k), ' = ', v)
        self.colors = [str(k) for k in self.colors.split(',')]
        self.figsize = [float(k) for k in self.figsize.split(',')]

    def gene_location(self, gff, lens, step):
        loc_gene, dict_chr, n = {}, {}, 0
        for i in lens.index:
            dict_chr[str(i)] = n
            n += float(lens.at[i, 1])
        for k in gff.index:
            if gff.loc[k, 0] not in dict_chr:
                continue
            loc = (float(dict_chr[gff.loc[k, 0]])+float(gff.loc[k, 5]))*step
            loc_gene[gff.loc[k, 1]] = loc
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
            n += float(lens.at[k, 1])
            mark_new = str(mark)+str(k)
            x = gl_start-float(n)*step
            mark_x = x+0.5*float(lens.at[k, 1])*step
            plt.plot([start_x, start_x+gl2], [x, x],
                     linestyle='-', color='black', linewidth=0.5)
            plt.text(mark_y, mark_x, mark_new, color='black',
                     fontsize=12, rotation=90, weight='semibold', **align)
        plt.plot([start_x, start_x+gl2], [gl_start, gl_start],
                 linestyle='-', color='black', linewidth=1)
        plt.text(mark_y-0.02, 0.5*(2*gl_start-gl), name, color='black',
                 fontsize=18, rotation=90, weight='semibold', **align)

    def plot_chr2(self, lens, gl, gl2, step, mark, name):
        gl_start, n, start_x = 0.05, 0, 0.95
        mark_y = 0.96
        align = dict(family='Times New Roman', style='normal',
                     horizontalalignment="center", verticalalignment="center")
        for k in lens.index:
            n += float(lens.at[k, 1])
            mark_new = str(mark)+str(k)
            x = gl_start+float(n)*step
            mark_x = x-0.5*float(lens.at[k, 1])*step
            plt.plot([x, x], [start_x, start_x-gl2],
                     linestyle='-', color='black', linewidth=0.5)
            plt.text(mark_x, mark_y, mark_new, color='black',
                     fontsize=12, rotation=0, weight='semibold', **align)
        plt.plot([gl_start, gl_start], [start_x, start_x-gl2],
                 linestyle='-', color='black', linewidth=1)
        plt.text(0.5*(2*gl_start+gl), mark_y+0.02, name, color='black',
                 fontsize=18, rotation=0, weight='semibold', **align)

    def run(self):
        fig = plt.figure(figsize=tuple(self.figsize))
        plt.axis('off')
        gff_1 = pd.read_csv(self.gff1, sep="\t", header=None)
        gff_2 = pd.read_csv(self.gff2, sep="\t", header=None)
        gff_1[0] = gff_1[0].astype('str')
        gff_2[0] = gff_2[0].astype('str')
        lens_1 = pd.read_csv(self.lens1, sep="\t", header=None, index_col=0)
        lens_2 = pd.read_csv(self.lens2, sep="\t", header=None, index_col=0)
        gl1, gl2 = 0.92, 0.92
        step1 = gl1 / float(lens_1[1].sum())
        step2 = gl2 / float(lens_2[1].sum())
        self.plot_chr1(lens_1, gl1, gl2, step1, '', self.genome1_name)
        self.plot_chr2(lens_2, gl1, gl2, step2, '', self.genome2_name)
        gene_loc_1 = self.gene_location(gff_1, lens_1, step1)
        gene_loc_2 = self.gene_location(gff_2, lens_2, step2)
        alignment = pd.read_csv(self.alignment, sep='\t',
                                header=None, index_col=0)
        alignment.replace('\s+', '', inplace=True)
        gene_loc_1 = self.gene_location(gff_1, lens_1, step1)
        gene_loc_2 = self.gene_location(gff_2, lens_2, step2)
        for k in alignment.columns:
            y, x = self.pair_positon(alignment[k], gene_loc_1, gene_loc_2)
            cols = [self.colors[k-1]]*len(x)
            plt.scatter(x, y, s=float(self.markersize), c=cols, alpha=0.5,
                        edgecolors=None, linewidths=0, marker=(9, 3, 30))
        plt.subplots_adjust(left=0.02, right=1, top=0.98, bottom=0)
        plt.savefig(self.savefile, dpi=500)
        sys.exit(0)
