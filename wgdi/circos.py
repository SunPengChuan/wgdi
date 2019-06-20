import re
import sys

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wgdi.base as base


class circos():
    def __init__(self, options):
        self.figsize = '10,10'
        self.position = 'order'
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.ring_width = float(self.ring_width)

    def plot_circle(self, loc_chr, radius, color='black', lw=1, alpha=1, linestyle='-'):
        for k in loc_chr:
            start, end = loc_chr[k]
            t = np.arange(start, end, 0.005)
            x, y = (radius) * np.cos(t), (radius) * np.sin(t)
            plt.plot(x, y, linestyle=linestyle,
                     color=color, lw=lw, alpha=alpha)

    def plot_labels(self, labels, loc_chr, radius, horizontalalignment="center", verticalalignment="center", fontsize=6,
                    color='black'):
        for k in loc_chr:
            loc = sum(loc_chr[k]) * 0.5
            x, y = radius * np.cos(loc), radius * np.sin(loc)
            if 1 * np.pi < loc < 2 * np.pi:
                loc += np.pi
            plt.text(x, y, labels[k], horizontalalignment=horizontalalignment, verticalalignment=verticalalignment,
                     fontsize=fontsize, color=color, rotation=loc * 180 / np.pi - 90)

    def Wedge(self, ax, loc, radius, start, end, width, color, alpha):
        p = mpatches.Wedge(loc, radius, start, end, width=width,
                           edgecolor=None, facecolor=color, alpha=alpha)
        ax.add_patch(p)

    def plot_bar(self, df, radius, length, lw, color, alpha):
        for k in df[df.columns[0]].drop_duplicates().values:
            if k in ['', np.nan]:
                continue
            df_chr = df.groupby(df.columns[0]).get_group(k)
            x1, y1 = radius * \
                np.cos(df_chr['rad']), radius * np.sin(df_chr['rad'])
            x2, y2 = (radius + length) * \
                np.cos(df_chr['rad']), (radius + length) * \
                np.sin(df_chr['rad'])
            x = np.array(
                [x1.values, x2.values, [np.nan] * x1.size]).flatten('F')
            y = np.array(
                [y1.values, y2.values, [np.nan] * x1.size]).flatten('F')
            plt.plot(x, y, linestyle='-',
                     color=color[str(k)], lw=lw, alpha=alpha)

    def chr_loction(self, lens, angle_gap, angle):
        start, end, loc_chr = 0, 0, {}
        for k in lens.index:
            end += angle_gap + angle * (float(lens[k]))
            start = end - angle * (float(lens[k]))
            loc_chr[k] = [float(start), float(end)]
        return loc_chr

    def deal_alignment(self, alignment, gff, lens, loc_chr, angle):
        alignment.replace('\s+', '', inplace=True)
        alignment.replace('.', '', inplace=True)
        # newalignment = alignment.replace(gff['chr'])
        print(alignment.info())
        # newalignment = alignment.applymap(gff['chr'])
        newalignment = alignment.copy()
        for i in  range(10):
            alignment[i] = alignment[i].astype(str)
            newalignment[i] = alignment[i].map(gff['chr'].to_dict())
            print(i)
        newalignment['loc'] = alignment[0].replace(gff[self.position])
        newalignment[0] = newalignment[0].astype('str')
        newalignment['loc'] = newalignment['loc'].astype('float')
        newalignment = newalignment[newalignment[0].isin(lens.index) == True]
        newalignment['rad'] = np.nan
        for name, group in newalignment.groupby([0]):
            if str(name) not in loc_chr:
                continue
            newalignment.loc[group.index, 'rad'] = loc_chr[str(
                name)][0]+angle * group['loc']
        return newalignment

    def run(self):
        fig = plt.figure(figsize=(tuple(self.figsize)))
        root = plt.axes([0, 0, 1, 1])
        mpl.rcParams['agg.path.chunksize'] = 100000000
        lens = base.newlens(self.lens1, self.position)
        radius, angle_gap = float(self.radius), float(self.angle_gap)
        angle = (2 * np.pi - (int(len(lens))) * angle_gap) / (int(lens.sum()))
        loc_chr = self.chr_loction(lens, angle_gap, angle)
        list_colors = [str(k).strip() for k in re.split(',|:', self.colors)]
        chr_color = dict(zip(list_colors[::2], list_colors[1::2]))
        for k in loc_chr:
            start, end = loc_chr[k]
            self.Wedge(root, (0.0, 0.0), radius + 0.03, start * 180 /
                       np.pi, end * 180 / np.pi, 0.03, chr_color[k], 0.9)
        gff = pd.read_csv(self.gff, sep='\t', header=None, index_col=1)
        gff.rename(columns={0: 'chr', 1: 'id', 2: 'start',
                            3: 'end', 5: 'order'}, inplace=True)
        alignment = pd.read_csv(self.alignment, sep='\t', header=None)
        newalignment = self.deal_alignment(
            alignment, gff, lens, loc_chr, angle)
        for k, v in enumerate(newalignment.columns[1:-2]):
            r = radius + self.ring_width*(k+1)
            self.plot_circle(loc_chr, r, lw=0.5, alpha=0.5, color='grey')
            self.plot_bar(newalignment[[v, 'rad']], r + self.ring_width *
                          0.15, self.ring_width*0.7, 0.15, chr_color, 0.9)
        labels = self.chr_label + lens.index
        labels = dict(zip(lens.index, labels))
        self.plot_labels(labels, loc_chr, radius - 0.03, fontsize=9)
        root.set_xlim(-1, 1)
        root.set_ylim(-1.05, 0.95)
        root.set_axis_off()
        plt.savefig(self.savefile, dpi=500)
        sys.exit(0)
