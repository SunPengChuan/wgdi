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
        self.label_size = 9
        self.label_radius = 0.015
        self.column_names = [None]*100
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.figsize = [float(k) for k in self.figsize.split(',')]
        self.ring_width = float(self.ring_width)
        if hasattr(self, 'legend_square'):
            self.legend_square = [float(k)
                                  for k in self.legend_square.split(',')]
        else:
            self.legend_square = 0.04, 0.04

    def plot_circle(self, loc_chr, radius, color='black', lw=1, alpha=1, linestyle='-'):
        for k in loc_chr:
            start, end = loc_chr[k]
            t = np.arange(start, end, 0.005)
            x, y = (radius) * np.cos(t), (radius) * np.sin(t)
            plt.plot(x, y, linestyle=linestyle,
                     color=color, lw=lw, alpha=alpha)

    def plot_labels(self, root, labels, loc_chr, radius, horizontalalignment="center", verticalalignment="center", fontsize=6,
                    color='black'):
        for k in loc_chr:
            loc = sum(loc_chr[k]) * 0.5
            x, y = radius * np.cos(loc), radius * np.sin(loc)
            self.Wedge(root, (x, y), self.label_radius, 0,
                       360, self.label_radius, 'white', 1)
            if 1 * np.pi < loc < 2 * np.pi:
                loc += np.pi
            plt.text(x, y, labels[k], horizontalalignment=horizontalalignment, verticalalignment=verticalalignment,
                     fontsize=fontsize, color=color, rotation=0)

    def Wedge(self, ax, loc, radius, start, end, width, color, alpha):
        p = mpatches.Wedge(loc, radius, start, end, width=width,
                           edgecolor=None, facecolor=color, alpha=alpha)
        ax.add_patch(p)

    def plot_bar(self, df, radius, length, lw, color, alpha):
        for k in df[df.columns[0]].drop_duplicates().values:
            if str(k) not in color.keys():
                color[str(k)] = 'black'
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
        start, end, loc_chr = 0, 0.2*angle_gap, {}
        for k in lens.index:
            end += angle_gap + angle * (float(lens[k]))
            start = end - angle * (float(lens[k]))
            loc_chr[k] = [float(start), float(end)]
        return loc_chr

    def deal_alignment(self, alignment, gff, lens, loc_chr, angle):
        alignment.replace('\s+', '', inplace=True)
        alignment.replace('.', '', inplace=True)
        newalignment = alignment.copy()
        for i in range(len(alignment.columns)):
            alignment[i] = alignment[i].astype(str)
            newalignment[i] = alignment[i].map(gff['chr'].to_dict())
        newalignment['loc'] = alignment[0].map(gff[self.position].to_dict())
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

    def deal_ancestor(self, alignment, gff, lens, loc_chr, angle, al):
        alignment.replace('\s+', '', inplace=True)
        alignment.replace('.', np.nan, inplace=True)
        newalignment = pd.merge(alignment, gff, left_on=0, right_on=gff.index)
        newalignment['rad'] = np.nan
        for name, group in newalignment.groupby(['chr']):
            if str(name) not in loc_chr:
                continue
            newalignment.loc[group.index, 'rad'] = loc_chr[str(
                name)][0]+angle * group[self.position]
        newalignment.index = newalignment[0]
        newalignment[0] = newalignment[0].map(newalignment['rad'].to_dict())
        data = []
        for index_al, row_al in al.iterrows():
            for k in alignment.columns[1:]:
                alignment[k] = alignment[k].astype(str)
                group = newalignment[(newalignment['chr'] == row_al['chr']) & (
                    newalignment['order'] >= row_al['start']) & (newalignment['order'] <= row_al['end'])].copy()
                group.loc[:, k] = group.loc[:, k].map(
                    newalignment['rad']).values
                group.dropna(subset=[k], inplace=True)
                group.index = group.index.map(newalignment['rad'].to_dict())
                group['color'] = row_al['color']
                group = group[group[k].notnull()]
                data += group[[0, k, 'color']].values.tolist()
        df = pd.DataFrame(data, columns=['loc1', 'loc2', 'color'])
        return df

    def plot_collinearity(self, data, radius, lw=0.02, alpha=1):
        for name, group in data.groupby('color'):
            x, y = np.array([]), np.array([])
            for index, row in group.iterrows():
                ex1x, ex1y = radius * \
                    np.cos(row['loc1']), radius*np.sin(row['loc1'])
                ex2x, ex2y = radius * \
                    np.cos(row['loc2']), radius*np.sin(row['loc2'])
                ex3x, ex3y = radius * (1-abs(row['loc1']-row['loc2'])/np.pi) * np.cos((row['loc1']+row['loc2'])*0.5), radius * (
                    1-abs(row['loc1']-row['loc2'])/np.pi) * np.sin((row['loc1']+row['loc2'])*0.5)
                x1 = [ex1x, 0.5*ex3x, ex2x]
                y1 = [ex1y, 0.5*ex3y, ex2y]
                step = .002
                t = np.arange(0, 1+step, step)
                xt = base.Bezier3(x1, t)
                yt = base.Bezier3(y1, t)
                x = np.hstack((x, xt, np.nan))
                y = np.hstack((y, yt, np.nan))
            plt.plot(x, y, color=name, lw=lw, alpha=alpha)

    def plot_legend(self, ax, chr_color, width, height):
        (x1, x2) = ax.get_xlim()
        (y1, y2) = ax.get_ylim()
        a = 1000
        for k, v in enumerate(chr_color.keys(), 0):
            h = y1-k//a*height*2
            k = k % a
            if x1 + width * k > x2-width:
                a = k
                h = y1-k//a*height*2
                k = k % a
            loc = [x1 + width * k, h]
            base.Rectangle(ax, loc, height, width, chr_color[v], 1)
            plt.text(loc[0] + width*0.382, h-0.618*height, v, fontsize=12)
        ax.set_ylim(h-2*height, y2)

    def run(self):
        fig, ax = plt.subplots(figsize=self.figsize)
        mpl.rcParams['agg.path.chunksize'] = 100000000
        lens = base.newlens(self.lens, self.position)
        radius, angle_gap = float(self.radius), float(self.angle_gap)
        angle = (2 * np.pi - (int(len(lens))+1.5)
                 * angle_gap) / (int(lens.sum()))
        loc_chr = self.chr_loction(lens, angle_gap, angle)
        list_colors = [str(k).strip() for k in re.split(',|:', self.colors)]
        chr_color = dict(zip(list_colors[::2], list_colors[1::2]))
        gff = base.newgff(self.gff)
        if hasattr(self, 'ancestor'):
            ancestor = pd.read_csv(self.ancestor, header=None)
            al = pd.read_csv(self.ancestor_location, sep='\t', header=None)
            al.rename(columns={0: 'chr', 1: 'start',
                               2: 'end', 3: 'color'}, inplace=True)
            al['chr'] = al['chr'].astype(str)
            data = self.deal_ancestor(ancestor, gff, lens, loc_chr, angle, al)
            self.plot_collinearity(data, radius, lw=0.1, alpha=0.8)

        if hasattr(self, 'alignment'):
            alignment = pd.read_csv(self.alignment, header=None)
            newalignment = self.deal_alignment(
                alignment, gff, lens, loc_chr, angle)
            if ',' in self.column_names:
                names = [str(k) for k in self.column_names.split(',')]
            else:
                names = [None]*len(newalignment.columns)
            n = 0
            align = dict(family='Arial', verticalalignment="center",
                         horizontalalignment="center")
            for k, v in enumerate(newalignment.columns[1:-2]):
                r = radius + self.ring_width*(k+1)
                self.plot_circle(loc_chr, r, lw=0.5, alpha=1, color='grey')
                self.plot_bar(newalignment[[v, 'rad']], r + self.ring_width *
                              0.15, self.ring_width*0.7, 0.15, chr_color, 1)
                if n % 2 == 0:
                    loc = 0.05
                    x, y = (r+self.ring_width*0.5) * \
                        np.cos(loc), (r+self.ring_width*0.5) * np.sin(loc)
                    plt.text(x, y, names[n], rotation=loc *
                             180 / np.pi, fontsize=self.label_size, **align)
                else:
                    loc = -0.08
                    x, y = (r+self.ring_width*0.5) * \
                        np.cos(loc), (r+self.ring_width*0.5) * np.sin(loc)
                    plt.text(x, y, names[n], fontsize=self.label_size,
                             rotation=loc * 180 / np.pi, **align)
                n += 1
        if hasattr(self, 'ancestor'):
            colors = al['color'].drop_duplicates().values.tolist()
            ancestor_chr_color = dict(zip(range(1, len(colors)+1), colors))
            self.plot_legend(ax, ancestor_chr_color,
                             self.legend_square[0], self.legend_square[1])
        if hasattr(self, 'alignment'):
            del chr_color['nan']
            self.plot_legend(
                ax, chr_color, self.legend_square[0], self.legend_square[1])
        labels = self.chr_label + lens.index
        labels = dict(zip(lens.index, labels))
        self.plot_labels(ax, labels, loc_chr, radius +
                         self.ring_width*0.3, fontsize=self.label_size)

        plt.axis('off')
        a = (ax.get_ylim()[1]-ax.get_ylim()[0]) / \
            (ax.get_xlim()[1]-ax.get_xlim()[0])
        fig.set_size_inches(self.figsize[0], self.figsize[0]*a, forward=True)
        plt.savefig(self.savefig, dpi=500)
        plt.show()
        sys.exit(0)
