import matplotlib.pyplot as plt
import pandas as pd
import wgdi.base as base

class retain:
    def __init__(self, options):
        self.position = 'order'
        
        # Initialize the options by setting attributes dynamically
        for k, v in options:
            setattr(self, str(k), v)
            print(f"{str(k)} = {v}")

        # Handle the ylim parameter, which defines the y-axis limits
        self.ylim = [float(k) for k in self.ylim.split(',')] if hasattr(self, 'ylim') else [0, 1]
        
        # Handle the colors and figsize parameters
        self.colors = [str(k) for k in self.colors.split(',')]
        self.figsize = [float(k) for k in self.figsize.split(',')]

    def run(self):
        # Load GFF and lens data
        gff = base.newgff(self.gff)
        lens = base.newlens(self.lens, self.position)
        
        # Filter GFF data based on lens chromosome index
        gff = gff[gff['chr'].isin(lens.index)]
        
        # Load alignment data and join with GFF
        alignment = pd.read_csv(self.alignment, header=None, index_col=0)
        alignment = alignment.join(gff[['chr', self.position]], how='left')
        
        # Perform alignment processing
        self.retain = self.align_chr(alignment)
        
        # Save the processed data to a file
        self.retain[self.retain.columns[:-2]].to_csv(self.savefile, sep='\t', header=None)
        
        # Create a figure for plotting
        fig, axs = plt.subplots(len(lens), 1, sharex=True, sharey=True, figsize=tuple(self.figsize))
        fig.add_subplot(111, frameon=False)
        
        align = dict(family='DejaVu Sans', verticalalignment="center", horizontalalignment="center")

        
        # Hide all the spines and ticks on the plot
        for spine in plt.gca().spines.values():
            spine.set_visible(False)
        plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
        
        # Group the retain data by chromosome and plot each chromosome's data
        groups = self.retain.groupby('chr')
        for i, chr_name in enumerate(lens.index):
            group = groups.get_group(chr_name)

            if len(lens) == 1:
                for j, col in enumerate(self.retain.columns[:-2]):
                    axs.plot(group['order'].values, group[col].values,
                                linestyle='-', color=self.colors[j], linewidth=1)
                axs.spines['right'].set_visible(False)
                axs.spines['top'].set_visible(False)
                axs.set_ylim(self.ylim)
                axs.tick_params(labelsize=12)                
            else:
                # Plot each column's data for the current chromosome
                for j, col in enumerate(self.retain.columns[:-2]):
                    axs[i].plot(group['order'].values, group[col].values,
                                linestyle='-', color=self.colors[j], linewidth=1)
            
                # Hide the right and top spines for each subplot
                axs[i].spines['right'].set_visible(False)
                axs[i].spines['top'].set_visible(False)
                axs[i].set_ylim(self.ylim)
                axs[i].tick_params(labelsize=12)

        for i, chr_name in enumerate(lens.index):
            if len(lens) == 1:
                x, y = axs.get_xlim()[1] * 0.90, axs.get_ylim()[1] * 0.8
                axs.text(x, y, f"{self.refgenome} {chr_name}", fontsize=14, **align)
            else:
                # Add a label for the reference genome and chromosome
                x, y = axs[i].get_xlim()[1] * 0.90, axs[i].get_ylim()[1] * 0.8
                axs[i].text(x, y, f"{self.refgenome} {chr_name}", fontsize=14, **align)
        
        # Adjust layout and save the figure as an image
        plt.ylabel(f"{self.ylabel}\n\n\n\n", fontsize=18, **align)
        plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05)
        plt.savefig(self.savefig, dpi=500)
        plt.show()

    def align_chr(self, alignment):
        """
        Perform the alignment processing for each chromosome by updating the values.
        """
        for i in alignment.columns[:-2]:
            # Update values: set '1' for valid values, '0' for invalid, and fill NaN with 0
            alignment.loc[alignment[i].str.contains(r'\w', na=False), i] = 1
            alignment.loc[alignment[i] == '.', i] = 0
            alignment.loc[alignment[i] == ' ', i] = 0
            alignment[i] = alignment[i].astype('float64').fillna(0)
            
            # Apply the moving average function to each group by chromosome
            for chr_name, group in alignment.groupby(['chr']):
                a = self.moving_average(group[i].values.tolist())
                alignment.loc[group.index, i] = a
        return alignment

    def moving_average(self, arr):
        """
        Calculate a moving average over a specified window size.
        This function smooths the input array using a sliding window.
        """
        a = []
        for i in range(len(arr)):
            # Define the window range
            start, end = max(0, i - int(self.step)), min(len(arr), i + int(self.step))
            ave = sum(arr[start:end]) / (end - start)
            a.append(ave)
        return a
