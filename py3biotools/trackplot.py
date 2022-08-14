import pyBigWig
import pyranges as pr
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from .utils import cluster_peaks


class FritoBaseTrack(object):
    """Track for Frito output"""
    def __init__(self, df, name, position="Pos", y='Score', color="grey", weight=1, ymin=None, ymax=None):
        self.data = df
        self.name = name
        self.position = position
        self.y = y
        self.color = color
        self.weight = weight
        self.ymin = ymin
        self.ymax = ymax

    def plot(self, ax, st, en):
        x = np.arange(st, en)
        dat = self.data.groupby(self.position).apply(lambda x: x.sort_values(by=self.y).iloc[0])
        y = dat.set_index(self.position).reindex(x)[self.y].fillna(np.inf).values
        sns.lineplot(ax=ax, x=x, y=y, color=self.color)
        y_valid = [i for i in y if np.abs(i) != np.inf]
        return np.min(y_valid), np.max(y_valid)
    

class FritoPeakTrack(object):
    """Track for Frito output"""
    def __init__(self, df, name, position="Pos", score_cols=['Neg.Pval', 'Neg.FC'], ci=2, color="grey", weight=1, ymin=None, ymax=None):
        x = df.loc[(df[score_cols[0]] >= ci) & (df[score_cols[1]] >= ci)][position].values
        self.peaks = cluster_peaks(x)
        self.name = name
        self.color = color
        self.weight = weight
        self.ymin = ymin
        self.ymax = ymax
    
    def plot(self, ax, height=1):
        for st, en in self.peaks:
            ax.add_patch(mpatches.Rectangle((st, 0), en-st, height, color=self.color))
        return 0, height

    
class RangeTrack(object):
    """Track for pyRanges"""
    def __init__(self, ranges, name, y="FC", hue='RBP', palette=None, weight=1, ymin=None, ymax=None):
        self.data = ranges
        self.name = name
        self.y = y
        self.hue= hue
        self.palette = palette
        self.weight = weight
        self.ymin = ymin
        self.ymax = ymax

    def plot(self, ax, chrom=None, strand=None, start=None, end=None):
        # Select data in the region
        if chrom is None:
            raise KeyError("No chrom specified!")
        else:
            if start is None or end is None or strand is None:
                self._pdata = self.data[chrom].df
            else:
                self._pdata = self.data[chrom, strand, start:end].df
        
        # Create palette
        if self.palette is None:
            self._labels =  sorted(np.unique(self._pdata[self.hue]))
            self._n_colors = len(self._labels)
            self._palette = dict(zip(self._labels, sns.color_palette("Spectral", self._n_colors)
            ))
        else:
            self._palette = self.palette
            self._labels = list(self._palette)
            self._n_colors = len(self._labels)
        
        # Plot
        self._pdata.apply(lambda x: self.add_patch(ax, x['Start'], x['End'], x[self.y], self._palette[x[self.hue]], .75), axis=1)
        return self._pdata[self.y].min(), self._pdata[self.y].max()
        
    @staticmethod
    def add_patch(ax, st, en, score, color, alpha):
        ax.add_patch(mpatches.Rectangle((st, 0), en-st, score, color=color, alpha=alpha))
        return 0
        

class BigwigTrack(object):
    """Track for bigWig"""
    def __init__(self, infile, name, color="grey", weight=1, ymin=None, ymax=None):
        self.infile = infile
        self.color = color
        self.name = name
        self.weight = weight
        self.ymin = ymin
        self.ymax = ymax
    
    def plot(self, ax, chrom, start, end):
        bw = pyBigWig.open(self.infile)
        self.cov = bw.stats(chrom, start, end, nBins=end-start)
        self.x = np.arange(start, end)
        self.y = [np.inf if i == 0 else i for i in self.cov]
        sns.lineplot(ax=ax, x=self.x, y=self.y, color=self.color, alpha=.5)
        ax.fill_between(x=self.x, y1=self.y, color=self.color, alpha=.75)
        bw.close()
        return np.min(self.cov), np.max(self.cov)
        
        
class GtfTrack(object):
    "Gene Structure"
    def __init__(self, gencode_gtf, name, gene_id, weight=2, ymin=None, ymax=None):
        self.gtf = gencode_gtf
        self.name = name
        self.weight = weight
        self.gene_id = gene_id
        self.ymin = ymin
        self.ymax = ymax

    def plot(self, ax, start, end):
        gene_tscps = self.gtf[(self.gtf.Feature == 'transcript') & (self.gtf.gene_id == self.gene_id)].transcript_name.values
        for idx, _tscp_id in enumerate(gene_tscps):
            last_x = None
            _exons = self.gtf[(self.gtf.Feature == 'exon') & (self.gtf.transcript_name == _tscp_id)]
            for st, en in _exons.df.apply(lambda x: [x['Start'], x['End']], axis=1).values.tolist():
                ax.add_patch(mpatches.Rectangle((st, -1-idx), en-st, .6, color="#0000b2"))
                if last_x is not None:
                    ax.plot([last_x, st], [-1-idx+.25, -1-idx+.25], color="#0000b2", linewidth=1)
                last_x = en
        return -1-idx-1, 0

    
class TrackPlot(object):
    def __init__(self, tracks, height_ratios=None):
        self.tracks = tracks
        if height_ratios is None:
            height_ratios = [i.weight for i in self.tracks]
        self.height_ratios = height_ratios

        
    def plot(self, figsize, chrom, strand, st, en, **kwargs):        
        self.fig, self.axes = plt.subplots(figsize=figsize, nrows=len(self.tracks), gridspec_kw={'height_ratios': self.height_ratios}, facecolor="white", sharex=False, **kwargs)
        for idx, track in enumerate(self.tracks):
            ax = self.axes[idx]
            if isinstance(track, RangeTrack):
                ymin, ymax = track.plot(ax, chrom, strand, st, en)
                ymin = np.floor(ymin) if track.ymin is None else track.ymin
                ymax = np.ceil(ymax) if track.ymax is None else track.ymax
                ax.set(xlim=[st, en], xticks=[], ylim=[0, ymax], yticks=[0, ymax])
                if idx == 0:
                    ax.legend(
                        handles=[mpatches.Patch(facecolor=track._palette[x], edgecolor=track._palette[x], label=x) for x in track._labels], 
                        ncol=9, loc="lower left", bbox_to_anchor=(0, 1.15), borderaxespad=0., frameon=False, fontsize=14
                    )
                sns.despine(ax=ax, bottom=True, top=True, left=False, right=True, trim=True)
            elif isinstance(track, FritoBaseTrack):
                ymin, ymax = track.plot(ax, st, en)
                ymin = ymin if track.ymin is None else track.ymin
                ymax = ymax if track.ymax is None else track.ymax
                ax.set(xlim=[st, en], xticks=[], ylim=[ymin, ymax], yticks=[ymin, ymax])
                sns.despine(ax=ax, bottom=True, top=True, left=False, right=True, trim=True)
            elif isinstance(track, FritoPeakTrack):
                ymin, ymax = track.plot(ax, height=1)
                ymin = ymin if track.ymin is None else track.ymin
                ymax = ymax if track.ymax is None else track.ymax
                ax.set(xlim=[st, en], xticks=[], ylim=[ymin, ymax], yticks=[])
                sns.despine(ax=ax, bottom=True, top=True, left=True, right=True, trim=True)
            elif isinstance(track, BigwigTrack):
                ymin, ymax = track.plot(ax, chrom, st, en)
                ymin = ymin if track.ymin is None else track.ymin
                ymax = ymax if track.ymax is None else track.ymax
                ax.set(xlim=[st, en], xticks=[], ylim=[ymin, ymax], yticks=[ymin, ymax])
                sns.despine(ax=ax, bottom=True, top=True, left=False, right=True, trim=True)  
            elif isinstance(track, GtfTrack):
                ymin, ymax = track.plot(ax, st, en)
                ymin = ymin if track.ymin is None else track.ymin
                ymax = ymax if track.ymax is None else track.ymax
                ax.set(xlim=[st, en], ylim=[ymin, ymax], yticks=[])
                sns.despine(ax=ax, bottom=False, top=True, left=True, right=True, trim=False) 
            else:
                pass        
        
            ax.yaxis.set_label_position("left")
            ax.set_ylabel(track.name, labelpad=15, fontsize=14, rotation=0, va='center', ha='right')
