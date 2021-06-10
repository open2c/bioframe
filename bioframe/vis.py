import itertools
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

from .core import arrops

DEFAULT_FACECOLOR = "skyblue"
DEFAULT_EDGECOLOR = "dimgray"

__all__ = ["plot_intervals"]


def _plot_interval(
    start, end, level, facecolor=None, edgecolor=None, height=0.6, ax=None
):
    facecolor = DEFAULT_FACECOLOR if facecolor is None else facecolor
    edgecolor = DEFAULT_EDGECOLOR if edgecolor is None else edgecolor

    ax = plt.gca() if ax is None else ax
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (start, level - height / 2),
            end - start,
            height,
            facecolor=facecolor,
            edgecolor=edgecolor,
        )
    )


def plot_intervals_arr(
    starts,
    ends,
    levels=None,
    labels=None,
    colors=None,
    xlim=None,
    show_coords=False,
    figsize=(10, 2),
):
    """
    Plot a collection of intervals.

    Parameters
    ----------
    starts, ends : np.ndarray
        A collection of intervals.

    levels : iterable or None
        The level of each interval, i.e. the y-coordinate at which the interval
        must be plotted. If None, it will be determined automatically.

    labels : str or iterable or None
        The label of each interval.

    colors : str or iterable or None.
        The color of each interval.

    xlim : (float, float) or None
        The x-span of the plot.

    show_coords : bool
        If True, plot x-ticks.

    figsize : (float, float) or None.
        The size of the figure. If None, plot within the current figure.

    """
    starts = np.asarray(starts)
    ends = np.asarray(ends)

    if figsize is not None:
        plt.figure(figsize=figsize)

    if levels is None:
        levels = arrops.stack_intervals(starts, ends)
    else:
        levels = np.asarray(levels)

    if isinstance(colors, str) or (colors is None):
        colors = itertools.cycle([colors])
    else:
        colors = itertools.cycle(colors)

    if isinstance(labels, str) or (labels is None):
        labels = itertools.cycle([labels])
    else:
        labels = itertools.cycle(labels)

    for i, (start, end, level, color, label) in enumerate(
        zip(starts, ends, levels, colors, labels)
    ):

        _plot_interval(start, end, level, facecolor=color)
        if label is not None:
            plt.text(
                (start + end) / 2,
                level,
                label,
                horizontalalignment="center",
                verticalalignment="center",
            )

    plt.ylim(-0.5, np.max(levels) + 0.5)
    if xlim is None:
        plt.xlim(-0.5, np.max(ends) + 0.5)
    else:
        plt.xlim(xlim[0], xlim[1])
    plt.gca().set_aspect(1)

    plt.gca().set_frame_on(False)
    plt.yticks([])
    if show_coords:
        pass
    else:
        plt.xticks([])


def plot_intervals(
    df,
    levels=None,
    labels=None,
    colors=None,
    xlim=None,
    show_coords=False,
    figsize=(10, 2),
):
    """
    Plot a collection of intervals, one plot per chromosome.

    Parameters
    ----------
    df : pandas.DataFrame
        A collection of intervals.

    levels : iterable or None
        The level of each interval, i.e. the y-coordinate at which the interval
        must be plotted. If None, it will be determined automatically.

    labels : str or iterable or None
        The label of each interval.

    colors : str or iterable or None.
        The color of each interval.

    xlim : (float, float) or None
        The x-span of the plot.

    show_coords : bool
        If True, plot x-ticks.

    figsize : (float, float) or None.
        The size of the figure. If None, plot within the current figure.

    """
    chrom_gb = df.groupby("chrom")
    for chrom, chrom_df in chrom_gb:
        if isinstance(levels, (list, pd.core.series.Series, np.ndarray)):
            chrom_levels = np.asarray(levels)[chrom_gb.groups[chrom].values]
        elif levels is None:
            chrom_levels = None
        else:
            raise ValueError(f"Unknown type of levels: {type(levels)}")

        if isinstance(labels, (list, pd.core.series.Series, np.ndarray)):
            chrom_labels = np.asarray(labels)[chrom_gb.groups[chrom].values]
        elif labels is None:
            chrom_labels = None
        else:
            raise ValueError(f"Unknown type of labels: {type(levels)}")

        if isinstance(colors, (list, pd.core.series.Series, np.ndarray)):
            chrom_colors = np.asarray(colors)[chrom_gb.groups[chrom].values]
        elif colors is None or isinstance(colors, str):
            chrom_colors = colors
        else:
            raise ValueError(f"Unknown type of colors: {type(colors)}")

        plot_intervals_arr(
            chrom_df.start,
            chrom_df.end,
            levels=chrom_levels,
            labels=chrom_labels,
            colors=chrom_colors,
            xlim=xlim,
            show_coords=show_coords,
            figsize=figsize,
        )
        plt.title(chrom)
