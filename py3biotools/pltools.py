import os
import numpy as np
import pandas
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.path import Path
from matplotlib.patches import PathPatch


def change_width(ax, new_value):
    """
    Change width of seaborn barplot
    """
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)


def setAlpha(ax, a):
    """
    Change alpha of seaborn violinplot
    """
    from matplotlib.collections import PolyCollection
    for art in ax.get_children():
        if isinstance(art, PolyCollection):
            art.set_alpha(a)


def customized_cmap(colors, intervals, N=256):
    """Generate customized cmap
    Generate continuous cmap from given color intervals and length
    :param colors: List of Matplotlib color specifications, or an equivalent Nx3 or Nx4 floating point array (N rgb or rgba values).
    :param intervals: List of lengths between each pair of adjacent colors.
    :param N: Number of entries in the map.
    """
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap

    x_sum = np.sum(intervals)
    cm = []
    for c1, c2, x in zip(colors[:-1], colors[1:], intervals):
        tmp_n = max(int(N*x/x_sum), 1)
        tmp_cm = LinearSegmentedColormap.from_list("", [c1, c2], N=tmp_n)
        cm.append(tmp_cm(np.linspace(0, 1, tmp_n)))

    cm = ListedColormap(np.concatenate(cm))
    return cm


def customize_palette(row_val, palette=None, size=50):
    mycolor = sns.color_palette('Greys', size) if palette is None else palette
    x = np.linspace(0, max(row_val), size)
    return [mycolor[pd.Series(abs(i-x)).idxmin()] for i in row_val]


def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    """
    From matplotlib examples.

    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)