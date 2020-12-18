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


def setAlpha(ax,a):
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
    import numpy as np
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap

    x_sum = np.sum(intervals)
    cm = []
    for c1, c2, x in zip(colors[:-1], colors[1:], intervals):
        tmp_n = max(int(N*x/x_sum), 1)
        tmp_cm = LinearSegmentedColormap.from_list("", [c1, c2], N=tmp_n)
        cm.append(tmp_cm(np.linspace(0, 1, tmp_n)))

    cm = ListedColormap(np.concatenate(cm))
    return cm
