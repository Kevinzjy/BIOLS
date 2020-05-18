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
