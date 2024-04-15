"""Bio specific plots, could be generalised."""

from io import StringIO

from bokeh.models import Range1d
from bokeh.models.tickers import FixedTicker
from bokeh.plotting import figure
import numpy as np
import pandas as pd

from aplanat import util

_chrom_data_ = """
chrom length
1 249250621
2 243199373
3 198022430
4 191154276
5 180915260
6 171115067
7 159138663
8 146364022
9 141213431
10 135534747
11 135006516
12 133851895
13 115169878
14 107349540
15 102531392
16 90354753
17 81195210
18 78077248
19 59128983
20 63025520
21 48129895
22 51304566
X 155270560
Y 59373566
"""
_chrom_data_ = pd.read_csv(StringIO(_chrom_data_), delimiter=" ")


@util.plot_wrapper
def karyotype(
        x_datas, y_datas, names=None, colors=None, alpha=0.2,
        chrom_data=_chrom_data_, **kwargs):
    """Create a heatmap from three columns.

    :param x: x-axis coordinates (chromosome position)
    :param y: y-axis coordinates (chromosome)
    :param name: title for z-axis.
    :param colors: colours to plot tics for each data series.
    :param alpha: alpha level (transparency) with which to plot tics.
    :params chrom_data: a dataframe (or ndarray) with fields
        `chrom`, and `length` used to plot chromosome outline
        rectangles. Chromosomes will be plotted vertically (top to bottom)
        in the order given by this input.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    """
    if names is None:
        names = [None] * len(x_datas)
    if colors is None:
        colors = [None] * len(x_datas)
    defaults = {
        "output_backend": "webgl",
        "height": 300, "width": 600}
    defaults.update(kwargs)
    p = figure(**defaults)
    width = 0.8
    chrom_order = pd.Series(range(1, len(chrom_data) + 1))
    p.rect(
        chrom_data["length"] // 2, chrom_order,
        chrom_data["length"], width,
        fill_color='white', line_color='black')

    plot_order = dict(zip(chrom_data['chrom'], chrom_order))
    for x, y, name, color in zip(x_datas, y_datas, names, colors):
        y = np.array([plot_order[i] for i in y])
        kw = {'alpha': alpha}
        if name is not None:
            kw['legend_label'] = name
        if color is not None:
            kw['color'] = color
        p.segment(
            x, y - width / 2, x, y + width / 2, **kw)

    # set up the axes
    p.xaxis.formatter.use_scientific = False
    p.yaxis.ticker = FixedTicker()
    p.yaxis.ticker.ticks = chrom_order
    p.yaxis.major_label_overrides = {v: k for k, v in plot_order.items()}
    p.x_range = Range1d(
        start=0, end=max(chrom_data['length']),
        bounds=(0, max(chrom_data['length'])))
    miny, maxy = 0, len(chrom_data) + 1
    p.y_range = Range1d(
        start=miny, end=maxy, bounds=(miny, maxy))
    return p
