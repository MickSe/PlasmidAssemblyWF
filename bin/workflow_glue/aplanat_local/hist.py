"""Histogram-like plots."""

from bokeh.models import Range1d
from bokeh.plotting import figure
import numpy as np

from aplanat import util
from aplanat.util import Limiter


@util.plot_wrapper
def histogram(
        datas, weights=None, names=None, colors=None,
        normalize=False, bins=30, binwidth=None,
        xlim=(None, None), ylim=(None, None), style='bars',
        **kwargs):
    """Create a histogram plot.

    :param datas: a list of datasets, each item containing the data
        for that set.
    :params weights: a list of weight vectors (one for each dataset).
    :param names: the name of each dataset.
    :param colors: color used to plot each dataset.
    :param normalize: normalize histogram counts by total weight (across all
        datasets).
    :param bins: number of histogram bins.
    :param binwidth: overrides `bins`, provide a precise binwidth.
    :param xlim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param ylim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param style: one of: `bars`, `lines`. `bars` will plot the classical
        bar histogram with a bar depicting counts per fixed interval. `lines`
        adapts this to plot a line per dataset rather than bars.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    After creating this plot, you may wish to change the following:

        * p.xaxis.axis_label
        * p.yaxis.axis_label
    """
    # TODO: slacken this
    if not isinstance(datas, list):
        raise ValueError(
            "datas and weights should be lists, one item per line to plot.")
    items = len(datas)
    if weights is None:
        weights = [None] * items
    if names is None:
        names = [None] * items
    if colors is None:
        colors = ['blue'] * items
    if any(len(x) != items for x in (weights, names, colors)):
        raise IndexError(
            "Lengths of datas, weights, names, and colors should be equal.")

    x_lim = Limiter()
    y_lim = Limiter()
    for data in datas:
        x_lim.accumulate(data)
    if binwidth is not None:
        bins = np.arange(x_lim.min, x_lim.max + 1e-10, binwidth)
    else:
        bins = np.linspace(x_lim.min, x_lim.max, num=bins)
    if normalize:
        total_weight = sum(sum(x) for x in weights)

    defaults = {
        "output_backend": "webgl",
        "height": 300, "width": 600}
    defaults.update(kwargs)

    p = figure(**defaults)
    for data, weight, name, color in zip(datas, weights, names, colors):
        counts, edges = np.histogram(data, weights=weight, bins=bins)
        if normalize:
            counts = counts.astype(float) / total_weight
        y_lim.accumulate(counts)
        kw = {}
        if name is not None:
            kw['legend_label'] = name
        if color is not None:
            kw['color'] = color
        if style == 'bars':
            p.quad(
                top=counts, bottom=0, left=edges[:-1], right=edges[1:],
                alpha=0.6, **kw)
        elif style == 'line':
            p.line(x=edges[:-1], y=counts, line_width=1.5, **kw)
    x_lim.fix(*xlim)
    y_lim.fix(*ylim)

    # limit the display range
    p.x_range = Range1d(
        start=x_lim.min, end=x_lim.max, bounds=(x_lim.min, x_lim.max))
    p.y_range = Range1d(
        start=y_lim.min, end=y_lim.max, bounds=(y_lim.min, y_lim.max))
    return p
