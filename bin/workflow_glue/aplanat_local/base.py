"""Base functions for simple plots."""

from bokeh.models import Range1d
from bokeh.plotting import figure

from aplanat import util
from aplanat.util import Limiter


@util.plot_wrapper
def simple(
        x_datas, y_datas, names=None, colors=None,
        xlim=(None, None), ylim=(None, None), style='line',
        **kwargs):
    """Create a simple line or scatter plot.

    :param x_datas: a list of datasets, each item containing the
        data for that set.
    :params y_datas: a list of weight vectors (one for each dataset).
    :param names: the name of each dataset.
    :param colors: color used to plot each dataset.
    :param xlim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param ylim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param style: one of `line` or `points`.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    After creating this plot, you may wish to change the following:

        * p.xaxis.axis_label
        * p.yaxis.axis_label
    """
    mode = None
    if 'mode' in kwargs:
        mode = kwargs['mode']
        del kwargs['mode']
    # TODO: slacken this
    if not isinstance(x_datas, list) or not isinstance(y_datas, list):
        raise ValueError(
            "x_datas and y_datas should be lists, one item per line to plot.")
    items = len(x_datas)
    if names is None:
        names = [None] * items
    if colors is None:
        colors = ['blue'] * items
    if any(len(x) != items for x in (y_datas, names, colors)):
        raise IndexError(
            "Lengths of x_datas, y_datas, names, and colors should be equal.")

    defaults = {
        "output_backend": "webgl",
        "height": 300, "width": 600}
    defaults.update(kwargs)

    p = figure(**defaults)
    x_lim = Limiter()
    y_lim = Limiter()
    for x, y, name, color in zip(x_datas, y_datas, names, colors):
        # TODO: add colours and names
        x_lim.accumulate(x)
        y_lim.accumulate(y)
        kw = {}
        if name is not None:
            kw['legend_label'] = name
        if color is not None:
            kw['color'] = color
        if style == 'line':
            p.line(x=x, y=y, line_width=1.5, **kw)
        elif style == 'points':
            p.circle(x=x, y=y, alpha=0.4, **kw)
        elif style == 'steps':
            p.step(x=x, y=y, line_width=1.5, mode=mode, **kw)
        else:
            raise ValueError('Unknown plot style: "{}"'.format(style))
    x_lim.fix(*xlim)
    y_lim.fix(*ylim)

    # limit the display range
    p.x_range = Range1d(
        start=x_lim.min, end=x_lim.max, bounds=(x_lim.min, x_lim.max))
    p.y_range = Range1d(
        start=y_lim.min, end=y_lim.max, bounds=(y_lim.min, y_lim.max))
    return p
