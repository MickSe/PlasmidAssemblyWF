"""Creation of line plots."""

from aplanat import util
from aplanat.base import simple


@util.plot_wrapper
def line(
        x_datas, y_datas, names=None, colors=None,
        xlim=(None, None), ylim=(None, None),
        **kwargs):
    """Create a simple line plot.

    :param x_datas: a list of datasets, each item containing the
        data for that set.
    :params y_datas: a list of weight vectors (one for each dataset).
    :param names: the name of each dataset.
    :param colors: color used to plot each dataset.
    :param xlim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param ylim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    After creating this plot, you may wish to change the following:

        * p.xaxis.axis_label
        * p.yaxis.axis_label
    """
    return simple(
        x_datas, y_datas, names=names, colors=colors,
        xlim=xlim, ylim=ylim, style='line',
        **kwargs)


@util.plot_wrapper
def steps(
        x_datas, y_datas, names=None, colors=None,
        xlim=(None, None), ylim=(None, None), mode='before',
        **kwargs):
    """Create a simple step-line plot.

    :param x_datas: a list of datasets, each item containing the
        data for that set.
    :params y_datas: a list of weight vectors (one for each dataset).
    :param names: the name of each dataset.
    :param colors: color used to plot each dataset.
    :param xlim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param ylim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param mode: one of 'before', 'center', 'after' indicating where the step
        (change in y value) should be placed.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    After creating this plot, you may wish to change the following:

        * p.xaxis.axis_label
        * p.yaxis.axis_label
    """
    return simple(
        x_datas, y_datas, names=names, colors=colors,
        xlim=xlim, ylim=ylim, mode=mode, style='steps',
        **kwargs)
