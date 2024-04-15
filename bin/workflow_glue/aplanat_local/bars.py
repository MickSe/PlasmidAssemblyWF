"""Creation of bar-like plots."""

from bokeh.models import Range1d
from bokeh.models.ranges import FactorRange
from bokeh.plotting import figure
import numpy as np
import pandas as pd

from aplanat import util


@util.plot_wrapper
def single_hbar(values, classes, colors, **kwargs):
    """Create a plot with a single set of stacked horizontal bars.

    :param values: values to stack.
    :param classes: names of items.
    :param colors: colours for each item.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    After creating this plot, you may wish to change the following:

        * p.xaxis.axis_label: the axis label.
        * p.x_range: the x range of the plot.
    """
    if len(values) != len(classes) or len(values) != len(colors):
        raise ValueError(
            '`values`, `classes`, and `colors` must be of equal length.')
    # describe our data
    data = dict(zip(classes, ([x] for x in values)))
    data['y_value'] = ['']

    defaults = {
        'output_backend': 'webgl',
        'y_range': data['y_value'],
        'plot_height': 150, 'plot_width': 600,
        'toolbar_location': None, 'tools': 'hover',
        'tooltips': '$name: @$name'}
    defaults.update(kwargs)

    p = figure(**defaults)
    p.hbar_stack(
        classes, y='y_value',
        height=0.9, alpha=0.7,
        color=colors, source=data,
        legend_label=classes)

    # hide some plotting artefacts
    p.xgrid.grid_line_color = None
    p.x_range = Range1d(0, int(1.2*sum(values)))
    p.yaxis.visible = False
    return p


def simple_bar(
        groups, counts, colors=util.Colors.cerulean, **kwargs):
    """Create a simple barplot.

    :param groups: the grouping variable (the x-axis values).
    :param counts: the data for bars are drawn (the y-axis values).
    :param colors: the bar color.
    :param kwargs: kwargs for bokeh figure.

    """
    # see https://docs.bokeh.org/en/latest/docs/user_guide/categorical.html
    # for how boxplots can get complicated fast!
    defaults = {
        'output_backend': 'webgl',
        'plot_height': 300, 'plot_width': 600}
    defaults.update(kwargs)
    p = figure(
        x_range=groups,
        **defaults)
    p.vbar(
        x=groups, top=counts,
        fill_color=colors, line_color=colors, width=0.9)
    p.xgrid.grid_line_color = None
    p.y_range.start = 0
    return p


def simple_hbar(groups, counts, colors=util.Colors.cerulean, **kwargs):
    """Create a simple horizontal barplot.

    :param groups: the grouping variable (the y-axis values).
    :param counts: the data for bars are drawn (the x-axis values).
    :param colors: the bar color.
    :param kwargs: kwargs for bokeh figure.

    """
    defaults = {
        'output_backend': 'webgl',
        'plot_height': 300, 'plot_width': 600}
    defaults.update(kwargs)
    fig = figure(**defaults)
    y = list(range(len(groups)))
    fig.hbar(y, right=counts, height=0.5, line_color=colors, fill_color=colors)
    # Override the numerical labels with categories.
    # Setting labels this way, rather than with figure(y_range= ...)
    # ensures that category labels align with bars.
    fig.yaxis.ticker = y
    mapper = {k: v for (k, v) in zip(y, groups)}
    fig.yaxis.major_label_overrides = mapper

    return fig


@util.plot_wrapper
def boxplot_series(
        groups, values, xlim=(None, None), ylim=(None, None),
        **kwargs):
    """Create a (e.g. time-) series of boxplots for a variable.

    :param groups: the grouping variable (the x-axis values). The function
        will handle also non-numeric, categorical grouping variables (though
        the sort order is not controllable).
    :param values: the data for boxplots are drawn (the y-axis values).
    :param xlim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param ylim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    After creating this plot, you may wish to change the following:

        * p.xaxis.axis_label = 'Read Length / bases'
        * p.yaxis.axis_label = 'Number of reads'
    """
    df = pd.DataFrame(
        dict(group=groups, value=values))
    uniq = df.group.unique()
    # numeric or categorical
    if not np.issubdtype(uniq.dtype, np.number):
        x_range = FactorRange(factors=uniq)
    else:
        xlim = util.Limiter(util.pad(uniq)).fix(*xlim)
        x_range = Range1d(
            start=xlim.min, end=xlim.max, bounds=(xlim.min, xlim.max))

    # find the quartiles and IQR for each category
    groups = df.groupby('group')
    quantiles = groups.quantile([0.25, 0.5, 0.75])
    quantiles.index.names = ['group', 'quantile']
    uniq = groups.apply(lambda x: x.name).tolist()
    q1 = quantiles.xs(0.25, level='quantile')
    q2 = quantiles.xs(0.50, level='quantile')
    q3 = quantiles.xs(0.75, level='quantile')
    iqr = q3 - q1
    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr

    defaults = {
        "output_backend": "webgl",
        "height": 300, "width": 600}
    defaults.update(kwargs)

    ylim = util.Limiter().accumulate(df['value']).fix(*ylim)
    y_range = Range1d(
        start=ylim.min, end=ylim.max, bounds=(ylim.min, ylim.max))

    p = figure(**defaults, x_range=x_range, y_range=y_range)

    # stems
    p.segment(uniq, upper.value, uniq, q3.value, line_color="black")
    p.segment(uniq, lower.value, uniq, q1.value, line_color="black")

    # boxes
    for low, high in ((q2.value, q3.value), (q1.value, q2.value)):
        p.vbar(
            uniq, 0.8, low, high, line_color='black')

    return p
