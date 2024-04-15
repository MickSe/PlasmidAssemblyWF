"""Plotting of spatial/3D data."""

from bokeh.models import ColorBar, LinearColorMapper, Range1d
from bokeh.palettes import Blues9
from bokeh.plotting import figure
from bokeh.transform import linear_cmap, log_cmap
from bokeh.util.hex import hexbin
import numpy as np
import pandas as pd

from aplanat import util


@util.plot_wrapper
def heatmap(x, y, z, name=None, **kwargs):
    """Create a heatmap from three columns.

    :param x: x-axis coordinates.
    :param y: y-axis coordinates.
    :param z: z-axis coordinates (e.g. counts or frequencies).
    :param name: title for z-axis.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    """
    data = {'x': x, 'y': y, 'z': z}

    defaults = {
        "output_backend": "webgl",
        "height": 300, "width": 600}
    defaults.update(kwargs)
    p = figure(**defaults)

    # define a colour map
    colors = Blues9[::-1]
    mapper = LinearColorMapper(
        palette=colors, low=min(data['z']), high=max(data['z']))

    # plot heatmap rectangles
    p.rect(
        source=data, x="x", y="y", width=1, height=1,
        fill_color={'field': 'z', 'transform': mapper},
        line_color=None)

    # add colour scale
    color_bar = ColorBar(
        title=name, color_mapper=mapper, label_standoff=10,
        location=(0, 0))
    p.add_layout(color_bar, 'right')

    # hide some plotting artefacts
    p.x_range.range_padding = 0
    p.y_range.range_padding = 0
    p.xaxis.visible = False
    p.yaxis.visible = False
    # set up the axes
    xbounds = util.pad(data['x'])
    ybounds = util.pad(data['y'])
    p.x_range = Range1d(start=xbounds[0], end=xbounds[1], bounds=xbounds)
    p.y_range = Range1d(start=ybounds[0], end=ybounds[1], bounds=ybounds)
    return p


@util.plot_wrapper
def heatmap2(
        x, y, name=None, x_bins=50, y_bins=50, log=False,
        xlim=(None, None), ylim=(None, None), zlim=(None, None),
        **kwargs):
    """Create a heatmap from two columns, using integrate binning.

    In contrast to `heatmap` which takes coordinates and counts, this function
    takes only coordinates and will count observations within bins.

    :param x: x-axis coordinates.
    :param y: y-axis coordinates.
    :param name: title for z-axis.
    :param x_bins: number of bins in x-direction.
    :param y_bins: number of bins in y-direction.
    :param xlim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param ylim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param zlim: tuple for colour map limits (start, end). A value None will
        trigger calculation from the data.
    :param log: use a log-scaled colour map.

    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    """
    hist, x_edges, y_edges = np.histogram2d(x, y, bins=(x_bins, y_bins))
    table = np.empty(
        x_bins * y_bins,
        dtype=[('row', float), ('col', float), ('value', float)])
    table['value'] = hist.flatten()
    table['row'] = np.repeat(x_edges[:x_bins], y_bins)
    table['col'] = np.tile(y_edges[:y_bins], x_bins)
    table = pd.DataFrame(table)

    defaults = {
        "output_backend": "webgl",
        "height": 300, "width": 600}
    defaults.update(kwargs)
    p = figure(**defaults)
    # define a colour map
    cmap = log_cmap if log else linear_cmap
    cmin, cmax = zlim
    if cmin is None:
        cmin = min(table['value'])
    if cmax is None:
        cmax = max(table['value'])
    mapper = cmap('counts', 'Viridis256', cmin, cmax)
    # plot heatmap rectangles
    width = (x_edges[-1] - x_edges[0]) / x_bins
    height = (y_edges[-1] - y_edges[0]) / y_bins
    p.rect(
        x="row", y="col", source=table,
        height=height * 1.05, width=width * 1.05,  # fudge to hide artefacts
        fill_color={'field': 'value', 'transform': mapper['transform']},
        line_color=None)
    # add colour scale
    color_bar = ColorBar(
        title=name, color_mapper=mapper['transform'], label_standoff=10,
        location=(0, 0))
    p.add_layout(color_bar, 'right')
    # hide some plotting artefacts
    p.x_range.range_padding = 0
    p.y_range.range_padding = 0
    # set up the axes
    xbounds = util.pad(x_edges[:x_bins])
    ybounds = util.pad(y_edges[:y_bins])
    x_lim = util.Limiter().accumulate(xbounds).fix(*xlim)
    y_lim = util.Limiter().accumulate(ybounds).fix(*ylim)
    p.x_range = Range1d(
        start=x_lim.min, end=x_lim.max, bounds=(x_lim.min, x_lim.max))
    p.y_range = Range1d(
        start=y_lim.min, end=y_lim.max, bounds=(y_lim.min, y_lim.max))
    return p


def hexmap(
        x, y, name=None, log=False, match_aspect=True, kernel_size=0.1,
        xlim=(None, None), ylim=(None, None), **kwargs):
    """Create a heatmap from two columns.

    :param x: x-axis coordinates.
    :param y: y-axis coordinates.
    :param name: title for z-axis.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh plot.

    """
    bins = hexbin(x, y, kernel_size)

    defaults = {
        "output_backend": "webgl",
        "height": 400, "width": 500,
        "background_fill_color": "#440154"}
    defaults.update(kwargs)
    p = figure(**defaults)

    cmap = log_cmap if log else linear_cmap
    mapper = cmap('counts', 'Viridis256', 0, max(bins.counts))
    p.hex_tile(
        q="q", r="r", size=kernel_size, line_color=None, source=bins,
        fill_color=mapper)

    color_bar = ColorBar(
        title="observations", color_mapper=mapper['transform'],
        label_standoff=10, location=(0, 0))
    p.add_layout(color_bar, 'right')
    p.grid.visible = False

    x_lim = util.Limiter().accumulate(x).fix(*xlim)
    y_lim = util.Limiter().accumulate(y).fix(*ylim)
    # limit the display range
    p.x_range = Range1d(
        start=x_lim.min, end=x_lim.max, bounds=(x_lim.min, x_lim.max))
    p.y_range = Range1d(
        start=y_lim.min, end=y_lim.max, bounds=(y_lim.min, y_lim.max))
    return p, bins
