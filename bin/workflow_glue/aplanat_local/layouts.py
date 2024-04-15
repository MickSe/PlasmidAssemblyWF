"""Creation of more complex plotting layouts."""

import itertools

from bokeh.layouts import gridplot, layout, Row, Spacer
from bokeh.models import Range1d
from bokeh.models.annotations import Label
from bokeh.plotting import figure

from aplanat import util


def facet_grid(
        df, aes, plot_func,
        xlim=(None, None), ylim=(None, None),
        facet=(None, None),
        x_axis_label="", y_axis_label="",
        x_facet_heading="", y_facet_heading="",
        transform=(None, None), link_axes=True,
        **kwargs):
    """Mainly for the creation of facet-grid plots from a Dataframe.

    :param df: pandas dataframe
    :param aes: a dictionary with keys 'x', 'y', 'col' corresponding to the
        data to plot along the x-axis, y-axis, and optional colouring. These
        should be column names of the input dataframe.
    :param plot_func: aplanat plotting function to use.
    :param xlim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param ylim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.
    :param facet: a tuple indicating the variables by which to create (x, y)
        facets of the data.
    :param x_axis_label: label text for the x axis.
    :param y_axis_label: label text for the y axis
    :param x_facet_heading: facet name for the columns (x facet)
    :param y_facet_heading: facet name for the rows (y facet)
    :param transform: a tuple of functions to apply to the (x, y) data of each
        facet before plotting.
    :param link_axes: link axes of subplots such that panning and zooming
        effects all plots.
    :param kwargs: kwargs for bokeh figure.

    :returns: a bokeh layout.
    """
    if 'x' not in aes or 'y' not in aes:
        raise KeyError("'aes' argument must contain keys: 'x', 'y'.")
    # find unique values of facets
    facet_x, facet_y = facet
    facet_x_values, facet_y_values = [None], [None]
    # and colours
    col_values = [None]
    if 'col' in aes:
        col_values = df[aes['col']].unique()
        cols = util.choose_palette(len(col_values))
        cols_dict = dict(zip(col_values, cols))

    if facet_x is not None:
        facet_x_values = df[facet_x].unique()
    if facet_y is not None:
        facet_y_values = df[facet_y].unique()
    facet_products = list(itertools.product(
        enumerate(facet_x_values), enumerate(facet_y_values)))

    plots = list()
    for (px, fx), (py, fy) in facet_products:
        d = df
        if fx is not None and fy is not None:
            d = df.loc[(df[facet_x] == fx) & (df[facet_y] == fy)]
        elif fx is not None:
            d = df.loc[df[facet_x] == fx]
        elif fy is not None:
            d = df.loc[df[facet_y] == fy]

        x_data = list()
        y_data = list()
        colors = list()
        for col in col_values:
            if col is not None:
                dcol = d.loc[d[aes['col']] == col]
            else:
                dcol = d
            if len(dcol) == 0:
                continue
            colors.append(cols_dict[col])
            # create plot vectors with transform if required
            x = dcol[aes['x']]
            if transform[0] is not None:
                x = transform[0](x)
            x_data.append(x)
            if 'y' in aes:
                y = dcol[aes['y']]
                if transform[1] is not None:
                    y = transform[1](y)
                y_data.append(y)
        # make the plot
        if len(x_data) == 0:
            continue  # no data for this facet
        plot = plot_func(
            x_data, y_data, colors=colors,
            xlim=xlim, ylim=ylim,
            **kwargs)
        plots.append(((px, py), plot))

    # fix limits of all plots to the same
    for attr in ['x_range', 'y_range']:
        zlims = [
            (getattr(p, attr).start, getattr(p, attr).end)
            for _, p in plots]
        zmin, zmax = min(z[0] for z in zlims), max(z[1] for z in zlims)
        zrange = Range1d(start=zmin, end=zmax, bounds=(zmin, zmax))
        for _, p in plots:
            if not link_axes:
                # nb: copy.deepcopy doesn't do as intended
                zrange = Range1d(start=zmin, end=zmax, bounds=(zmin, zmax))
            setattr(p, attr, zrange)

    # make the main panel of plots
    ncols, nrows = len(facet_x_values), len(facet_y_values)
    grid = []
    for _ in range(nrows):
        grid.append([None]*ncols)
    for (col, row), p in plots:
        grid[row][col] = p
    plotgrid = gridplot(grid)
    # figure out the height and width that were used
    bh, bw = plots[0][1].plot_height, plots[0][1].plot_width

    # get labels and headers
    x_header, x_axis_label, *x_labels = make_facet_labels(
        x_facet_heading, facet_x_values, x_axis_label,
        base_height=bh, base_width=bw)
    y_header, y_axis_label, *y_labels = make_facet_labels(
        y_facet_heading, facet_y_values, y_axis_label, facet='y',
        base_height=bh, base_width=bw)

    # add legend
    # bokeh doesn't support legends outside a figure
    # https://github.com/bokeh/bokeh/issues/7607
    # so make a whole new figure
    leg_plot = plot_func(
        [[0, 0.1]]*len(col_values), [[0, 0.1]]*len(col_values),
        colors=cols, names=col_values,
        title=None, toolbar_location=None,
        height=bh*len(facet_y_values), width=150)
    leg_plot.axis.visible = False
    leg_plot.grid.visible = False
    leg_plot.outline_line_color = None
    # fudge plot to be outside where data is
    leg_plot.x_range = Range1d(start=-1.0, end=-0.5)
    leg_plot.y_range = Range1d(start=-1.0, end=-0.5)

    # create the plot layout, row by row.
    plot_layout = []
    # toobox, aligned right
    plot_layout.append(
        [Spacer(sizing_mode='stretch_width'), plotgrid.children[0]])
    if len(facet_x_values) > 1:
        # x facet header
        plot_layout.append(
            [Spacer(width=40), x_header, Spacer(sizing_mode='stretch_width')])
        # x facet labels
        plot_layout.append(
            [Spacer(width=40)] + x_labels +
            [Spacer(sizing_mode='stretch_width')])
    # y axis label, plots, y facet labels, y facet header
    row = [
        y_axis_label, Row(plotgrid.children[1])]
    if len(facet_y_values) > 1:
        row.extend([y_labels, y_header])
    if len(col_values) > 1:
        row.append(leg_plot)
    plot_layout.append(row)
    # x axis label
    plot_layout.append(
        [Spacer(width=40), x_axis_label, Spacer(sizing_mode='stretch_width')])
    # make plot
    plot = layout(plot_layout)
    return plot


def make_facet_labels(
        facet_name, facet_labels, axis_label,
        facet='x', base_height=300, base_width=400):
    """Create labels for a facet grid plot.

    :param facet_name: the header name for the facet.
    :param facet_labels: the names of the factors of the facet.
    :param axis_label: the axis label (as would appear on a single plot)
    :param facet: either 'x' or 'y'.
    :param base_height: the height of the corresponding plots.
    :param base_width: the width of the corresponding plots.

    :returns: [facet header, axis label, facet 1 label, ...]

    """
    plots = list()
    h, w = 40, base_width
    head_h, head_w = 40, len(facet_labels) * w
    angle = 0
    axis_angle = 0
    plots = list()
    if facet == 'y':
        h, w = base_height, 40
        head_h, head_w = len(facet_labels) * h, 40
        angle = -3.14/2
        axis_angle = 3.14/2

    def _label_plot(label, h, w, angle=angle, **kwargs):
        # a figure that simply holds a label
        p = figure(
            output_backend='webgl',
            plot_height=h, plot_width=w,
            title=None, toolbar_location=None)
        p.axis.visible = False
        p.grid.visible = False
        p.outline_line_color = None
        p.toolbar.active_drag = None
        p.x_range = Range1d(start=0.0, end=1.0)
        p.y_range = Range1d(start=0.0, end=1.0)
        if label is None:
            label = ""
        p.add_layout(
            Label(
                x=0.5, y=0.5, text=label, text_align='center',
                text_baseline='middle', angle=angle, **kwargs))
        return p

    # create header and axis label
    plots.append(_label_plot(facet_name, h=head_h, w=head_w))
    plots.append(
        _label_plot(
            axis_label, h=head_h, w=head_w, angle=axis_angle,
            text_font_style='italic'))
    # label for each facet (column/row)
    for label in facet_labels:
        plots.append(_label_plot(label, h=h, w=w))
    return plots
