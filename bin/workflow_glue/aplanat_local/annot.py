"""Add annotations to plots."""

from bokeh.models import Span, Title
import numpy as np


def marker_vline(
        plot, x, label=None, color='black', width=1.5, text_baseline='bottom'):
    """Add a vertical line to a plot with optional label.

    :param plot: a bokeh plot instance.
    :param x: the x-axis value where the line should be added.
    :param label: text label.
    :param color: colour for line.
    :param width: line width.
    :param text_baseline: one of `top`, `middle`, or `bottom` - where to place
        the label with respect to the line. `top` corresponds to the right,
        `bottom` to the left of the drawn line.

    :returns: updated plot.

    """
    vline = Span(
        location=x, dimension='height',
        line_color=color, line_width=width)
    plot.renderers.extend([vline])
    if label is not None:
        ymax = plot.y_range.end
        if ymax is None:
            raise ValueError(
                "Please set plot.y_range before adding a vline with text.")
        plot.text(
            x=[x], y=[0.5 * ymax], angle=np.pi / 2, text=[label],
            text_baseline=text_baseline)
    return plot


def marker_hline(
        plot, y, label=None, color='black', width=1.5, text_baseline='bottom'):
    """Add a horizontal line to a plot with optional label.

    :param plot: a bokeh plot instance.
    :param y: the y-axis value where the line should be added.
    :param label: text label.
    :param color: colour for line.
    :param width: line width.
    :param text_baseline: one of `top`, `middle`, or `bottom` - where to place
        the label with respect to the line.

    :returns: updated plot.

    """
    hline = Span(
        location=y, dimension='width',
        line_color=color, line_width=width)
    plot.renderers.extend([hline])
    if label is not None:
        xmax = plot.x_range.end
        if xmax is None:
            raise ValueError(
                "Please set plot.x_range before adding a hline with text.")
        plot.text(
            x=[0.5 * xmax], y=[y], angle=0, text=[label])
    return plot


def subtitle(plot, subtitle):
    """Add a subtitle to a plot.

    :param subtitle: the title to add

    :returns: updated plot.
    """
    plot.add_layout(
        Title(text=subtitle, text_font_style='italic'), 'above')
    return plot
