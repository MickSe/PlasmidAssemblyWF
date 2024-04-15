"""Simple bokeh plotting API."""

import argparse
import importlib
import json
import warnings

from bokeh.colors import RGB
from bokeh.embed.util import OutputDocumentFor, standalone_docs_json
import bokeh.io as bkio
from bokeh.layouts import gridplot
from bokeh.plotting import Figure

__version__ = "0.6.10"

# we don't run a comprehensive test suite and mostly in notebooks,
# so show warnings all the time.
warnings.simplefilter('always', DeprecationWarning)


class Grid(list):
    """Create a list of lists of fixed length."""

    def __init__(self, width=4):
        """Initialize the Grid."""
        self.width = width
        self.append(list())

    def add(self, item):
        """Add an item to the grid.

        :param item: item to add

        """
        self[-1].append(item)
        if len(self[-1]) == self.width:
            self.append(list())

    def extend(self, items):
        """Add multiple items.

        :param items: items to add.

        """
        for item in items:
            self.add(item)


def all_children(plotlike):
    """Recursively find all Figures in a plotting layout.

    :param plotlike: a bokeh layout or plot.
    """
    my_children = []
    if hasattr(plotlike, 'children'):
        for child in plotlike.children:
            if isinstance(child, tuple):
                child = child[0]
            my_children.extend(all_children(child))
    if isinstance(plotlike, Figure):
        my_children.append(plotlike)
    return my_children


def show(plot, background=None):
    """Show a plot in a notebook.

    :param background: a fill colour for the plot background, hex string or RGB
        tuple.

    """
    children = all_children(plot)
    orig_colours = None
    if background is not None:
        if isinstance(background, tuple):
            background = RGB(*background)
        elif not (isinstance(background, str) and background.startswith('#')):
            raise TypeError(
                "`background` should be a RGB tuple or hex-colour string.")
        orig_colours = list()
        for child in children:
            orig_colours.append((
                child.background_fill_color,
                child.border_fill_color))
            child.background_fill_color = background
            child.border_fill_color = background
    bkio.output_notebook(hide_banner=True)
    bkio.show(plot)
    if background is not None:
        for colours, child in zip(orig_colours, children):
            child.background_fill_color = colours[0]
            child.border_fill_color = colours[1]


def json_item(plot, target=None, theme=None, always_new=True):
    """Export a plot to a JSON doc.

    :param plot: the Bokeh object to embed.
    :param target: a div id to embed the model into.
    :param theme: applies a specified theme.
    :param always_new: force creation of a new document.

    This is a reimplementation of bokeh.embed.json_item that exposes
    the `always_new` option and by default sets it to obtain a
    truly independent plot document.
    """
    with OutputDocumentFor(
            [plot], apply_theme=theme, always_new=always_new) as doc:
        doc.title = ""
        docs_json = standalone_docs_json([plot])

    doc_json = list(docs_json.values())[0]
    root_id = doc_json['roots']['root_ids'][0]

    return {'target_id': target, 'root_id': root_id, 'doc': doc_json}


def dump_json(plot, target=None, theme=None, always_new=True):
    """Create a JSON string representing a plot.

    :param plot: the Bokeh object to embed.
    :param target: a div id to embed the model into.
    :param theme: applies a specified theme.
    :param always_new: force creation of a new document.
    """
    data = json_item(plot, target=None, theme=None, always_new=always_new)
    return json.dumps(data)


def export_jsx(plot, fname):
    """Export plot to a JSX (react) file.

    :param plot: bokeh plot object.
    :param fname: export filename.
    """
    with open(fname, "w") as fh:
        fh.write("const plotJson = ")
        fh.write(dump_json(plot))
        fh.write("\n")
        fh.write("export default plotJson")


def grid(plots, ncol=4, display=True, **kwargs):
    """Show a grid of plots in a notebook.

    :param plots: a list of bokeh plots.
    :param ncol: number of columns.
    :param kwargs: kwargs for bokeh `gridplot`.

    """
    warnings.warn(
        "Please use `bokeh.layout.gridplot` directly with an optional "
        "call to `.show()`", DeprecationWarning)

    grid = Grid(width=ncol)
    grid.extend(plots)
    plot = gridplot(grid, **kwargs)
    if display:
        show(plot)
    else:
        return plot


def InfoGraphItems():
    """Cumulatively create items for an infographic."""
    from aplanat import graphics
    warnings.warn(
        "This class has been moved to `aplanat.graphics.InfoGraphItems()`",
        DeprecationWarning)
    return graphics.InfoGraphItems()


def infographic(items, **kwargs):
    """Create an infographic 'plot'."""
    from aplanat import graphics
    warnings.warn(
        "This function has been moved to `aplanat.graphics.infographic()`",
        DeprecationWarning)
    return graphics.infographic(items, **kwargs)


def cli():
    """Run aplanat entry point."""
    parser = argparse.ArgumentParser(
        'aplanat',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-v', '--version', action='version',
        version='%(prog)s {}'.format(__version__))

    subparsers = parser.add_subparsers(
        title='subcommands', description='valid commands',
        help='additional help', dest='command')
    subparsers.required = True

    # add reporting modules
    modules = [
        'demo', 'bcfstats', 'mapula', 'nextclade',
        'fastcat', 'simple', 'depthcoverage']
    for module in modules:
        mod = importlib.import_module('aplanat.components.{}'.format(module))
        p = subparsers.add_parser(module, parents=[mod.argparser()])
        p.set_defaults(func=mod.main)

    args = parser.parse_args()
    args.func(args)
