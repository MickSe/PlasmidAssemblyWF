"""Plots which are more graphics than plots."""

import os
import tempfile

from bokeh.layouts import gridplot
from bokeh.models import Range1d
from bokeh.models.annotations import Label
from bokeh.plotting import figure
from icon_font_to_png.icon_font import IconFont
import numpy as np
from PIL import Image
import pkg_resources
from si_prefix import si_format
import sigfig


class IconRGBA:
    """Wrapper to icon_font_to_png for multiple ttfs."""

    def __init__(self):
        """Initialize the class."""
        self.css_file = pkg_resources.resource_filename(
            __package__, 'data/fontawesome.css')
        self.ttf_files = [
            pkg_resources.resource_filename(__package__, 'data/{}'.format(x))
            for x in ('fa-regular-400.ttf', 'fa-solid-900.ttf')]
        self.icon_sets = [
            IconFont(self.css_file, ttf) for ttf in self.ttf_files]
        self.icons = dict()
        for i, icon_set in enumerate(self.icon_sets):
            self.icons.update(
                {k: (i, icon_set) for k in icon_set.css_icons.keys()})
        self.tmpdir = tempfile.TemporaryDirectory()

    def rgba(self, icon, size, color='black', scale='auto', uint32=True):
        """Create RGBA array for icon.

        :param icon: valid icon name
        :param size: icon size in pixels
        :param color: color name or hex value
        :param scale: scaling factor between 0 and 1,
            or 'auto' for automatic scaling
        :param uint32: return a two-dimension uint32-packed array (as
            required by bokeh image_rgba).
        """
        fname = icon + '.png'
        tmpfile = os.path.join(self.tmpdir.name, fname)
        try:
            icon_set = self.icon_sets[self.icons[icon][0]]
        except KeyError:
            raise KeyError('Unknown icon: {}.'.format(icon))
        icon_set.export_icon(
            icon, size, color, scale, filename=fname,
            export_dir=self.tmpdir.name)
        if not os.path.isfile(tmpfile):
            raise RuntimeError('Image not produced.')
        image = np.ascontiguousarray(
            np.array(Image.open(tmpfile))[::-1, :, :])
        if uint32:
            image = image.view(dtype=np.uint32).reshape(size, size)
        return image


fa_icons = IconRGBA()


class InfoGraphItems(dict):
    """Helper class to cumulatively create items for and infographic."""

    def __init__(self, *args, **kwargs):
        """Initialize the helper, and bootstrap fontawesome."""
        super().__init__(*args, **kwargs)

    def append(self, label, value, icon, unit=''):
        """Add an item.

        :param label: infographic item label.
        :param value: numerical value of headline number (without SI units).
        :param icon: font-awesome icon to use.
        :param unit: additional suffix after SI unit suffix, e.g. "bases".

        """
        try:
            fa_icons.icons[icon]
        except KeyError:
            raise KeyError("'{}' is not a known icon.".format(icon))
        self[label] = (label, value, icon, unit)

    def extend(self, items):
        """Add multiple items at once.

        :param items: iterable of 4-tuples, as required by `.add()`.
        """
        for i in items:
            self.append(*i)


def infographic(items, **kwargs):
    """Create and infographic 'plot'.

    :param items: 3-tuples of (label, value, unit, icon); the label should be
        a one or two word description, the value the headline number, and the
        icon the name of a fontawesome icon. If `value` is  numeric, it
        will be normalised by use of an SI suffix for display after which
        `unit` will be appended. If `value` is a string it will be used as is.
    :param kwargs: kwargs for bokeh gridplot.

    ..note:: If `bootstrap_fontawesome` has not already been called, the
        function will load the required fonts, however they will not display
        the first time an Jupyter labs cell is run. If using the
        `InfoGraphItems` helper class, this wrinkle will be taken care of
        provided the helper is initiated in a previous cell.

    """
    plots = list()
    seen = set()
    for label, value, icon, unit in items:
        if label in seen:
            continue
        if not isinstance(value, str):
            if unit == '%':
                value = "{}%".format(sigfig.round(value, 3))
            else:
                value = si_format(value) + unit
        seen.add(label)
        width, height = 175, 100
        aspect = height / width
        p = figure(
            output_backend='webgl',
            plot_width=width, plot_height=height,
            title=None, toolbar_location=None)
        p.axis.visible = False
        p.grid.visible = False
        p.outline_line_color = None
        p.rect([0.5], [0.5], [1.0], [1.0], fill_color="#6205A8")
        p.x_range = Range1d(start=0.1, end=0.9, bounds=(0.1, 0.9))
        p.y_range = Range1d(start=0.1, end=0.9, bounds=(0.1, 0.9))
        p.add_layout(
            Label(
                x=0.15, y=0.45, text=value, text_color="#F3EAF9",
                text_font_size="24px"))
        p.add_layout(
            Label(
                x=0.15, y=0.2, text=label, text_color="#F3EAF9",
                text_font_size="16px"))
        image = fa_icons.rgba(icon, 75, color='#BC61FF')
        p.image_rgba(image=[image], x=0.6, y=0.4, dw=0.25, dh=0.25/aspect)
        plots.append(p)
    defaults = {'toolbar_location': None, "ncols": len(items)}
    defaults.update(kwargs)
    return gridplot(plots, **defaults)
