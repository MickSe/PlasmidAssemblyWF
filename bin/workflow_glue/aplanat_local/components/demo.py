"""Aplanat demo."""

import argparse
from copy import deepcopy

from bokeh.layouts import gridplot
import numpy as np
import pandas as pd

import aplanat
from aplanat import \
    annot, bars, bio, graphics, hist, lines, points, spatial, util
from aplanat.report import HTMLReport


def main(args):
    """Aplanat demo entry point."""
    util.set_basic_logging()
    logger = util.get_named_logger("Aplanat Demo")

    x = np.random.normal(size=2000)
    y = np.random.normal(size=2000)
    sorted_xy = [np.sort(x), np.sort(x)]

    # Start a report
    logger.info("Start a report")
    report = HTMLReport(
        title="Aplanat Demo",
        lead="A brief demonstration of Aplanat's API.",
        require_keys=False)  # set True to require keys on item addition

    # The report is an ordered dictionary, so we can add placeholders to
    # delay the addition of items. Markdown can be used for text items.
    report.markdown("placeholder", key="simple preamble")

    # We can also add additional sections. These function as their own
    # mini-report and each section is rendered in its entirety in the order it
    # was added.  (The main HTMLReport is logically the first section).
    logger.info("Adding gallery")
    gallery = report.add_section(key="additional_section")

    # Adding a plot
    logger.info("Adding points plot")
    pnt_plot = points.points(sorted_xy, sorted_xy[::-1])
    report.plot(pnt_plot, key="points_plot")

    # Test we can output a plot to json in a different doc/layout
    # using embed.json_item() fails in this case
    logger.info("Testing json export")
    report.markdown("A simple line plot.")
    line_plot = lines.line(sorted_xy, sorted_xy[::-1])
    json_plot = aplanat.dump_json(line_plot)
    logger.info("Some json: {}...".format(json_plot[:40]))
    report.plot(gridplot([[line_plot]]), key="lines_plot")

    # Using placeholder is more explicit, and checks will be made before
    # rendering that the item has be assigned a real value
    report.placeholder("histogram preamble")

    # There's no need to provide key (unless require_keys is set). Items added
    # without a key cannot be replaced however.
    logger.info("Adding histogram")
    h = hist.histogram([x - 1, y + 1], colors=['red', 'green'])
    h = annot.marker_vline(
        h, np.mean(x) - 1, label='x values - 1', color='red',
        text_baseline='bottom')
    h = annot.marker_vline(
        h, np.mean(y) + 1, label='y values + 1', color='green',
        text_baseline='top')
    report.plot(h)

    # To delete an item, just delete the key
    logger.info("Deleting an item for fun")
    report.markdown("Garbage", key='garbage')
    del report['garbage']

    # Add more plots
    logger.info("Adding heatmap")
    report.placeholder("heatmap preamble")
    report.plot(spatial.heatmap2(x, y))

    # Add a data table
    logger.info("Adding table")
    report.placeholder("table preamble")
    df = pd.DataFrame({'x': x, 'y': y})
    report.markdown(
        "Here's a table:")
    df = df[0:5]
    report.table(
        df, key='Table with auto_height', pagination='false',
        sortable='false')
    report.markdown(
        "Here's a second table with more data, pagination, and searchable")
    df = pd.DataFrame({'x': y, 'y': x})
    report.table(df, key='Table with more data and pagination', height=200)

    # boxplot series with discretised x
    logger.info("Adding boxplot")
    report.placeholder("boxplot preamble")
    x_discrete = np.around(x, 0)
    x_str = [str(x) for x in np.abs(x_discrete)]
    report.plot(
        gridplot([
            bars.boxplot_series(x_discrete, y, width=300, title='continuous'),
            bars.boxplot_series(x_str, y, width=300, title='categorical')],
            ncols=2))

    # Gallery
    gallery.placeholder("gallery preamble")
    logger.info("Adding some infographics")
    exec_summary = graphics.InfoGraphItems()
    exec_summary.append('Example', 0.0051, 'angle-up', '%')
    exec_summary.append('Total reads', 1000000, 'angle-up', '')
    exec_summary.append('Total yield', 1e9, 'signal', 'b')
    exec_summary.append('Mean read length', 50e3, 'align-center', 'b')
    exec_summary.append('Mean qscore (pass)', 14, 'thumbs-up', '')
    plot = graphics.infographic(exec_summary.values())
    gallery.plot(plot)

    logger.info("Adding karyogram")
    chrom_data = deepcopy(bio._chrom_data_)
    chrom_data['chrom'] = chrom_data['chrom'].apply(lambda x: 'chr' + x)
    chroms = np.random.choice(chrom_data['chrom'], size=10000)
    positions = list()
    for chrom in chroms:
        length = chrom_data.loc[chrom_data['chrom'] == chrom, 'length']
        positions.append(np.random.randint(length)[0])
    plot = bio.karyotype([positions], [chroms], chrom_data=chrom_data)
    gallery.plot(plot)

    report.placeholder("failed plot preamble")
    logger.info(
        "Adding failed plot -- an exception trace will be logged below.")
    group = ("x1", "x2", "x3", 4, "x5")
    count = ('wrong data')

    @util.plot_wrapper
    def decorator_example(groups, counts):
        """Exemplify a decorator."""
        p = bars.simple_bar(groups, counts)
        return p

    plot = decorator_example(group, count)
    report.plot(plot)
    logger.info("End of plot that will fail.")

    # Trying to render now will raise ValueError because the placeholders are
    # not filled in
    logger.info("Checking empty pladeholder validation")
    try:
        report.render()
    except ValueError:
        logger.exception("Caught error as expected:")
    logger.info("End of expected exception log.")

    # Fill in the placeholders
    report.markdown("""
    ### Simple plots

    A simple points plot:
    """, "simple preamble")

    report.markdown("""
    ### Histograms

    Multi-variate histograms:
    """, "histogram preamble")

    report.markdown("""
    ### Heatmaps

    Heatmaps from `x-y` values:
    """, "heatmap preamble")

    report.markdown("""
    ### Tables

    Data tables:
    """, "table preamble")

    report.markdown("""
    ### Boxplots

    A series of boxplots using either a continuous or categorical grouping:
    """, "boxplot preamble")

    report.markdown("""
    ### Failed Plot

    When plot fails:
    """, "failed plot preamble")

    gallery.markdown("""
    ### Gallery

    Assortment of possibilities:
    """, "gallery preamble")

    # write the output, implicitely calls .render()
    logger.info("Rendering report")
    report.write(args.output)


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(
        "aplanat demo",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--output", default="aplanat_demo_report.html",
        help="Output HTML file.")
    return parser
