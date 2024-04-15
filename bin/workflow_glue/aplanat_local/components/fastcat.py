"""Report components for displaying information from bcftools stats."""

import argparse

from bokeh.layouts import layout
import numpy as np
import pandas as pd

from aplanat_local import annot, bars, hist
from aplanat_local.report import _maybe_new_report, HTMLReport
from aplanat_local.util import Colors, read_files


def read_length_plot(
        seq_summary, min_len=None,
        max_len=None, xlim=(0, None)):
    """Create a read length plot.

    :param seq_summary: summary data from fastcat.
    :param min_len: minimum length.
    :param max_len: maximum length.
    :param xlim: tuple for plotting limits (start, end). A value None will
        trigger calculation from the data.

    The minimum and maximum lengths are used only to annotate the plot
    (not filter the data).
    """
    total_bases = seq_summary['read_length'].sum()
    try:
        mean_length = total_bases / len(seq_summary)
    except ZeroDivisionError:
        mean_length = 0
    median_length = np.median(seq_summary['read_length'])
    datas = [seq_summary['read_length']]
    length_hist = hist.histogram(
        datas, colors=[Colors.HES_purple_d], bins=100,
        title="Read length distribution.",
        x_axis_label='Read Length / bases',
        y_axis_label='Number of reads',
        xlim=xlim)
    if min_len is not None:
        length_hist = annot.marker_vline(
            length_hist, min_len,
            label="Min: {}".format(min_len), text_baseline='bottom',
            color='grey')
    if max_len is not None:
        length_hist = annot.marker_vline(
            length_hist, max_len,
            label="Max: {}".format(max_len), text_baseline='top')
    length_hist = annot.subtitle(
        length_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_length, median_length))
    return length_hist


def read_quality_plot(seq_summary):
    """Create read quality plot.

    :param seq_summary: summary data from fastcat.
    """
    datas = [seq_summary['mean_quality']]
    mean_q, median_q = np.mean(datas[0]), np.median(datas[0])
    q_hist = hist.histogram(
        datas, colors=[Colors.HES_purple_d], bins=100,
        title="Read quality score",
        x_axis_label="Quality score",
        y_axis_label="Number of reads",
        xlim=(4, 25))
    q_hist = annot.subtitle(
        q_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_q, median_q))
    return q_hist


def sample_read_counts(seq_summary, min_len=None, max_len=None):
    """Create bar plot counting unique samples.

    :param seq_summary: summary data from fastcat.
    :param min_len: minimum length.
    :param max_len: maximum length.
    """
    if min_len is not None:
        seq_summary = seq_summary.loc[
            (seq_summary['read_length'] > min_len)]
    if max_len is not None:
        seq_summary = seq_summary.loc[
            (seq_summary['read_length'] < max_len)]
    sample_counts = (
        pd.DataFrame(seq_summary['sample_name'].value_counts())
        .sort_index()
        .reset_index()
        .rename(
            columns={'index': 'sample', 'sample_name': 'count'})
        )

    title = 'Number of reads per barcode'
    if min_len is not None or max_len is not None:
        t0, t1 = "", ""
        if min_len is not None:
            t0 = "{} < ".format(min_len)
        if max_len is not None:
            t1 = " > {}".format(max_len)
        title += " (filtered by {}length{}).".format(t0, t1)

    plot = bars.simple_bar(
        sample_counts['sample'].astype(str), sample_counts['count'],
        colors=[Colors.cerulean]*len(sample_counts),
        title=title,
        plot_width=None)
    plot.xaxis.major_label_orientation = 3.14/2
    return plot


_full_report_header = """
### Read summary

The following tables and figures are derived from
the output of [fastcat](https://github.com/epi2me-labs/fastcat).
"""


def full_report(
        stats, header=_full_report_header, report=None,
        sample_counts=False, min_len=None, max_len=None):
    """Create a report section from the output of fastcat.

    :param stats: one or more outputs from `fastcat`.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.
    :param sample_counts: include barplot of sample counts.
    :param max_len: maximum read length.
    :param max_len: minimum read length.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    stats = read_files(stats)
    report = _maybe_new_report(report)
    report.markdown(header)

    read_length = read_length_plot(stats, min_len=min_len, max_len=max_len)
    read_qual = read_quality_plot(stats)
    plots = [[read_length, read_qual]]
    if sample_counts:
        plots.append(
            [sample_read_counts(stats, min_len=min_len, max_len=max_len)])
    report.plot(layout(plots, sizing_mode="stretch_width"))

    return report


def main(args):
    """Entry point to create a report from bcftools stats."""
    report = full_report(
        args.stats_files, report=HTMLReport(),
        min_len=args.min_len, max_len=args.max_len,
        sample_counts=args.samples)
    report.write(args.output)


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(
        "fastcat per-read report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--samples", action='store_true',
        help="Show sample counts.")
    parser.add_argument(
        "--min_len", type=int,
        help="Minimum length filter.")
    parser.add_argument(
        "--max_len", type=int,
        help="Maximum length filter.")
    parser.add_argument(
        "--output", default="fastcat_report.html",
        help="Output HTML file.")
    parser.add_argument(
        "stats_files", nargs='+',
        help=(
            "One or more files containing the per-read output from "
            "`fastcat`."))
    return parser
