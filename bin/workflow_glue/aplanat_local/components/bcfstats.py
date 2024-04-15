"""Report components for displaying information from bcftools stats."""

import argparse

from bokeh.models import LinearColorMapper
from bokeh.palettes import Blues9
from bokeh.plotting import figure
import pandas as pd

from aplanat import hist
from aplanat.parsers.bcfstats import parse_bcftools_stats_multi
from aplanat.report import _maybe_new_report, HTMLReport
from aplanat.util import Colors


def variant_counts_table(
        bcf_stats, samples_as_columns=False,
        header="**Variant counts:**", report=None):
    """Create a report section contains variant counts.

    :param bcf_stats: one or more outputs from `bcftools stats`.
    :param sample_as_columns: transpose table to put data for each
        sample into a column.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    report = _maybe_new_report(report)
    report.markdown(header)

    df = bcf_stats['SN'].drop(columns='samples')
    index = False
    if samples_as_columns:
        df = df.set_index('sample').transpose()
        index = True
    report.table(df, index=index)
    return report


def trans_counts(
        bcf_stats, samples_as_columns=False,
        header="**Transitions and tranversions:**", report=None):
    """Create a report section with transition and transversion counts.

    :param bcf_stats: one or more outputs from `bcftools stats`.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    report = _maybe_new_report(report)
    report.markdown(header)

    df = bcf_stats['TSTV']
    index = False
    if samples_as_columns:
        df = df.set_index('sample').transpose()
        index = True
    report.table(df, index=index)
    return report


_sub_header = """
**Substitution types**

Base substitutions aggregated across all samples (symmetrised by pairing)
"""


def sub_matrix(bcf_stats, header=_sub_header, report=None):
    """Create a report section with a base substitution matrix.

    :param bcf_stats: one or more outputs from `bcftools stats`.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    report = _maybe_new_report(report)
    report.markdown(header)

    sim_sub = {
        'G>A': 'C>T', 'G>C': 'C>G', 'G>T': 'C>A',
        'T>A': 'A>T', 'T>C': 'A>G', 'T>G': 'A>C'}

    def canon_sub(sub):
        b1 = sub[0]
        if b1 not in {'A', 'C'}:
            return canon_sub(sim_sub[sub])
        else:
            return b1, sub[2]

    df = bcf_stats['ST']
    df['canon_sub'] = df['type'].apply(canon_sub)
    df['original'] = df['canon_sub'].apply(lambda x: x[0])
    df['substitution'] = df['canon_sub'].apply(lambda x: x[1])
    df['count'] = df['count'].astype(int)
    df = df[['original', 'substitution', 'count']] \
        .groupby(['original', 'substitution']) \
        .agg(count=pd.NamedAgg(column='count', aggfunc='sum')) \
        .reset_index()

    colors = Blues9[::-1]
    mapper = LinearColorMapper(
        palette=colors, low=min(df['count']), high=max(df['count']))
    p = figure(
        y_range=['C', 'A'], x_range=['A', 'C', 'G', 'T'],
        x_axis_location="above",
        x_axis_label='alternative base',
        y_axis_label='reference base',
        tools="save", toolbar_location='below',
        output_backend="webgl",
        height=225, width=300,
        tooltips=[('sub', '@original>@substitution'), ('count', '@count')])
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.rect(
        source=df, y="original", x="substitution", width=1, height=1,
        fill_color={'field': 'count', 'transform': mapper},
        line_color=None)
    report.plot(p)
    return report


_indel_length_header = """
**Indel lengths**

Insertion and deletion lengths aggregated across all samples.
"""


def indel_lengths(
        bcf_stats, header=_indel_length_header, report=None,
        color=Colors.light_cornflower_blue):
    """Create a report section containing an indel length chart.

    :param bcf_stats: one or more outputs from `bcftools stats`.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.
    :param color: color of bars in chart.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    report = _maybe_new_report(report)
    report.markdown(header)
    try:
        df = bcf_stats['IDD']
    except KeyError:
        # If there are no indels, bcftools doesn't contain the table
        report.markdown("*No indels to report.*")
    else:
        df['nlength'] = df['length (deletions negative)'].astype(int)
        df['count'] = df['number of sites'].astype(int)
        # pad just to pull out axes by a minimum
        pad = pd.DataFrame({'nlength': [-10, +10], 'count': [0, 0]})
        counts = df.groupby('nlength') \
            .agg(count=pd.NamedAgg(column='count', aggfunc='sum')) \
            .reset_index().append(pad)
        plot = hist.histogram(
            [counts['nlength']], weights=[counts['count']],
            colors=[color], binwidth=1,
            title='Insertion and deletion variant lengths',
            x_axis_label='Length / bases (deletions negative)',
            y_axis_label='Count')
        report.plot(plot)
    return report


_full_report_header = """
### Variant call summaries

The following tables and figures are derived from
the output of `bcftools stats`.
"""


def full_report(bcf_stats, header=_full_report_header, report=None):
    """Create a report section from the output of bcftools stats.

    :param bcf_stats: one or more outputs from `bcftools stats`.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    if not isinstance(bcf_stats, (list, tuple)):
        bcf_stats = [bcf_stats]
    bcf_stats = parse_bcftools_stats_multi(bcf_stats)

    report = _maybe_new_report(report)
    report.markdown(header)
    report = variant_counts_table(bcf_stats, report=report)
    report = trans_counts(bcf_stats, report=report)
    report = sub_matrix(bcf_stats, report=report)
    report = indel_lengths(bcf_stats, report=report)
    return report


def main(args):
    """Entry point to create a report from bcftools stats."""
    report = full_report(args.bcfstats_files, report=HTMLReport())
    report.write(args.output)


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(
        'bcftools stats report',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--output", default="bcf_stats_report.html",
        help="Output HTML file.")
    parser.add_argument(
        "bcfstats_files", nargs='+',
        help="One or more files containing output from `bcftools stats.`")
    return parser
