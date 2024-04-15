"""Simple report components."""

import argparse
import json
import os

import pandas as pd

from aplanat.report import _maybe_new_report, HTMLReport
from aplanat.util import get_named_logger


_version_header = """
### Software versions

The table below highlights versions of key software used within the analysis.
"""


def version_table(
        versions, header=_version_header, report=None, th_color='#0084A9'):
    """Create a software version report from CSVs with information.

    :param versions: directory containing headerless CSVs with
        'software,version', or a single file.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    logger = get_named_logger('VerReport')
    report = _maybe_new_report(report)
    report.markdown(header)

    if os.path.isdir(versions):
        versions = [os.path.join(versions, x) for x in os.listdir(versions)]
    elif os.path.isfile(versions):
        versions = [versions]
    else:
        raise IOError('`versions` should be a file or directory.')

    verdata = list()
    for fname in versions:
        logger.info(f"Reading versions from file: {fname}.")
        try:
            with open(fname, 'r') as fh:
                for line in fh.readlines():
                    name, version = line.strip().split(',')
                    verdata.append((name, version))
        except Exception:
            logger.warning(f"Failed to read versions from: {fname}.")
    verdata = pd.DataFrame(verdata, columns=('Name', 'Version'))
    report.table(verdata, index=False, th_color=th_color)
    return report


_params_header = """
### Workflow parameters

The table below highlights values of the main parameters used in this analysis.
"""


def params_table(
        params, header=_params_header, report=None, th_color='#0084A9'):
    """Create a workflow parameter report from a JSON file.

    :param versions: Flat JSON file containing key, value pairs
        corresponding to workflow parameter 'key/name,value' pairs.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    report = _maybe_new_report(report)
    report.markdown(header)

    if not os.path.isfile(params):
        raise IOError('`params` should be a JSON file.')

    with open(params) as f:
        data = json.load(f)
        params = [(k, v) for k, v in data.items()]
        df_params = pd.DataFrame(params, columns=['Key', 'Value'])
        report.table(df_params, index=False, th_color=th_color)
        return report


def main(args):
    """Entry point to create demo report of functionality in this module."""
    report = HTMLReport(
        title="Simple component demo.",
        lead="A demonstration of small, simple reporting components.")
    report.add_section(section=version_table(args.versions))
    report.add_section(section=params_table(args.params))
    report.write(args.output)


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(
        "Simple report elements demo.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--versions",
        help=(
            "Headerless CSV containing 'software,version', "
            "or directory of such files."))
    parser.add_argument(
        "--params",
        help=(
            "A JSON file containing the workflow parameter "
            "key/values."))
    parser.add_argument(
        "--output", default="simple_components_report.html",
        help="Output HTML file.")
    return parser
