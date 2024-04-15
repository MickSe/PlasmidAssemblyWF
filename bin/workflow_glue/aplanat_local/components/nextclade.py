"""Nextclade reporting section."""

import argparse

from jinja2 import Template
import pkg_resources

from aplanat.report import HTMLReport, HTMLSection


class NextClade(HTMLSection):
    """A nextclade report component."""

    def __init__(self, json, add_title=True):
        """Initialize a report component.

        :param json: json output file from nextclade CLI.
        """
        super().__init__()
        template = Template(
            """\
            <div>
            <nxt-table>
                <script defer="">
                const data = {{ data }}
                var nxt = document.querySelector('nxt-table')
                nxt.data = data.results
                </script>
            </nxt-table>
            </div>
            """  # noqa
        )
        if add_title:
            self.markdown('''
### NextClade analysis
The following view is produced by the
[nextclade](https://clades.nextstrain.org/) software.
''')
        with open(json, encoding='utf8') as fh:
            self._add_item(
                template.render(data=fh.read()))
        script = pkg_resources.resource_filename(
            'aplanat', 'data/nextclade.html')
        with open(script, encoding='utf8') as fh:
            self._script = fh.read()

    def _extra_components(self):
        """Return script and div for report."""
        # out main content div will be taken care of by the parent class
        return [self._script], dict()


def main(args):
    """Entry point to create a report from nextclade."""
    report = HTMLReport()
    report.add_section(section=NextClade(args.json))
    report.write(args.output)


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(
        'nextclade report',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--output", default="nextclade_report.html",
        help="Output HTML file.")
    parser.add_argument(
        "json",
        help="JSON output file from nextclade CLI.")
    return parser
