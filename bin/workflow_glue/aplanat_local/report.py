"""Report building from multiple items."""

from collections import OrderedDict
import textwrap
import uuid

from bokeh.embed import components
from bokeh.model import Model
from bokeh.models.sources import ColumnDataSource
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.resources import INLINE
from jinja2 import Template
import markdown
import pkg_resources

JS_RESOURCES = [
    'simple-datatables_latest.js']


def _maybe_new_report(section, require_keys=False):
    """Create a new report section, or returns input."""
    if section is None:
        section = HTMLSection(require_keys=require_keys)
    else:
        if not isinstance(section, HTMLSection):
            raise TypeError("`section` should be an `HTMLSection`.")
    return section


class HTMLSection(OrderedDict):
    """A section of a report.

    The class can be extended o provide novel report components by
    overriding `._extra_components()`.
    """

    def __init__(self, require_keys=False):
        """Initialize the report item collection.

        :param title: report title.
        :param lead: report strapline, shown below title.
        :param require_keys: require keys when adding items.
        """
        self.require_keys = require_keys
        self.plots = list()
        self.md = markdown.Markdown()

    def _add_item(self, item, key=None):
        """Add an item to the report.

        :param item: item to add.
        :param key: unique key for item.
        """
        if key is None:
            if self.require_keys:
                raise ValueError('A key is required.')
            else:
                key = str(uuid.uuid4())
        self[key] = item

    def placeholder(self, key):
        """Add a placeholder to be filled in later."""
        self._add_item(None, key=key)

    def plot(self, plot, key=None):
        """Add a plot to the report.

        :param plot: bokeh plot instance.
        :param key: unique key for item.
        """
        self._add_item(plot, key=key)

    def table(
            self, data_frame, index=False, key=None, searchable=True,
            paging=True, sortable=True, **kwargs):
        """Add a pandas dataframe to the report.

        :param df: pandas dataframe instance.
        :param index: include dataframe index in output.
        :param key: unique key for item.
        :param shrink: shrink wrap the table to avoid whitespace.
        :param kwargs:
                       https://github.com/fiduswriter/Simple-DataTables/wiki/API
        """
        table = Table(
            data_frame, index=index, searchable=searchable,
            paging=paging, sortable=sortable, **kwargs)
        self._add_item(table, key=key)

    def markdown(self, text, key=None):
        """Add markdown formatted text to the report.

        :param text: markdown formatted text. The text will be dedented before
            use so it it safe to use triple-quoted strings indented to match
            code indentation.
        :param key: unique key for item.
        """
        if text is None or text == '':
            return
        text = textwrap.dedent(text)
        html = self.md.convert(text)
        self.md.reset()
        self._add_item(html, key=key)

    def alert(self, title, text, level, key=None):
        """Add an alert to the report.

        :param title: plain text - emphasis used on report
        :param text: plain text - the text of the warning
        :param level: the level of the warning; danger, warning, success, \
            or info
        :param key: unique key for item.
        """
        levels = ['danger', 'warning', 'success', 'info']
        if level not in levels:
            raise ValueError(
                'The alert level must be one of danger, warning, success, \
                    or info')

        if text is None or text == '':
            return

        html = '<div class="alert alert-'+level + \
            '"><p><strong>'+title+'</strong></p>'+text+'</div>'
        self._add_item(html, key=key)

    def _plot_components(self):
        """Return html script and div tags for bokeh plots."""
        # handle bokeh plots
        scripts = list()
        divs = dict()
        plots = {k: v for k, v in self.items() if isinstance(v, Model)}
        if len(plots) > 0:
            scripts, divs = components(plots)
            scripts = [scripts]
        return scripts, divs

    def _table_components(self):
        """Return html script and div for tables."""
        # TODO: reimplement `.table()` to not use bokeh?
        scripts = list()
        divs = dict()
        tables = {k: v for k, v in self.items() if isinstance(v, Table)}
        if len(tables) > 0:
            for k, v in tables.items():
                divs[k] = v.div
        return scripts, divs

    def _extra_components(self):
        """Return additional script and divs."""
        # hook for derived classes:
        # return: list of scripts, dict of divs keyed on self.keys
        return list(), dict()

    def components(self):
        """Fetch script and div tags for report."""
        scripts = list()
        specials = dict()
        # retrieve bokeh plot divs
        plot_scripts, plot_divs = self._plot_components()
        scripts.extend(plot_scripts)
        specials.update(plot_divs)
        # TODO: add special handing for new tables?
        table_scripts, table_divs = self._table_components()
        scripts.extend(table_scripts)
        specials.update(table_divs)
        # retrieve any extra specials
        extra_scripts, extra_divs = self._extra_components()
        scripts.extend(extra_scripts)
        specials.update(extra_divs)

        # put things together in order
        divs = list()
        for k in self.keys():
            try:
                divs.append(specials[k])
            except KeyError:
                if self[k] is None:
                    raise ValueError(
                        "Placeholder `{}` was not assigned "
                        "a value.".format(k)
                    )
                # fall back to putting items in directly
                divs.append(self[k])
        return scripts, divs


class HTMLReport(HTMLSection):
    """Generate HTML Report from a series of bokeh figures.

    Items added to the report take an optional key argument, adding items
    with the same key allows an update in place whilst maintaining the order
    in which items were added. Items can be grouped into sections for easier
    out of order addition.
    """

    def __init__(self, title="", lead="", require_keys=False, style='hes'):
        """Initialize the report item collection.

        :param title: report title.
        :param lead: report strapline, shown below title.
        :param require_keys: require keys when adding items.
        """
        super().__init__(require_keys=require_keys)
        self.title = title
        self.lead = lead
        self.sections = OrderedDict()
        self.sections['main'] = self
        self.style = style


        self.CSS_RESOURCES = dict(
                bootstrap='bootstrap.min.css',
                datatables='simple-datatables_latest.css',
                ont='custom-ont.css',
                hes='custom-hes.css',       # ADDED AS A TEST TO MODIFY REPORT ASPECT (document is in repo_clone/bin/workflow_glue/aplanat_local/data/custom-hes.css)
                ond='custom-ond.css',
                epi2me='custom-epi2me.css'
        )

        template = pkg_resources.resource_filename(
            __package__, 'data/report_template.html')
        with open(template, 'r', encoding="UTF-8") as fh:
            template = fh.read()
        self.template = Template(template)

    def add_section(self, key=None, section=None, require_keys=False):
        """Add a section (grouping of items) to the report.

        :param key: unique key for section.
        :param section: `HTMLSection` to add rather than creating anew.

        :returns: the report section.

        """
        if key is None:
            key = str(uuid.uuid4())
        section = _maybe_new_report(section, require_keys=require_keys)
        self.sections[key] = section
        return self.sections[key]

    def render(self):
        """Generate HTML report containing figures."""
        bokeh_resources = INLINE.render()

        libs = []
        css_files = ['bootstrap', 'datatables', self.style]
        CSS_RESOURCES = [self.CSS_RESOURCES[file] for file in css_files]
        for resources, stub in (
                [CSS_RESOURCES, "<style>{}</style>"],
                [JS_RESOURCES, '<script type="text/javascript">{}</script>']):
            for res in resources:
                fn = pkg_resources.resource_filename(
                    __package__, 'data/{}'.format(res))
                with open(fn, encoding="UTF-8") as fh:
                    libs.append(stub.format(fh.read()))

        all_scripts = list()
        all_divs = list()
        for sec_name, section in self.sections.items():
            scripts, divs = section.components()
            all_scripts.extend(scripts)
            all_divs.extend(divs)

        script = '\n'.join(all_scripts)
        divs = '\n'.join(all_divs)
        return self.template.render(
            title=self.title, lead=self.lead,
            bokeh_resources=bokeh_resources, resources="\n".join(libs),
            script=script, div=divs)

    def write(self, path):
        """Write html report to file."""
        with open(path, "w", encoding='utf8') as outfile:
            outfile.write(self.render())


class WFReport(HTMLReport):
    """Report template for epi2me-labs/wf* workflows."""

    def __init__(
            self, title, workflow, commit=None, revision=None,
            require_keys=False, about=True, style='ont'):
        """Initialize the report item collection.

        :param workflow: workflow name (e.g. wf-hap-snps)
        :param title: report title.
        :param require_keys: require keys when adding items.
        """
        if commit is None:
            commit = "unknown"
        if revision is None:
            revision = "unknown"
        self.wf_commit = commit
        self.wf_revision = revision
        self.workflow = workflow
        self.about = about
        self.style = style

        lead = (
            "Results generated through the {} Nextflow workflow "
            "provided by Oxford Nanopore Technologies.".format(workflow))
        super().__init__(
            title=title, lead=lead, require_keys=require_keys, style=style)
        self.tail_key = str(uuid.uuid4())

    def render(self):
        """Generate HTML report containing figures."""
        # delete and re-add the tail (in case we are called twice,
        # and something was added).
        try:
            del self.sections[self.tail_key]
        except KeyError:
            pass

        if self.about:
            self.add_section(key=self.tail_key).markdown("""
    ### About

    This report was produced using the
    [epi2me-labs/{0}](https://github.com/epi2me-labs/{0}).  The
    workflow can be run using `nextflow epi2me-labs/{0} --help`

    **Version details** *Revision*: {1} *Git Commit*: {2}


    ***Oxford Nanopore Technologies products are not intended for use for
    health assessment or to diagnose, treat, mitigate, cure or prevent any
    disease or condition.***


    ---
    """.format(self.workflow, self.wf_revision, self.wf_commit))
        return super().render()


def bokeh_table(df, index=True, **kwargs):
    """Generate a bokeh table from a pandas dataframe."""
    columns = [TableColumn(field=x, title=x) for x in df.columns]
    if index:
        columns = [TableColumn(field="index", title="Feature")] + columns
    return DataTable(
        columns=columns, source=ColumnDataSource(df), index_position=None,
        **kwargs)


class Table():
    """A table report component."""

    # color modified (previously #0084A9)
    def __init__(self, data_frame, index, th_color='#6205A8', **kwargs):
        """Initialize table component.

        :param dataframe: dataframe to turn in to simple table.
        """
        template = Template(
            """\
            <style>
            #{{ table_id }}{
            font-family: Arial, Helvetica, sans-serif;
            border-collapse: collapse;
            word-wrap: break-word;
            white-space: pre-wrap;
            width: 100%;
            }
            #{{ table_id }} td, #{{ table_id }} th {
            border: 1px solid #ddd;
            padding: 4px;
            border-collapse: collapse;
            word-wrap: break-word;
            white-space: pre-wrap;
            }
            #{{ table_id }} tr:nth-child(even){background-color: #F3EAF9;}
            #{{ table_id }} tr:hover {background-color: #D4A0FF;}
            #{{ table_id }} th {
            padding-top: 4px;
            padding-bottom: 4px;
            text-align: left;
            background-color: """+th_color+""" ;
            color: white;
            border-collapse: collapse;
            word-wrap: break-word;
            white-space: pre-wrap;
            }
            </style>

                    <body>
                        {{ dataframe }}
                    </body>
                   <script defer type="text/javascript">
                        let {{ table_id|safe }}=new simpleDatatables.DataTable(
                            "#{{ table_id }}",{
                         {% for k, v in kwargs.items() %}
                            {{ k }}: {{ v }},
                        {% endfor %}
                        });
                    </script>

            """)
        for key, val in kwargs.items():
            if isinstance(val, bool):
                kwargs[key] = str(val).lower()
        key = str(uuid.uuid4()).replace('-', '')
        key = "a{}".format(key[1:])
        # html id must not begin with numeric

        self.div = template.render(
            dataframe=data_frame.to_html(table_id=key, index=index),
            table_id=key, kwargs=kwargs)
