#!/usr/bin/env python

"""Report component for displaying information from wf-clone-validation."""
import argparse
import json
import os
import sys
import time
import re

# import the local version of the aplanat module -> done in main() and fastcat_report_tab(), as the needed path is obtained through args

from .util import wf_parser  # noqa: ABS101

from bokeh.layouts import layout
from bokeh.models import Panel, Tabs
import pandas as pd
from plannotate.bokeh_plot import get_bokeh



def tidyup_status_file(status_sheet, annotations):
    """Tidy up the sample status file."""
    sample_status = pd.read_csv(status_sheet[0], header=None)
    unique_samples = sample_status[0].unique()
    pass_fail_dic = {}
    for sample in unique_samples:
        pass_fail_dic[sample] = 'Completed successfully'
    filter_pass = sample_status[sample_status[1] != 'Completed successfully']
    failures = dict(zip(filter_pass[0], filter_pass[1]))
    completed_annotations = list(annotations)
    success = sample_status[sample_status[1] == 'Completed successfully']
    no_annotations = success[~success[0].isin(completed_annotations)]
    for sample in list(no_annotations[0]):
        failures[sample] = 'Completed but no annotations found in the database'
    all_sample_names = unique_samples.tolist()
    all_sample_names.sort()
    passed_list = unique_samples.tolist()
    for k, v in failures.items():
        pass_fail_dic[k] = v
        if v != 'Completed but failed to reconcile':
            passed_list.remove(k)
    passed_list.sort()
    status_df = pd.DataFrame(
        pass_fail_dic.items(), columns=['Sample', 'pass/failed reason'])
    sort_df = status_df['Sample'].astype(str).argsort()
    status_df = status_df.iloc[sort_df]
    return (status_df, passed_list, all_sample_names, pass_fail_dic)


def read_files(summaries, sep='\t'):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep=sep))
    return pd.concat(dfs)


def fastcat_report_tab(file_name, tab_name, aplanat_rep):

    sys.path.insert(0, aplanat_rep)
    # import the local version of the aplanat module
    import aplanat_local as aplanat_local
    from aplanat_local import bars, report
    from aplanat_local.components import fastcat
    from aplanat_local.components import simple as scomponents
    import aplanat_local.graphics
    from aplanat_local.util import Colors

    """Read fastcat dataframe and create a tab with qual and len plots."""
    df = pd.read_csv(file_name, sep='\t')
    depth = len(df.index)
    min_length = df["read_length"].min()
    max_length = df["read_length"].max()
    mean_length = df["read_length"].mean()
    mean_quality = df["mean_quality"].mean()
    lengthplot = fastcat.read_length_plot(
        df,
        min_len=min_length,
        max_len=max_length)
    qstatplot = fastcat.read_quality_plot(df)
    exec_summary = aplanat_local.graphics.InfoGraphItems()
    exec_summary.append(
        'No. reads',
        str(depth),
        "angle-up", '')
    exec_summary.append( 
        'Mean read length',
        str(round(mean_length)),
        "align-center", '')
    exec_summary.append(
        'Mean read quality',
        str(round(mean_quality)),
        "signal", '')
    exec_plot = aplanat_local.graphics.infographic(
        exec_summary.values(), ncols=4)
    tab = Panel(child=layout(
        [[exec_plot], [lengthplot], [qstatplot]],
        aspect_ratio="auto",
        sizing_mode='stretch_width'),
        title=tab_name)
    return tab


def create_fastcat_dic(sample_names, raw, hostfilt, downsampled):
    """Create dictionary using sample names and fastcat files available."""
    per_sample_dic = {}
    #lists = {'raw': raw, 'hostfilt': hostfilt, 'downsampled': downsampled}
    lists = {'raw reads': raw, 'hostfilt reads': hostfilt, 'downsampled reads': downsampled}
    for sample in sample_names:
        new_dic = {}
        item_search = '/' + sample + '.'
        for list_name, fc_list in lists.items():
            #  find index of item that contains the sample name as a substring
            item_index = [i for i, s in enumerate(fc_list) if item_search in s]
            if item_index:
                indice = item_index[0]
                new_dic[list_name] = fc_list[indice]
            else:
                pass
        per_sample_dic[sample] = new_dic
    return per_sample_dic


def main(args):

    aplanat_rep = str(args.aplanat)
    sys.path.insert(0, aplanat_rep)

    # import the local version of the aplanat module
    import aplanat_local as aplanat_local
    from aplanat_local import bars, report
    from aplanat_local.components import fastcat
    from aplanat_local.components import simple as scomponents
    import aplanat_local.graphics
    from aplanat_local.util import Colors


    # We defer this until processing through all the samples in the loop below
    plannotate_annotations = json.load(open(args.plannotate_json)).keys()
    lengths = json.load(open(args.lengths))
    pass_fail = tidyup_status_file(args.status, plannotate_annotations)

    passed_samples = pass_fail[1]
    sample_names = pass_fail[2]

    # if variant calling was performed, full paths to vcf files are extracted 
    # and stored in a list for both sniffles and medaka
    vcf_paths_sniffles = []
    if args.VCFsniffles:
        for vcf_string in args.VCFsniffles:
            
            vcf_string = vcf_string.strip('[]')
            paths = vcf_string.split(',')
            
            for path in paths:
                vcf_paths_sniffles.append(path.strip())

    vcf_paths_medaka = []
    if args.VCFmedaka:
        for vcf_string in args.VCFmedaka:
            
            vcf_string = vcf_string.strip('[]')
            paths = vcf_string.split(',')
            
            for path in paths:
                vcf_paths_medaka.append(path.strip())

    # HES and ONT encoded as base64 strings and stored in annexes txt files 
    # are written in an appropriate variable
    file = open(aplanat_rep + "/logo_hes.txt", "r")
    logo_hes = file.read()
    file.close()

    file = open(aplanat_rep + "/logo_ont.txt", "r")
    logo_ont = file.read()
    file.close()


    # iterates through all samples to generate a new html report for each
    for item in sample_names:

        """Entry point to create a wf-clone-validation report."""
        report_doc = report.HTMLReport(
            "Plasmid Assembly Report : " + item,
            "Draft report for clone validation workflow",
            style='hes')
        

        # Report passed analysis
        if item in passed_samples:

            # Per sample details
            section = report_doc.add_section()

            section.markdown('''Status :
                             <font color='green'>*{}*
                             </font>'''.format(pass_fail[3][item]))

            section.markdown("""
            ## Assembly
            A [pLannotate plot](https://github.com/barricklab/pLannotate) is displayed
            as well as a feature table providing description of the annotated sequence.

            Unfilled features on the plannotate plots are incomplete features; the sequence
            match in the plasmid covers less than 95% of the full length of the feature in
            the database. These elements may be leftover fragments from earlier cloning
            steps used to create a plasmid. If they include only a small fraction of the
            feature, they likely do not still have the annotated function. However, even
            small feature fragments may affect plasmid function if they result in cryptic
            gene expression or are inadvertently combined with other elements during later
            cloning steps.
                                

            """)

            host_ref_stats = args.host_filter_stats
            downsampled_stats = args.downsampled_stats
            initial_stats = args.per_barcode_stats
        
            if ('host_filter_stats/OPTIONAL_FILE' in host_ref_stats):
                host_filt = host_ref_stats.remove('host_filter_stats/OPTIONAL_FILE')
                if host_filt is None:
                    host_filt = []
            else:
                host_filt = host_ref_stats
        
            if ('host_filter_stats/OPTIONAL_FILE' in downsampled_stats):
                summary_stats = downsampled_stats.remove(
                    'downsampled_stats/OPTIONAL_FILE')
                if summary_stats is None:
                    summary_stats = []
            else:
                summary_stats = downsampled_stats

            fast_cat_dic = create_fastcat_dic(
                sample_names, initial_stats, host_filt, summary_stats)
            plannotate = json.load(open(args.plannotate_json))
            sample_stats = []

            # stats graphs and plannotate where appropriate
            tup_dic = plannotate[item]
            annotations = pd.read_json(tup_dic['annotations'])
            plot = get_bokeh(pd.read_json(tup_dic['plot']), False)

            fast_cat_tabs = fast_cat_dic[item]
            alltabs = []

            for key, value in fast_cat_tabs.items():
                alltabs.append(fastcat_report_tab(value, key,aplanat_rep))
            cover_panel = Tabs(tabs=alltabs)
            plasmidplot = [plot]
            # sizing_mode="stretch_width" was added 
            section.plot(layout([[plasmidplot]], sizing_mode="stretch_width"))
            
            section.markdown('''<div style="text-align: left; margin-top: 60px;">
                                </div>''')
            section.markdown('''
                                    ## Annotation table
                             ''')

            section.table(annotations.drop(columns=['Plasmid length']),
                        index=False, key="table"+str(item))


            section = report_doc.add_section()
            section.markdown('''<div style="text-align: left; margin-top: 60px;">
                                </div>''')
            section.markdown('''
                                    ## Quality control
                             ''')

            section.plot(cover_panel)

            stats_df = pd.DataFrame(sample_stats, columns=["Sample", "Length"])
            merged_inner = pd.merge(pass_fail[0], stats_df)

            merged_inner.to_csv('sample_status.txt', index=False)



            # check if sample name is found in args.VCFsniffles if yes, use that path
            for path in vcf_paths_sniffles:

                if item in path:
                    VCF = path

                    section = report_doc.add_section()
                    section.markdown('''<div style="text-align: left; margin-top: 60px;">
                                    </div>''')
                    section.markdown('''
                                    ## Variant Calling
                                     
                                    The tables below contain all detected single nucleotide polymorphisms, short and structural variants. The columns show respectively the position of the variation on the reference sequence, the reference base or bases at that position as well as the alternative allele found, the type of variant (INS = insertion, DEL = deletion, DUP = duplication, INV = inversion & BND = breakend, translocation), its lengths and finally the number of reads supporting the structural variation. Sequences longer than 50 nucleotides are not displayed here for readability reasons. For more informations, consider the corresponding vcf file.

                                    ### Structural variants:
                                    ''')
                    
                    # Note: 'try' is used to avoid an error if VCF file is empty
                    try:
                        Variants = pd.read_csv(VCF, sep=' ', header=None)
                    except:
                        section.markdown('''
                                            No structural variants where found by comparing with the given reference.
                                    ''')
                    else:

                        Variants = pd.read_csv(VCF, sep=' ', header=None)

                        Variants = Variants.drop(columns=Variants.columns[[0, 2, 6, 9]])
                        Variants.columns = ['Position', 'Reference', 'Alternative', 'Quality', 'Type', 'Length']
                        Variants[['Reference', 'Alternative']] = Variants[['Reference', 'Alternative']].applymap(lambda x: 'sequence longer than 50 nucleotides, see corresponding vcf file' if isinstance(x, str) and len(x) > 50 else x)

                        section.table(Variants,
                                    index=False, key="vcf_table"+str(item))


            # check if sample name is found in args.VCFmedaka if yes, use that path
            for path in vcf_paths_medaka:

                if item in path:
                    VCF = path

                    section = report_doc.add_section()
                    section.markdown('''<div style="text-align: left; margin-top: 60px;">
                                    </div>''')
                    section.markdown('''
                                    ### Short variants:
                                    ''')
                    
                    # Note: 'try' is used to avoid an error if VCF file is empty
                    try:
                        Variants = pd.read_csv(VCF, sep=' ', header=None)
                    except:
                        section.markdown('''
                                            No short variants where found by comparing with the given reference
                                    ''')
                    else:

                        Variants = pd.read_csv(VCF, sep=' ', header=None)
                        Variants = Variants.drop(columns=Variants.columns[[0, 2, 6]])
                        Variants.columns = ['Position', 'Reference', 'Alternative', 'Quality']
                        Variants[['Reference', 'Alternative']] = Variants[['Reference', 'Alternative']].applymap(lambda x: 'sequence longer than 50 nucleotides, see corresponding vcf file' if isinstance(x, str) and len(x) > 50 else x)

                        section.table(Variants,
                                    index=False, key="vcf_table"+str(item))


        # Report failed analysis (in case Flye failed)
        else:
            section = report_doc.add_section()
            section.markdown('''Status :
                             <font color='red'>*{}*
                             </font>'''.format(pass_fail[3][item]))
            section.markdown('''<div style="text-align: left; margin-top: 60px;">
                                </div>''')
            section.markdown('''
                                    ## Quality control
                             ''')
            
            fast_cat_dic = create_fastcat_dic(
                sample_names, initial_stats, host_filt, summary_stats)
            sample_stats = []
            
            fast_cat_tabs = fast_cat_dic[item]
            alltabs = []
            
            for key, value in fast_cat_tabs.items():
                alltabs.append(fastcat_report_tab(value, key, aplanat_rep))
            cover_panel = Tabs(tabs=alltabs)
            section.plot(cover_panel)
            stats_df = pd.DataFrame(sample_stats, columns=["Sample", "Length"])
            merged_inner = pd.merge(pass_fail[0], stats_df)
            merged_inner.to_csv('sample_status.txt', index=False)

        from datetime import date
        YYY_MM_DD = str(date.today())
        section = report_doc.add_section()
        section.markdown('''<div style="text-align: left; margin-top: 60px;">
                                </div>''')
        section.markdown('''
                                #### About this report:
                         ''')
        section.markdown(f"""
        This report was generated on {YYY_MM_DD} using an extended version of the epi2me-labs/wf-clone-validation pipeline ([nf-core framework](https://doi.org/10.1038/s41587-020-0439-x)).
                         

        <div style="text-align: center; margin-top: 100px;">
            <img src="{logo_hes}" style="display: inline-block; margin-right: 20px; max-width: 20%;">
                         
            <img src="{logo_ont}" style="display: inline-block; max-width: 20%;">
        </div>           
        """)


        # write report in html format
        # report_doc.write(args.report_name + item + '.html')
        report_doc.write(item + '.report.html')



def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")

    parser.add_argument(
        "--downsampled_stats",
        nargs='+',
        required=True
    )
    parser.add_argument(
        "--revision", default='unknown',
        help="revision number")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit number")
    parser.add_argument(
        "--status", nargs='+',
        help="status")
    parser.add_argument(
        "--per_barcode_stats", nargs='+',
        help="fastcat stats file for each sample before filtering")
    parser.add_argument(
        "--host_filter_stats", nargs='+',
        help="fastcat stats file after host filtering")
    parser.add_argument(
        "--params", default=None,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--versions",
        help="directory contained CSVs containing name,version.")
    parser.add_argument(
        "--plannotate_json",
        help="Plannotate Json.")
    parser.add_argument(
        "--inserts_json",
        help="inserts Json.")
    parser.add_argument(
        "--report_name",
        help="report name")
    parser.add_argument(
        "--lengths",
        help="report name")
    parser.add_argument(
        "--aplanat", default=None,
        help="path to aplanat library location")
    parser.add_argument('--VCFsniffles', nargs='+', help='Paths to VCF files separated by semicolons', default=None)
    parser.add_argument('--VCFmedaka', nargs='+', help='Paths to VCF files separated by semicolons', default=None)
    return parser


if __name__ == "__main__":
    args = argparse().parse_args()
    main(args)
