#!/usr/bin/env python

"""Report component for displaying information from wf-clone-validation."""
import argparse
import json
import os

from aplanat import json_item
import numpy as np
import pandas as pd
from plannotate.annotate import annotate
from plannotate.bokeh_plot import get_bokeh
import pysam
from .util import wf_parser  # noqa: ABS101


def run_plannotate(fasta, linear=False):
    """Run annotate and create Bokeh plot."""
    with pysam.FastxFile(fasta) as fh:
        seq = next(fh).sequence
    df = annotate(
        seq, is_detailed=True, linear=linear, yaml_file="plannotate.yaml")
    plot = get_bokeh(df, linear=linear)
    plot.xgrid.grid_line_color = None
    plot.ygrid.grid_line_color = None
    old_df = df.copy(deep=True)
    clean_df = clean_results(df)
    old_df.reset_index(drop=True, inplace=True)
    return plot, clean_df, old_df


def clean_results(df):
    """Clean-up annotation dataframe for display."""
    rename = {
        'Feature': 'Feature',
        'db': 'Database',
        'pident': 'Identity', 'abs percmatch': 'Match Length',
        'Description': 'Description',
        'qstart': 'Start Location', 'qend': 'End Location', 'length': 'Length',
        'sframe': 'Strand', 'db': 'Database',
        'qlen': 'qlen'}
    numeric_columns = ['Identity', 'Match Length']
    display_columns = [
        'Feature',
        'Database',
        'Identity', 'Match Length',
        'Description',
        'Start Location', 'End Location', 'Length',
        'Strand', 'qlen']
    df = df.rename(columns=rename)[display_columns]
    df['Plasmid length'] = df.iloc[0]['qlen']
    df = df.drop(columns='qlen')
    df[numeric_columns] = np.round(df[numeric_columns], 1).astype(str) + "%"
    df.loc[df['Database'] == "infernal", 'Identity'] = "-"
    df.loc[df['Database'] == "infernal", 'Match Length'] = "-"
    df = df.set_index("Feature", drop=True).reset_index()
    return df


#####################################################
#####################################################
# added function:
# generate a .gbk file from .fasta and .bed format

def fa2gbk(sample_file, bed, item):
    
    #preparing input files
    input_fasta_file = open("{}".format(sample_file))
    input_bed_file = bed
  
    #fasta file
    input_fasta = input_fasta_file.read()
    header = input_fasta.splitlines()[0]
    fasta_seq = input_fasta.splitlines()[1].lower()
    fasta = []
 
    fasta.append(fasta_seq)
    fasta = fasta[0]
    result0 = []
    for i in range(0, len(fasta), 10):
        result0.append(fasta[i:i+10])

    #bed file
    bed_lines = input_bed_file.to_numpy().tolist()
    result_bed = []
    for b in bed_lines:
        result_bed.append(b)

    #writing GBK file into a variable
    result1 = []
    result1.append("LOCUS       {0}               {1} bp ds-DNA     linear".format(header[1:], len(fasta)))
    result1.append("DEFINITION  {}".format(item))
    result1.append("ACCESSION   .")
    result1.append("VERSION     .")
    result1.append("KEYWORDS    {}".format(header[1:]))
    result1.append("SOURCE      natural DNA sequence")
    result1.append("  ORGANISM  {}".format('.'))
    result1.append("REFERENCE   1  (bases 1 to {})".format(len(fasta)))
    result1.append("  AUTHORS   {}".format('Clone validation workflow'))
    result1.append("  TITLE     {} locus GenBank formatted file".format(header[1:]))
    result1.append("  JOURNAL   Unpublished")
    result1.append("FEATURES             Location/Qualifiers")
    result1.append("     source          1..{}".format(len(fasta)))
    result1.append("                     /organism=\"{}\"".format('.'))
    result1.append("                     /mol_type=\"{}\"".format('genomic DNA'))
    for i in result_bed:
        # start position +1 compared to the bed file... is that normal ???
        result1.append("     misc_feature    {0}..{1}".format(str(int(i[1])+1), i[2]))
        result1.append("                     /label=\"{}\"".format(i[3]))

    result1.append("ORIGIN")
    for h, i in zip(range(1, len(fasta), 60), range(0, len(result0), 6)):
            result1.append("{:>9}".format(str(h)) + " " + " ".join(result0[i:i+6]))
    result1.append("//")
    #printing GBK variable as a single formatted string-file
    print(*result1, sep='\n')

    genbank_file = open(header[1:] + ".annotations.gbk", "w")
    genbank_file.writelines('\n'.join([*result1]))
    genbank_file.close()

#############################################################
#############################################################

def bed_file(item, df):
    """Bed format for annotations."""
    display_columns = [
        'Start Location', 'End Location',
        'Feature',
        'Strand']
    df = df[display_columns]
    df["Strand"] = df["Strand"].apply(pd.to_numeric)
    df.loc[df['Strand'] == 0, 'Strand'] = "-"
    df.loc[df['Strand'] == 1, 'Strand'] = "+"
    df.insert(0, 'Name', value=str(item))
    df.to_csv(
        str(item)+'.annotations.bed', sep="\t", header=False, index=False)
    return df


def per_assembly(sample_file, item):
    """Run plannotate for a sample.

    :param database:
    :param sample_file:
    :param item: the sample
    """
    plot, annotations, clean_df = run_plannotate(sample_file)
    bed = bed_file(item, annotations)
    fa2gbk(sample_file, bed, item)
    with pysam.FastxFile(sample_file) as fh:
        seq_len = len(next(fh).sequence)
    tup = {
        'sample_name': item,
        'plot': plot,
        'annotations': annotations,
        'seq_len': seq_len}
    return tup, clean_df


def output_feature_table(data):
    """Build feature table text file or if no data output empty file."""
    if data:
        df = data[0]['annotations']
        sample_column = data[0]['sample_name']
        df.insert(0, 'Sample_name', sample_column)
        df.to_csv('feature_table.txt', mode='a', header=True, index=False)
        for sample in data[1:]:
            df = sample['annotations']
            sample_column = sample['sample_name']
            df.insert(0, 'Sample_name', sample_column)
            df.to_csv('feature_table.txt', mode='a', header=False, index=False)
    else:
        # If no samples passed create empty feature_table file
        feature_file = open("feature_table.txt", "w")
        feature_file.write(
"""Sample_name,Feature,Uniprot ID,Database,Identity,Match Length,Description,Start Location,End Location,Length,Strand,Plasmid length""") # noqa
        feature_file.close()


def make_yaml(database):
    """Create a yaml file for plannotate."""
    plannotate_yaml = """
Rfam:
  details:
    compressed: false
    default_type: ncRNA
    location: None
  location: {0}
  method: infernal
  priority: 3
  version: release 14.5
fpbase:
  details:
    compressed: false
    default_type: CDS
    location: Default
  location: {0}
  method: diamond
  parameters:
  - -k 0
  - --min-orf 1
  - --matrix BLOSUM90
  - --gapopen 10
  - --gapextend 1
  - --algo ctg
  - --id 75
  priority: 1
  version: downloaded 2020-09-02
snapgene:
  details:
    compressed: false
    default_type: None
    location: Default
  location: {0}
  method: blastn
  parameters:
  - -perc_identity 95
  - -max_target_seqs 20000
  - -culling_limit 25
  - -word_size 12
  priority: 1
  version: Downloaded 2021-07-23
swissprot:
  details:
    compressed: true
    default_type: CDS
    location: Default
  location: {0}
  method: diamond
  parameters:
  - -k 0
  - --min-orf 1
  - --matrix BLOSUM90
  - --gapopen 10
  - --gapextend 1
  - --algo ctg
  - --id 50
  priority: 2
  version: Release 2021_03
        """.format(database)

    with open("plannotate.yaml", "w") as text_file:
        text_file.write(plannotate_yaml)


def attempt_annotation(sample_file, name):
    """Create annotation dictionary for report and EPI2ME."""
    tup_dic, clean_df = per_assembly(sample_file, name)
    plasmid_len = tup_dic['annotations']['Plasmid length'][0]
    feature_dic = tup_dic['annotations'].drop(
                    ['Plasmid length'], axis=1)
    features = feature_dic.to_dict('records')
    output_json = json_item(tup_dic['plot'])
    plannotate_dic = {
        "reflen": float(plasmid_len),
        "features": features,
        "plot": output_json['doc']}
    report = {}
    report['sample_name'] = name
    report['plot'] = clean_df.to_json()
    report['annotations'] = tup_dic['annotations'].to_json()
    report['seq_len'] = tup_dic['seq_len']
    return tup_dic, report, plannotate_dic


def main(args):
    """Entry point to create a wf-clone-validation report."""
    final_samples = []
    report_dic = {}
    plannotate_collection = {}
    make_yaml(args.database)
    json_file = open("plannotate.json", "a")
    if args.sequences:
        for filename in os.listdir(args.sequences):
            name = str(filename).split('.final')[0]
            file = os.path.join(args.sequences, filename)
            try:
                tup_dic, report, plannotate_dic = attempt_annotation(
                    file, name)
                final_samples.append(tup_dic)
                plannotate_collection[name] = plannotate_dic
                report_dic[name] = report
            except KeyError:
                with pysam.FastxFile(file) as fh:
                    seq_len = len(next(fh).sequence)
                plannotate_dic = {
                    "reflen": float(seq_len),
                    "features": [],
                    "plot": {}}
                plannotate_collection[name] = plannotate_dic
                continue

    # outputs for epi2me
    output_feature_table(final_samples)
    json_object = json.dumps(plannotate_collection, indent=4)
    json_file.write(json_object)
    json_file.close()

    # outputs for report
    json_file = open("plannotate_report.json", "a")
    json_object = json.dumps(report_dic)
    json_file.write(json_object)
    json_file.close()


def argparser():
    """Argument parser for entry point."""
    parser = wf_parser("run_plannotate")
    parser.add_argument(
        "--database", default='unknown',
        help="database to use, directory containing BLAST et. al. files.")
    parser.add_argument(
        "--sequences",
        help="sequences in directory to run plannotate on.",
        required=False)
    return parser


if __name__ == "__main__":
    args = argparse().parse_args()
    main(args)
