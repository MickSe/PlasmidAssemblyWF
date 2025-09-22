The attached workflow is a pipeline for reconstructing plasmids from fastq files (possibly barcoded), generating an assembly in fasta format, an annotation in both bed and gbk formats, a report in html as well as a variant table in both text and vcf format if a reference sequence was provided. It is a modified version of the workflow used by the EPI2MElabs software proposed by oxford nanopore technologies.

The original version of this workflow is available here:

https://github.com/epi2me-labs/wf-clone-validation/releases/tag/v0.3.1

Note that it is no longer the most recent version, and important differences have appeared in the meantime, notably the use of the aplanat library, which is no longer required for report generation in v0.4.0. Also, report aspect has been modified in v0.5.0 and annotation in gbk format is generated in v0.5.1.

The main modifications introduced to this initial workflow (v0.3.1) are presented below.


## GENERAL USAGE


The installation of nextflow and its dependencies can be tested using this command:

> wf='/path/to/the/workflow/repository/main.nf'
> nextflow run $wf --help


The workflow is used by passing input arguments as presented in the typical example below:

> wf='/path/to/the/workflow/repository/main.nf'
> nextflow run $wf \\
>         --fastq /path/to/directory/containing/input/fastq \\
>         --db\_directory /path/to/database/built/during/installation \\
>         --threads 4
 

Using the command: nextflow run $wf --help
will provide additional information on available arguments

- The path given to the --fastq parameters can contain subdirectories, which will be considered as different samples (barcodes).

- parameter --approx_size_sheet can be used to specify a different plasmid size estimation for each barcode and must give the path to a text file.

- parameter --ref_seq will take as input the path to a reference fasta sequence and allow the workflow to perform variant calling and return a vcf table. If this argument is used with several samples (barcodes), the same reference will be applied to each barcodes.

- parameter --ref_seq_sheet takes as input the path to a text file containing a list of barcode names and the path to their corresponding reference, as in the example below:

| Barcode   | Reference Path      |
|-----------|---------------------|
| barcode01 | /path/to/reference1 |
| barcode02 | /path/to/reference2 |
| barcode03 |                     |
| barcode05 | /path/to/reference5 |


note: if a barcode does not appear in the txt file, or if there is no attributed reference, variant calling is ignored for that barcode (here barcode 03 and 04 for 	instance). Barcode names and paths are tab separated.
WARNING: If a ref_seq_sheet is given but no corresponding sample ID is found, report will not be generated.

All output files are generated in the directory where the command was executed, separated in a different subdirectory for each barcode.


## PRINCIPAL MODIFICATIONS

### Generation of an annotation file in genbank format

A function fa2gbk() was added in the run_plannotate python script at line 65 and called line 150 of the same document

position:
 	/EPI2ME_pipeline/repo_clone/bin/workflow_glue/run_plannotate.py

This function takes as input the fasta file and the annotation bed file of a sample (barcode) and returns the annotation contained in the bed file in a genbank format. the gbk file is generated in the same output directory as all other output files.

### Variant calling process

Four processes were added in main.nf (for variant calling and VCF filtering with Sniffles and Medaka) and are called in the workflow pipeline { } part.If neither --ref_seq_sheet or --ref_seq is used, the process is not called and variant calling does not happen (but a path channel is still generated, pointing to "$projectDir/data/OPTIONAL_FILE" to avoid the workflow to stop.

position:
 	/EPI2ME_pipeline/repo_clone/main.nf

### Feature table

A SplitFeatureTable{} process was added in main.nf to allow splitting the output of the PLannotate process (one huge feature table) into different csv files for each barcodes, later generated in the corresponding subdirectory.

position:
 	/EPI2ME_pipeline/repo_clone/main.nf

### File output

The output{} process was modified to create a subdirectory for each barcode and generate there only .gbk / .fasta / .html and .vcf files
All other files are generated in an annex subdirectory called "other", with all barcodes.
The original output{} process was commented out and left in place as a backup, if needed.

position:
 	/EPI2ME_pipeline/repo_clone/main.nf

### Report general aspect modification

The general aspect of the final html report was modified to fit better the lab's needs. Modifications were done to the report python script.
position:
 	/EPI2ME_pipeline/repo_clone/bin/workflow_glue/report.py

Text can be simply added or modified in the main() part using section.markdown("""
 								## This is some title
 								This is an exemple of text that would then appear in the report
 							    """)

One different and independant report is now generated for each sample (barcode), reordering elements in the for loop line 160. The report now consists of the following sections: Assembly, Annotation table, Quality control, Variant calling (if a reference was given) and a general information part. An if-else statement was introduced lines 170 and 294 to generate a different record in case of a failed analysis, showing only the quality control section.

HES and ONT logos where added at the top and bottom of the report. They are encoded as character strings in base64 in annex txt files, readable by the python script.
position of logos in txt files: /EPI2ME_pipeline/repo_clone/bin/workflow_glue

Some infographics were added in the quality control section, defined in the fastcat_report_tab() function line 68. Colors of the report were also modified are can be changed at these locations:

Temporary color choice for the report:
purple_d : #6205A8
purple_m : #BC61FF
purple_l : #DFB4FF

- Main banner:
/EPI2ME_pipeline/repo_clone/bin/workflow_glue/aplanat_local/data/custom_hes.css
#6205A8

- Table:
/EPI2ME_pipeline/repo_clone/bin/workflow_glue/aplanat_local/report.py
header: #6205A8
even lines: #F3EAF9
hover: #D4A0FF

- infographics:
/EPI2ME_pipeline/repo_clone/bin/workflow_glue/aplanat_local/graphics.py
shape: #6205A8
text: #F3EAF9
icon:#BC61FF

- Plots QC :
/EPI2ME_pipeline/repo_clone/bin/workflow_glue/aplanat_local/components/fastcat.py
line 36 and 65 for read length and read quality plot respectively. 'Colors.HES_purple_d' is called, which can be modified in the _colors class defined line 48 of this document:
/EPI2ME_pipeline/repo_clone/bin/workflow_glue/aplanat_local/util.py
