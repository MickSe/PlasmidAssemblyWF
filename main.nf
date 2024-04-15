#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'


process checkIfEnoughReads {
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(meta),
            path("input.fastq.gz"),
            path("per-read-stats.tsv"),
            val(approx_size)
        val extra_args
    output:
        tuple val(meta.alias), path("${meta.alias}.fastq.gz"), val(approx_size),
            optional: true, emit: sample
        path "${meta.alias}.stats", emit: stats
        tuple val(meta.alias), env(STATUS), emit: status
    script:
        def expected_depth = "$params.assm_coverage"
        // a little heuristic to decide if we have enough data
        int value = (expected_depth.toInteger()) * 0.8
        int bgzip_threads = task.cpus == 1 ? 1 : task.cpus - 1
    """
    STATUS="Failed due to insufficient reads"
    mv per-read-stats.tsv ${meta.alias}.stats
    fastcat -s ${meta.alias} -r ${meta.alias}.interim $extra_args input.fastq.gz \
    | bgzip -@ $bgzip_threads > interim.fastq.gz
    if [[ "\$(wc -l < "${meta.alias}.interim")" -ge "$value" ]]; then
        mv interim.fastq.gz ${meta.alias}.fastq.gz
        STATUS="Completed successfully"
    fi
    """
}

process filterHostReads {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample_id), path(fastq), val(approx_size)
        path reference
        path regions_bedfile
    output:
        tuple val(sample_id), path("*.filtered.fastq"), val(approx_size), optional: true, emit: unmapped
        path "*.stats", optional: true, emit: host_filter_stats
        tuple val(sample_id), env(STATUS), emit: status
    script:
        def name = sample_id
        def regs = regions_bedfile.name != 'NO_REG_BED' ? regs : 'none'
    """
    STATUS="Failed due to filtered host reads"
    (minimap2 -t $task.cpus -y -ax map-ont $reference $fastq \
        | samtools sort -o ${name}.sorted.aligned.bam -
    samtools index ${name}.sorted.aligned.bam 
    samtools view -b -f 4  ${name}.sorted.aligned.bam > unmapped.bam
    samtools view -b -F 4  ${name}.sorted.aligned.bam > mapped.bam
    samtools fastq unmapped.bam > ${name}.filtered.fastq
    fastcat -s ${name} -r ${name}.interim ${name}.filtered.fastq > /dev/null
    if [[ "\$(wc -l <"${name}.stats")" -ge "1" ]];  then
        mv ${name}.interim ${name}.stats
    fi
    if [[ -f "$regs" ]]; then
        bedtools intersect -a mapped.bam -b $regs -wa \
            | samtools view -bh - > retained.bam
        samtools fastq retained.bam >> ${name}.filtered.fastq
    fi ) && STATUS="Completed successfully"
    """
}

process assembleCore {
    errorStrategy = {task.attempt <= 4 ? 'retry' : 'ignore'}
    maxRetries 4
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample_id), path(fastq), val(approx_size)
    output:
        tuple val(sample_id), path("*.reconciled.fasta"), optional: true, emit: assembly
        tuple val(sample_id), path("*.downsampled.fastq"), optional: true, emit: downsampled
        tuple val(sample_id), env(STATUS), emit: status
    script:
        name = sample_id
        cluster_dir = "trycycler/cluster_001"
        int target = params.assm_coverage * 3
        int min_dep = (params.assm_coverage / 3) * 2
        int min_len = 100
        int max_len = approx_size.toInteger() * 1.2
        int min_q = 7
        int exit_number = task.attempt <= 4 ? 1 : 0
    """

    ############################################################
    # Trimming
    ############################################################
    STATUS="Failed to trim reads"
    (seqkit subseq -j $task.cpus -r $params.trim_length:-$params.trim_length $fastq | \
        seqkit subseq -j $task.cpus -r 1:$max_len | \
        seqkit seq -j $task.cpus -m $min_len -Q $min_q -g > ${name}.trimmed.fastq) \
        && STATUS="Failed to downsample reads" &&

    ############################################################
    # Downsampling
    ############################################################


    (rasusa \
        --coverage $target \
        --genome-size $approx_size \
        --input ${name}.trimmed.fastq > ${name}.downsampled.fastq) \
        && STATUS="Failed to Subset reads" &&

    ############################################################
    # Subsetting
    ############################################################

    (trycycler subsample \
        --count 3 \
        --min_read_depth $min_dep \
        --reads ${name}.downsampled.fastq \
        --out_dir sets \
        --genome_size $approx_size) \
        && STATUS="Failed to assemble using Flye" &&

    ############################################################
    # Assembly
    ############################################################
    (for SUBSET in \$(ls sets/sample_*.fastq)
    do
        SUBSET_NAME=\$(basename -s .fastq \$SUBSET)
        flye \
            --${params.flye_quality}\
            \${SUBSET} \
            --threads $task.cpus \
            --genome-size $approx_size \
            --out-dir "assm_\${SUBSET_NAME}" \
            --meta
             
        mv assm_sample_0*/assembly.fasta "assm_\${SUBSET_NAME}/\${SUBSET_NAME}_assembly.fasta" 
    done) && STATUS="Failed to trim Assembly" &&

    ############################################################
    # Trim assemblies
    ############################################################

    (for assembly in \$(ls assm_sample_0*/*assembly.fasta)
    do  
        echo \$assembly
        assembly_name=\$(basename -s .fasta \$assembly)
        ass_stats=\$(dirname \$assembly)/assembly_info.txt
        workflow-glue deconcatenate \
            \$assembly \
            -o \${assembly_name}.deconcat.fasta
    done
    ls *.deconcat.fasta > /dev/null 2>&1) \
    && STATUS="Failed to reconcile assemblies" &&


    ############################################################
    # Reconciliation
    ############################################################

    (trycycler cluster \
        --assemblies *.deconcat.fasta \
        --reads ${name}.downsampled.fastq \
        --out_dir trycycler) &&
    (trycycler reconcile \
        --reads ${name}.downsampled.fastq \
        --cluster_dir $cluster_dir \
        --max_trim_seq_percent 20 \
        --max_add_seq_percent 10) &&
    (trycycler msa --cluster_dir $cluster_dir) &&
    (trycycler partition --reads ${name}.downsampled.fastq --cluster_dirs $cluster_dir) &&
    (trycycler consensus --cluster_dir $cluster_dir)

    ############################################################
    # Exit handling
    ############################################################

    if [ ! -f "${cluster_dir}/7_final_consensus.fasta" ]; then
        if ls ${cluster_dir}/1_contigs/*.fasta 1> /dev/null 2>&1; then
            STATUS="Completed but failed to reconcile"
            (seqkit sort ${cluster_dir}/1_contigs/*.fasta --by-length \
                | seqkit head -n 1 > ${name}.reconciled.fasta) \
                && echo "Trycycler failed, outputting un-reconciled assembly"
        elif [ "$exit_number" == "1" ]; then
            echo \$STATUS
            echo "Assembly failed, retrying process"
            exit 1
        elif [ "$exit_number" == "0" ]; then
            echo \$STATUS
            echo "Failed final attempt"
        fi
    else
        mv ${cluster_dir}/7_final_consensus.fasta ${name}.reconciled.fasta
        STATUS="Completed successfully"
    fi
    """
}

process lookup_medaka_model {
    label "wfplasmid"
    input:
        path("lookup_table")
        val basecall_model
    output:
        stdout
    shell:
    '''
    medaka_model=$(workflow-glue resolve_medaka_model lookup_table '!{basecall_model}')
    echo $medaka_model
    '''
}

process medakaPolishAssembly {
    label "medaka"
    cpus params.threads
    input:
        tuple val(sample_id), path(draft), path(fastq), val(medaka_model)
    output:
        tuple val(sample_id), path("*.final.fasta"), emit: polished
        tuple val(sample_id), env(STATUS), emit: status
    script:
    def model = medaka_model

    """
    STATUS="Failed to polish assembly with Medaka"
    medaka_consensus -i "${fastq}" -d "${draft}" -t $task.cpus -f -o . -m $medaka_model
    echo ">${sample_id}" >> "${sample_id}.final.fasta"
    sed "2q;d" consensus.fasta >> "${sample_id}.final.fasta"
    mv consensus.fasta "${sample_id}.final.fastq"
    STATUS="Completed successfully"
    """
}

process downsampledStats {
    label "wfplasmid"
    cpus 1
    input:
        tuple val(sample_id), path(sample)
    output:
        path "*.stats", optional: true
    """
    fastcat -s ${sample_id} -r ${sample_id}.downsampled $sample > /dev/null
    if [[ "\$(wc -l <"${sample_id}.downsampled")" -ge "2" ]];  then
        mv ${sample_id}.downsampled ${sample_id}.stats
    fi
    """
}

process findPrimers {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus 1
    input:
        path primers
        tuple val(sample_id), path(sequence)
    output:
        path "*.bed", optional: true
    shell:
    '''
    cat !{sequence} | seqkit amplicon -p !{primers} -m 3 -j !{task.cpus} --bed >> !{sample_id}.interim
    if [[ "$(wc -l <"!{sample_id}.interim")" -ge "1" ]];  then
        mv !{sample_id}.interim !{sample_id}.bed
    fi
    '''
}

process medakaVersion {
    label "medaka"
    output:
        path "medaka_version.txt"
    """
    medaka --version | sed 's/ /,/' >> "medaka_version.txt"
    """
}

process getVersions {
    label "wfplasmid"
    cpus 1
    input:
        path "input_versions.txt"
    output:
        path "versions.txt"
    script:
    """
    cat "input_versions.txt" >> "versions.txt"
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    seqkit version | sed 's/ /,/' >> versions.txt
    trycycler --version | sed 's/ /,/' >> versions.txt
    porechop --version | sed 's/^/porechop,/'  >> versions.txt
    bedtools --version | sed 's/ /,/' >> versions.txt
    flye --version | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    last --version | sed 's/ /,/' >> versions.txt
    rasusa --version | sed 's/ /,/' >> versions.txt
    python -c "import spoa; print(spoa.__version__)" | sed 's/^/spoa,/'  >> versions.txt
    """
}

process getParams {
    label "wfplasmid"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process runPlannotate {
    label "wfplasmid"
    cpus 1
    input:
        path annotation_database
        path "assemblies/*"
        path final_status
    output:
        path "feature_table.txt", emit: feature_table
        path "plannotate.json", emit: json
        path "*annotations.bed", optional: true, emit: annotations
        path "*annotations.gbk", optional: true, emit: gbk
        path "plannotate_report.json", emit: report
    script:
        def database =  annotation_database.name.startsWith('OPTIONAL_FILE') ? "Default" : "${annotation_database}"
    """
    if [ -e "assemblies/OPTIONAL_FILE" ]; then
        assemblies=""
    else
        assemblies="--sequences assemblies/"
    fi
    workflow-glue run_plannotate \$assemblies --database $database
    """
}

process inserts {
    label "wfplasmid"
    cpus 1
    input:
         path "primer_beds/*"
         path "assemblies/*"
         path align_ref
    output:
        path "inserts/*", optional: true, emit: inserts
        path "*.json", emit: json
    script:
        def ref =  align_ref.name.startsWith('OPTIONAL_FILE') ? '' : "--reference ${align_ref}"
    """
    if [ -e "primer_beds/OPTIONAL_FILE" ]; then
        inserts=""
    else
        inserts="--primer_beds primer_beds/*"
    fi
    workflow-glue find_inserts \$inserts $ref  
    """
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
// added processes - Mickaël - March 2024

// Variant calling for SV using medaka
process variantCallingSniffles {

    label "wfplasmid"
    cpus 1
    input:
        tuple val(sample_id), path(ref_seq), path(fastq), val(approx_size)

    output:
        tuple val(sample_id), path(ref_seq), path("${sample_id}.vcfTable.csv"), optional: true, emit: vcf_table
        path "${sample_id}.variantCalling.vcf", optional: true, emit: vcf_file
        
    shell:

    """
    echo \$PWD

    minimap2 -a !{ref_seq} !{fastq} > !{sample_id}.alignement.sam

    samtools view -h -bS !{sample_id}.alignement.sam > !{sample_id}.alignement.bam

    samtools sort -O bam -T !{sample_id}.alignement.sort -o !{sample_id}.alignement.sort.bam !{sample_id}.alignement.bam

    samtools index !{sample_id}.alignement.sort.bam 


    sniffles --reference !{ref_seq}  --allow-overwrite -i !{sample_id}.alignement.sort.bam  -v !{sample_id}.variantCalling.vcf


    awk -F '\t' '
        !/^#/ {
            split(\$8, info_fields, ";")
            svtype = svlen = support = ""
            for (i in info_fields) {
                split(info_fields[i], info_pair, "=")
                if (info_pair[1] == "SVTYPE") svtype = info_pair[2]
                else if (info_pair[1] == "SVLEN") svlen = info_pair[2]
                else if (info_pair[1] == "SUPPORT") support = info_pair[2]
            }
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, svtype, svlen, support
        }
    ' '!{sample_id}.variantCalling.vcf' > !{sample_id}.vcfTable.csv

    """
}

// variant calling for shorter variants using sniffles
process variantCallingMedaka {

    label "medaka"
    cpus 1
    input:
    tuple val(sample_id), path(ref_seq), path(fastq), val(approx_size), val(medaka_model)

    output:
    tuple val(sample_id), path(ref_seq), path("${sample_id}.vcfTable_SN.vcf"), optional: true, emit: vcf_table
    path "${sample_id}.vcffile_SN.vcf", optional: true, emit: vcf_file_SN
        
    shell:

    """

    medaka_haploid_variant -i !{fastq} \
    -r !{ref_seq} \
    -m !{medaka_model} \

    mv medaka/medaka.vcf !{sample_id}.vcffile_SN.vcf

    awk -F '\t' '
        !/^#/ {
            split(\$8, info_fields, ";")
            svtype = svlen = support = ""
            for (i in info_fields) {
                split(info_fields[i], info_pair, "=")
                if (info_pair[1] == "SVTYPE") svtype = info_pair[2]
                else if (info_pair[1] == "SVLEN") svlen = info_pair[2]
                else if (info_pair[1] == "SUPPORT") support = info_pair[2]
            }
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7
        }
    ' '!{sample_id}.vcffile_SN.vcf' > !{sample_id}.vcfTable_SN.vcf

    """
}

// filtering of Sniffles VCF output, with a min quality score and homopolymer threshold of 10 (hardcoded here)
process FilterVCF_Sniffles {

    label "wfplasmid"
    cpus 1
    input:
        tuple val(sample_id), path(ref_seq), path(VCF)

    output:
        path "*.filtered_sniffles.vcf"
        
    shell:

    """

    workflow-glue vcf_filter \
    --minscore 10 \
    --fasta !{ref_seq} \
    --type fasta \
    --vcf !{VCF} \
    --homopolymer_threshold 10 \
    --output_file !{sample_id}.filtered_sniffles.vcf

    """
}

// filtering of Medaka VCF output, with a min quality score and homopolymer threshold of 10 (hardcoded here)
process FilterVCF_Medaka {

    label "wfplasmid"
    cpus 1
    input:
        tuple val(sample_id), path(ref_seq), path(VCF)

    output:
        path "*.filtered_medaka.vcf"
        
    shell:

    """

    workflow-glue vcf_filter \
    --minscore 10 \
    --fasta !{ref_seq} \
    --type fasta \
    --vcf !{VCF} \
    --homopolymer_threshold 10 \
    --output_file !{sample_id}.filtered_medaka.vcf

    """
}

// Split CSV files based on sample_name (used to generate a separate feature table for each barcode)
process SplitFeatureTable {
    
    input:
    path csv_file
    
    output:
    path "*.FeatureTable.csv"
    
    shell:
    """
    mapfile -t csv_lines < $csv_file

    column_names="\${csv_lines[0]}"
    echo "Column names: \$column_names"

    for line in "\${csv_lines[@]}"; do

        first_column="\$(echo "\$line" | cut -d ',' -f 1)"

        if [ "\$first_column" != "Sample_name" ]; then
            output_file="\${first_column}.FeatureTable.csv"

            if [[ ! -f "\$output_file" ]]; then
                echo "column names was written at the beginning of the file for \$first_column"
                echo "\$column_names" > "\$output_file"
            fi
            echo "\$line" >> "\$output_file"
        fi
    done

    """
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

process report {
    label "wfplasmid"
    cpus 1
    input:
        path "downsampled_stats/*"
        path final_status
        path "per_barcode_stats/*"
        path "host_filter_stats/*"
        path "versions/*"
        path "params.json"
        path plannotate_json
        path inserts_json
        path lengths
        val rep_aplanat
        val vcfTable
        val vcfTable_SN
    output:
        path "*.report.html", emit: html
        path "sample_status.txt", emit: sample_stat
        path "inserts/*", optional: true, emit: inserts
    script:
        report_name = "wf-clone-validation-"
    """
    workflow-glue report \
    --downsampled_stats downsampled_stats/* \
    --revision $workflow.revision \
    --commit $workflow.commitId \
    --status $final_status \
    --per_barcode_stats per_barcode_stats/* \
    --host_filter_stats host_filter_stats/* \
    --params params.json \
    --versions versions \
    --report_name $report_name \
    --plannotate_json $plannotate_json \
    --lengths $lengths \
    --inserts_json $inserts_json \
    --aplanat $rep_aplanat \
    --VCFsniffles $vcfTable \
    --VCFmedaka $vcfTable_SN

    """
}

workflow pipeline {
    take:
        samples
        host_reference
        regions_bedfile
        database
        primers
        align_ref
        min_read_length
        max_read_length
        ref_seq_sheet
        ref_seq

    main:

        // remove samples that didn't have sequences (i.e. metamap entries without
        // corresponding barcode sub-directories) and get the per-read stats file from
        // the fastcat stats dir
        samples = samples
        | filter { it[1] }
        | map { it[2] = it[2].resolve("per-read-stats.tsv"); it }

        // Min/max filter reads
        fastcat_extra_args = "-a $min_read_length -b $max_read_length"
        
        // drop samples with too low coverage
        sample_fastqs = checkIfEnoughReads(samples, fastcat_extra_args)

        // Optionally filter the data, removing reads mapping to
        // the host or background genome
        if (host_reference.name != "NO_HOST_REF") {
            filtered = filterHostReads(
                    sample_fastqs.sample, host_reference, regions_bedfile)
            samples_filtered = filtered.unmapped
            updated_status = filtered.status
            filtered_stats = filtered.host_filter_stats.collect()
                             .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        }
        else {
            samples_filtered = sample_fastqs.sample
            updated_status = sample_fastqs.status
            filtered_stats = file("$projectDir/data/OPTIONAL_FILE")
        }
       
        // Core assembly and reconciliation
        assemblies = assembleCore(samples_filtered)
        
        named_drafts = assemblies.assembly.groupTuple()
        named_samples = assemblies.downsampled.groupTuple()
        named_drafts_samples = named_drafts.join(named_samples)

        if(params.medaka_model) {
            log.warn "Overriding Medaka model with ${params.medaka_model}."
            medaka_model = Channel.fromList([params.medaka_model])
        }
        else {
            // map basecalling model to medaka model
            lookup_table = Channel.fromPath("${projectDir}/data/medaka_models.tsv", checkIfExists: true)
            medaka_model = lookup_medaka_model(lookup_table, params.basecaller_cfg)
        }
        // Polish draft assembly
        polished = medakaPolishAssembly(named_drafts_samples.combine(medaka_model))
       
        // Concat statuses and keep the last of each
        final_status = sample_fastqs.status.concat(updated_status)
        .concat(assemblies.status).concat(polished.status).groupTuple()
        .map { it -> it[0].toString() + ',' + it[1][-1].toString() }
        final_status = final_status.collectFile(name: 'final_status.csv', newLine: true)
    
        downsampled_stats = downsampledStats(assemblies.downsampled)

        primer_beds = findPrimers(primers, polished.polished)
        medaka_version = medakaVersion()
        software_versions = getVersions(medaka_version)
        workflow_params = getParams()

        annotation = runPlannotate(
            database, polished.polished.map { it -> it[1] }.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            final_status)

        insert = inserts(primer_beds.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            polished.polished.map { it -> it[1] }.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            align_ref)

        //////////////////////////////////////////
        //             ADDED PROCESS            //
        //////////////////////////////////////////
        // note: if a barcode does not appear in the txt file (when using --ref_seq_sheet), or if there is no attributed
        // reference, variant calling is ignored for that barcode
        // WARNING: edge case if --ref_seq_sheet is given and not empty, but does not contain any of the expected sample ID, 
        // further processes are not called and report is not generated !

        VC = false
        if (ref_seq_sheet != null) { //check if --ref_seq_sheet is given
            Channel
                .fromPath( ref_seq_sheet )
                .splitCsv( header: false, sep: '\t' )
                .map { row ->
                    if (row.size() >= 2) {
                        tuple(row[0], row[1])
                    } 
                }
                .set { align_ref } // Channel with tuple val(sample_id) path(ref_seq)
            
            align_ref.join(samples_filtered).set {samplesAndRef}
            variant_calling_sniffles = variantCallingSniffles(samplesAndRef)
            variant_calling_sniffles_filtered = FilterVCF_Sniffles(variant_calling_sniffles.vcf_table)

            variant_calling_medaka = variantCallingMedaka(samplesAndRef.combine(medaka_model))
            variant_calling_medaka_filtered = FilterVCF_Medaka(variant_calling_medaka.vcf_table)

            VC = true

        } else if (ref_seq != null) { //if not check if --ref_seq is given and use the same ref for all barcodes

            samples_filtered.map { it -> tuple(it[0], ref_seq, it[1], it[2]) }
                            .set { samplesAndRef }
            variant_calling_sniffles = variantCallingSniffles(samplesAndRef)
            variant_calling_sniffles_filtered = FilterVCF_Sniffles(variant_calling_sniffles.vcf_table)

            variant_calling_medaka = variantCallingMedaka(samplesAndRef.combine(medaka_model))
            variant_calling_medaka_filtered = FilterVCF_Medaka(variant_calling_medaka.vcf_table)

            VC = true

        } else { // If not do not perform variant calling but Channel is still generated (point to an empty directory - Channel needs to exist for the report process)
            Channel.from("$projectDir/data/OPTIONAL_FILE").set {variant_calling_sniffles}
            Channel.from("$projectDir/data/OPTIONAL_FILE").set {variant_calling_medaka}
        }

        // Split the feature table csv file to get a different file for each sample ID (barcode)
        feature_table = SplitFeatureTable(annotation.feature_table.collect())


        if (VC == true) { // sends different channels to the report and output processes depending if variant calling was applied or not

            variant_calling_sniffles_filtered.collect().map {items -> 
                                                tuple(items)
                                                }.set { vcf_path_sniffles }

            variant_calling_medaka_filtered.collect().map {items -> 
                                                tuple(items)
                                                }.set { vcf_path_medaka }

            report = report(
                downsampled_stats.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
                final_status,
                sample_fastqs.stats.collect(),
                filtered_stats,
                software_versions.collect(),
                workflow_params,
                annotation.report,
                insert.json,
                annotation.json,
                "$projectDir/bin/workflow_glue",
                vcf_path_sniffles,
                vcf_path_medaka)

            results = polished.polished.map { it -> it[1] }.concat(
                report.html,
                report.sample_stat,
                feature_table,
                insert.inserts,
                annotation.json,
                annotation.annotations,
                annotation.gbk,
                variant_calling_sniffles_filtered,
                variant_calling_sniffles.vcf_file,
                variant_calling_medaka_filtered,
                variant_calling_medaka.vcf_file_SN,
                workflow_params)

        } else {
            report = report(
                downsampled_stats.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
                final_status,
                sample_fastqs.stats.collect(),
                filtered_stats,
                software_versions.collect(),
                workflow_params,
                annotation.report,
                insert.json,
                annotation.json,
                "$projectDir/bin/workflow_glue",
                variant_calling_sniffles,
                variant_calling_medaka)

            results = polished.polished.map { it -> it[1] }.concat(
                report.html,
                report.sample_stat,
                feature_table,
                insert.inserts,
                annotation.json,
                annotation.annotations,
                annotation.gbk,
                workflow_params)
        }
    emit:
        results
        telemetry = workflow_params


}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.

/* Original unmodified process (as a backup) - outputs all files in the same directory
process output {
    // publish inputs to output directory
    label "wfplasmid"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", saveAs: {
        f -> params.prefix ? "${params.prefix}-${f}" : "${f}" }
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"

    """
}
*/

// modified process, creates an output dir for each barcode and moves only fasta, gbk, csv and html files there.   
// All other output files are moved to a directory called other - Mickaël - March 2024
process output {
    label "wfplasmid"

    input:
        file(fname)
        path dir

    output:
        stdout
    
    shell:
    """
    for file in $fname; do
        full_name=\${file#[}
        full_name=\${full_name%,}
        full_name=\${full_name%]}
        sample_id=\$(echo "\$full_name" | awk -F'.' '{print \$1}')
        file_extension=\$(echo "\$full_name" | awk -F '.' '{print \$NF}')

        if [ "\$file_extension" == "gbk" ] || [ "\$file_extension" == "fasta" ] || [ "\$file_extension" == "html" ] || [ "\$file_extension" == "csv" ]; then

            mkdir -p $dir/output/\$sample_id
            cp \$full_name $dir/output/\$sample_id/\$full_name

        else

            mkdir -p $dir/output/other
            cp \$full_name $dir/output/other/\$full_name

        fi

    done
    """
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    // calculate min and max read length for filtering with `fastcat`
    List approx_sizes; int min_read_length, max_read_length
    if(params.approx_size_sheet) {
        // the file provided with `--approx_size_sheet` contains a size estimate for
        // each sample, but we will filter all samples with the same parameters)
        approx_sizes = file(params.approx_size_sheet).splitCsv(header:true)
        min_read_length = approx_sizes.collect { it["approx_size"].toInteger() }.min()
        max_read_length = approx_sizes.collect { it["approx_size"].toInteger() }.max()
    } else {
        // we only got a single size estimate --> take as max and min
        min_read_length = max_read_length = params.approx_size
    }
    // +/- 50% margins for read length thresholds
    min_read_length *= 0.5
    max_read_length *= 1.5

    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "fastcat_stats": true
        ])

    // add the size estimates to the channel with the samples
    if (params.approx_size_sheet) {
        samples = samples
        | map { [it[0]["alias"], *it] }
        | join(Channel.of(*approx_sizes) | map { [it["alias"], it["approx_size"]] } )
        | map { it[1..-1] }
        | ifEmpty {
            error "The sample aliases in the CSV file provided with " +
                "`--approx_size_sheet` don't match up with the sample names. Have " +
                "you made sure the 'alias' column contains the correct sample names? " +
                "You could have also forgotten to provide a sample sheet alongside " +
                "the approx. size sheet. In case you are not using a sample sheet, " +
                "make sure the 'alias' column in the size sheet matches the barcode " +
                "directory names."
        }
    } else {
        samples = samples | map { [*it, params.approx_size] }
    }

    // assigns params.host_reference to host_reference if it has a value, 'NO_HOST_REF' otherwise
    host_reference = params.host_reference ?: 'NO_HOST_REF'
    host_reference = file(host_reference, checkIfExists: host_reference == 'NO_HOST_REF' ? false : true)
    regions_bedfile = params.regions_bedfile ?: 'NO_REG_BED'
    regions_bedfile = file(regions_bedfile, checkIfExists: regions_bedfile == 'NO_REG_BED' ? false : true)
    primer_file = file("$projectDir/data/OPTIONAL_FILE")
    if (params.primers != null){
        primer_file = file(params.primers, type: "file")
    }
    align_ref = file("$projectDir/data/OPTIONAL_FILE")
    if (params.reference != null){
        align_ref = file(params.reference, type: "file")
    }
    database = file("$projectDir/data/OPTIONAL_FILE")
    if (params.db_directory != null){
         database = file(params.db_directory, type: "dir")

    }


    ref_seq_sheet = params.ref_seq_sheet
    if (ref_seq_sheet != null) {
        ref_seq_sheet = file(ref_seq_sheet, checkIfExists: ref_seq_sheet == null ? false : true)
    }

    ref_seq = params.ref_seq
    if (ref_seq != null) {
        ref_seq = file(ref_seq, type: 'file')
    }


    // Run pipeline
    results = pipeline(
        samples,
        host_reference,
        regions_bedfile,
        database,
        primer_file,
        align_ref,
        min_read_length,
        max_read_length,
        ref_seq_sheet,
        ref_seq)


    output(results[0],"$launchDir")
    // previously "$projectDir"


}


if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }
    
    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
