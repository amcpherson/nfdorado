nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Define inputs
params.modelQuality = 'dna_r9.4.1_e8.1_sup'
params.methyl_context = 'modbases_5mC_5hmC_CpG_v001'
params.models_dir = '/data1/shahs3/users/schrait/dorado/models'
params.convert_fast5 = true


workflow NFDORADO {

    take:
    samplesheet

    main:

    inputs = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row -> file(row.filename) }

    inputs.view()

    pod5_files = params.convert_fast5
        ? fast5_to_pod5(inputs)
        : inputs

    // // Split by extension
    // split_files = all_files
    //     .branch {
    //         is_pod5: { k, f -> f.name.endsWith('.pod5') }
    //         is_fast5: { k, f -> f.name.endsWith('.fast5') }
    //     }
    // // split_files.is_fast5.view()
    // // split_files.is_pod5.view()

    // // Convert fast5 to pod5
    // fast5_to_pod5(split_files.is_fast5)
    //     .concat(split_files.is_pod5)
    //     .groupTuple()
    //     .set { grouped_pod5_files }

    pod5_files.view()

    basecalling = dorado_basecalling(pod5_files, params.modelQuality, params.methyl_context, params.models_dir)

    sorting = samtools_sort(basecalling.basecalled_bam)

    samtools_stats(sorting.sorted_bam)
}


process fast5_to_pod5 {

    tag "Convert FAST5 to POD5 for ${sample_id}"

    input:
    path(fast5_file)

    output:
    path("converted.pod5")

    script:
    """
    pod5 convert fast5 ${fast5_file} converted.pod5
    """
}


process dorado_basecalling {

    tag "Dorado basecalling on ${pod5_files}"

    input:
    path pod5_files
    val modelQuality
    val methyl_context
    val models_dir

    output:
    path "basecalled.bam", emit: basecalled_bam

    script:
    """
    dorado basecaller --models-directory ${models_dir} ${modelQuality},${methyl_context} ./ --device cuda:all --recursive --verbose > basecalled.bam
    """
}


process samtools_sort {

    tag "Sorting BAM"

    input:
    path basecalled_bam

    output:
    path "sorted.bam", emit: sorted_bam

    script:
    """
    samtools sort -N ${basecalled_bam} > sorted.bam
    """
}


process samtools_stats {

    tag "Generating stats"

    input:
    path sorted_bam

    output:
    path "sorted.bam.stats"

    script:
    """
    samtools stats ${sorted_bam} > sorted.bam.stats
    """
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

