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
params.model_quality = 'sup@latest'
params.methyl_context = '5mC_5hmC@latest,6mA@latest'
params.models_dir = '/data1/shahs3/users/schrait/dorado/models' // TODO fix this
params.convert_fast5 = true

println "Running with the following parameters:"
println "Model Quality: ${params.model_quality}"
println "Methylation Context: ${params.methyl_context}"
println "Models Directory: ${params.models_dir}"
println "Convert FAST5: ${params.convert_fast5}"


// TODO: publishdir
// TODO: name the bam file


workflow NFDORADO {

    take:
    samplesheet
    sample_id

    main:

    inputs = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row -> file(row.filename) }

    pod5_files = params.convert_fast5
        ? fast5_to_pod5(inputs)
        : inputs

    basecalling = dorado_basecalling(pod5_files, params.model_quality, params.methyl_context, params.models_dir)

    merged = samtools_merge(basecalling.basecalled_bam.collect(), sample_id)

    sorting = samtools_sort(merged.merged_bam, sample_id)

    samtools_stats(sorting.sorted_bam, sample_id)
}


process fast5_to_pod5 {
    container 'quay.io/shahlab_singularity/ont_methylation'

    tag "Convert FAST5 to POD5"

    input:
    path fast5_file

    output:
    path "converted.pod5", emit: pod5_file

    script:
    """
    pod5 convert fast5 ./ -o converted.pod5
    """
}


process merge_pod5 {
    container 'quay.io/shahlab_singularity/ont_methylation'

    tag "Merge POD5"

    input:
    path pod5_file

    output:
    path "merged.pod5", emit: pod5_file

    script:
    """
    pod5 merge ./ -o merged.pod5
    """
}


process dorado_basecalling {
    container 'quay.io/shahlab_singularity/ont_methylation'

    tag "Dorado basecalling"

    input:
    path pod5_files
    val model_quality
    val methyl_context
    val models_dir

    output:
    path "*.bam", emit: basecalled_bam

    script:
    """
    dorado basecaller --models-directory ${models_dir} ${model_quality},${methyl_context} ./ --device cuda:all --recursive --verbose -o ./

    n_bams=\$(ls *.bam | wc -l)
    if [ "\$n_bams" -ne 1 ]; then
        echo "ERROR: Expected exactly one BAM file but found \$n_bams"
        exit 1
    fi
    """
}


process samtools_merge {
    container 'quay.io/shahlab_singularity/ont_methylation'

    tag "Merging BAM files"

    input:
    path basecalled_bams, stageAs: "?/*"
    val sample_id

    output:
    path "${sample_id}.bam", emit: merged_bam

    script:
    """
    samtools merge -f ${sample_id}.bam ${basecalled_bams}
    """
}


process samtools_sort {
    container 'quay.io/shahlab_singularity/ont_methylation'

    tag "Sorting BAM"
    publishdir "results"

    input:
    path basecalled_bam, name: 'basecalled.bam'
    val sample_id

    output:
    path "${sample_id}.bam", emit: sorted_bam

    script:
    """
    samtools sort -N ${basecalled_bam} > ${sample_id}.bam
    """
}


process samtools_stats {
    container 'quay.io/shahlab_singularity/ont_methylation'

    tag "Generating stats"
    publishdir "."

    input:
    path sorted_bam

    output:
    path "${sorted_bam}.stats" 

    script:
    """
    samtools stats ${sorted_bam} > ${sorted_bam}.stats
    """
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

