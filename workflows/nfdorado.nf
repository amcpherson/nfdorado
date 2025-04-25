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
params.modelQuality = 'sup@latest'
params.methyl_context = '5mC_5hmC@latest,6mA@latest'
params.models_dir = '/data1/shahs3/users/schrait/dorado/models'
params.convert_fast5 = true

println "Running with the following parameters:"
println "Model Quality: ${params.modelQuality}"
println "Methylation Context: ${params.methyl_context}"
println "Models Directory: ${params.models_dir}"
println "Convert FAST5: ${params.convert_fast5}"


workflow NFDORADO {

    take:
    samplesheet

    main:

    inputs = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row -> file(row.filename) }
        .collect()

    pod5_files = params.convert_fast5
        ? fast5_to_pod5(inputs)
        : merge_pod5(inputs)

    basecalling = dorado_basecalling(pod5_files, params.modelQuality, params.methyl_context, params.models_dir)

    sorting = samtools_sort(basecalling.basecalled_bam)

    samtools_stats(sorting.sorted_bam)
}


process fast5_to_pod5 {

    tag "Convert FAST5 to POD5"

    input:
    path fast5_file

    output:
    path "converted.pod5", emit: pod5_file

    script:
    """
    /home/mcphera1/micromamba/envs/pod5/bin/pod5 convert fast5 ./ -o converted.pod5
    """
}


process merge_pod5 {

    tag "Merge POD5"

    input:
    path pod5_file

    output:
    path "merged.pod5", emit: pod5_file

    script:
    """
    /home/mcphera1/micromamba/envs/pod5/bin/pod5 merge ./ -o merged.pod5
    """
}


process dorado_basecalling {

    tag "Dorado basecalling"

    input:
    path pod5_files
    val modelQuality
    val methyl_context
    val models_dir

    output:
    path "*.bam", emit: basecalled_bam

    script:
    """
    singularity run --nv --bind /data1/shahs3:/data1/shahs3 \
    /data1/shahs3/users/schrait/dorado/ont_methylation.sif \
    dorado basecaller --models-directory ${models_dir} ${modelQuality},${methyl_context} ./ --device cuda:all --recursive --verbose -o ./

    n_bams=\$(ls *.bam | wc -l)
    if [ "\$n_bams" -ne 1 ]; then
        echo "ERROR: Expected exactly one BAM file but found \$n_bams"
        exit 1
    fi
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

