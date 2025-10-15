process COMPRESS_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::pigz=2.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b8f6f9860663fb4ab74531715c96bb5f4fe84284:1c63de55bba297d99f73b7a5fd5112290f6064f0-0' :
        'quay.io/biocontainers/mulled-v2-b8f6f9860663fb4ab74531715c96bb5f4fe84284:1c63de55bba297d99f73b7a5fd5112290f6064f0-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.gz"), emit: reads

    script:
        def compression_cmds = reads.collect { read ->
            if (read.name.endsWith('.gz')) {
                def output_name = "${read.name}.linked.gz"
                return "ln -s $read $output_name"
            } else {
                def output_name = "${read.name}.gz"
                return "pigz -p 8 -c $read > $output_name"
            }
        }.join('\n')

        """
        $compression_cmds
        """
}

