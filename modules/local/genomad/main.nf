process GENOMAD {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::genomad=1.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.5.0--pyhdfd78af_0'':
        'quay.io/biocontainers/genomad:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)


}

//genomad end-to-end -t 40 --cleanup 9941.fna 
//genomad_output /mds_data/dbs/genomad_db