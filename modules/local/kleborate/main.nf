process KLEBORATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kleborate:3.2.4--pyhdfd78af_0' :
        'biocontainers/kleborate:3.2.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    kleborate \\
        $args \\
        -p kpsc \\
        --trim_headers \\
        -o ./ \\
        --assemblies $fastas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """

    stub:
    // Kleborate names its outputs after the preset, not the sample
    """
    touch klebsiella_pneumo_complex_output.txt
    touch klebsiella_pneumo_complex_hAMRonization_output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: 3.2.4
    END_VERSIONS
    """
}
