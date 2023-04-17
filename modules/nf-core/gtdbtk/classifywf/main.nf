process GTDBTK_CLASSIFYWF {
    tag "$meta.id"
    label 'process_gtdbtk'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::gtdbtk=2.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.1.1--pyhdfd78af_1' :
        'quay.io/biocontainers/gtdbtk:2.1.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fna, stageAs:"fna-temp/*")
    path db
    //path("database/*")

    output:
    path "gtdbtk.${meta.assembler}-${meta.id}.*.summary.tsv"        , emit: summary
    path "gtdbtk.${meta.assembler}-${meta.id}.*.classify.tree.gz"   , emit: tree
    path "gtdbtk.${meta.assembler}-${meta.id}.*.markers_summary.tsv", emit: markers
    path "gtdbtk.${meta.assembler}-${meta.id}.*.msa.fasta.gz"       , emit: msa
    path "gtdbtk.${meta.assembler}-${meta.id}.*.user_msa.fasta"     , emit: user_msa
    path "gtdbtk.${meta.assembler}-${meta.id}.*.filtered.tsv"       , emit: filtered
    path "gtdbtk.${meta.assembler}-${meta.id}.log"                  , emit: log
    path "gtdbtk.${meta.assembler}-${meta.id}.warnings.log"         , emit: warnings
    path "gtdbtk.${meta.assembler}-${meta.id}.failed_genomes.tsv"   , emit: failed
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //def pplacer_scratch = params.gtdbtk_pplacer_scratch ? "--scratch_dir pplacer_tmp" : ""

    """
    export GTDBTK_DATA_PATH="$db"

    gtdbtk classify_wf \\
        $args \\
        --genome_dir ./fna-temp \\
        -x fa \\
        --prefix "gtdbtk.${meta.id}" \\
        --out_dir "\${PWD}" \\
        --cpus $task.cpus

    gzip "gtdbtk.${meta.id}".*.classify.tree "gtdbtk.${meta.id}".*.msa.fasta
    mv gtdbtk.log "gtdbtk.${meta.id}.log"
    mv gtdbtk.warnings.log "gtdbtk.${meta.id}.warnings.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """

    stub:
    def VERSION = '2.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo "$VERSION")
    END_VERSIONS
    """
}
