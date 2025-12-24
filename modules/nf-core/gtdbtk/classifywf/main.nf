process GTDBTK_CLASSIFYWF {
    tag "batch_${genomes.size()}_genomes"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.3.2--pyhdfd78af_0' :
        'biocontainers/gtdbtk:2.3.2--pyhdfd78af_0' }"

    input:
    val(metas)
    path(genomes)
    tuple val(db_name), path(database)
    path(mash_db)

    output:
    path("gtdbtk.batch.*.summary.tsv")         , emit: summary
    path("gtdbtk.batch.*.classify.tree.gz")    , emit: tree, optional: true
    path("gtdbtk.batch.*.markers_summary.tsv") , emit: markers, optional: true
    path("gtdbtk.batch.*.msa.fasta.gz")        , emit: msa, optional: true
    path("gtdbtk.batch.*.user_msa.fasta.gz")   , emit: user_msa, optional: true
    path("gtdbtk.batch.*.filtered.tsv")        , emit: filtered, optional: true
    path("gtdbtk.batch.failed_genomes.tsv")    , emit: failed, optional: true
    path("gtdbtk.batch.log")                   , emit: log
    path("gtdbtk.batch.warnings.log")          , emit: warnings
    path("versions.yml")                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pplacer_scratch = params.gtdbtk_pplacer_scratch ? true : false
    def mash_mode = mash_db ? "--mash_db ${mash_db}" : "--skip_ani_screen"
    def extension = params.gtdbtk_extension ?: "fa"

    """
    # Create bins directory and stage all genomes
    mkdir -p bins
    for genome in ${genomes}; do
        ln -s "\$(readlink -f "\${genome}")" bins/
    done

    export GTDBTK_DATA_PATH="\${PWD}/${database}"
    if [ "${pplacer_scratch}" = "true" ] ; then
        mkdir pplacer_tmp
    fi

    gtdbtk classify_wf \\
        $args \\
        --genome_dir bins \\
        --prefix "gtdbtk.batch" \\
        --out_dir "\${PWD}" \\
        --cpus $task.cpus \\
        --extension ${extension} \\
        $mash_mode \\
        ${pplacer_scratch == true ? "--scratch_dir pplacer_tmp" : ""} \\
        --min_perc_aa $params.gtdbtk_min_perc_aa \\
        --min_af $params.gtdbtk_min_af

    mv classify/* .

    mv identify/* .

    mv align/* .

    mv gtdbtk.log "gtdbtk.batch.log"

    mv gtdbtk.warnings.log "gtdbtk.batch.warnings.log"

    find -name gtdbtk.batch.*.classify.tree | xargs -r gzip # do not fail if .tree is missing

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """

    stub:
    def VERSION = '2.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch gtdbtk.batch.stub.summary.tsv
    touch gtdbtk.batch.stub.classify.tree.gz
    touch gtdbtk.batch.stub.markers_summary.tsv
    touch gtdbtk.batch.stub.msa.fasta.gz
    touch gtdbtk.batch.stub.user_msa.fasta.gz
    touch gtdbtk.batch.stub.filtered.tsv
    touch gtdbtk.batch.log
    touch gtdbtk.batch.warnings.log
    touch gtdbtk.batch.failed_genomes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo "$VERSION")
    END_VERSIONS
    """
}
