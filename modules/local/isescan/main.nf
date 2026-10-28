process ISESCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isescan:1.7.3--hab16a5f_1' :
        'biocontainers/isescan:1.7.3--hab16a5f_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.isescan.tsv")              , emit: tsv
    tuple val(meta), path("${prefix}.isescan.sum")    , optional: true, emit: summary
    tuple val(meta), path("${prefix}.isescan.gff")    , optional: true, emit: gff
    tuple val(meta), path("${prefix}.isescan.is.fna") , optional: true, emit: is_fna
    tuple val(meta), path("${prefix}.isescan.orf.fna"), optional: true, emit: orf_fna
    tuple val(meta), path("${prefix}.isescan.orf.faa"), optional: true, emit: orf_faa
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    isescan.py \\
        $args \\
        --seqfile ${fasta} \\
        --output isescan \\
        --nthread ${task.cpus}

    # ISEScan names its results after the input file (<fasta>.tsv, .sum, ...) and
    # writes none when it finds no IS element
    for result in tsv sum gff is.fna orf.fna orf.faa; do
        found=\$(find isescan -maxdepth 2 -type f -name "${fasta}.\${result}" | head -n 1)
        if [ -n "\$found" ]; then
            mv "\$found" ${prefix}.isescan.\${result}
        fi
    done
    # a header-only table when no IS element was found, so the sample is reported
    # with none rather than as not run
    if [ ! -f ${prefix}.isescan.tsv ]; then
        printf '%s\\n' "${tsv_header()}" > ${prefix}.isescan.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isescan: \$( isescan.py --version 2>&1 | sed 's/isescan //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '%s\\n' "${tsv_header()}" > ${prefix}.isescan.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isescan: 1.7.3
    END_VERSIONS
    """
}

// The columns of ISEScan 1.7.3's <fasta>.tsv, tab-separated
def tsv_header() {
    return [
        'seqID', 'family', 'cluster', 'isBegin', 'isEnd', 'isLen', 'ncopy4is', 'start1', 'end1', 'start2', 'end2',
        'score', 'irId', 'irLen', 'nGaps', 'orfBegin', 'orfEnd', 'strand', 'orfLen', 'E-value', 'E-value4copy',
        'type', 'ov', 'tir'
    ].join('\t')
}
