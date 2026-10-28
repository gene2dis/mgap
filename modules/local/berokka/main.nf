process BEROKKA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/berokka:0.2.3--hdfd78af_3' :
        'biocontainers/berokka:0.2.3--hdfd78af_3' }"

    input:
    tuple val(meta), path(assembly), path(assembly_info)

    output:
    tuple val(meta), path("${prefix}.berokka.fasta"), emit: fasta
    tuple val(meta), path("${prefix}.berokka.tsv")  , emit: results
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # gzip -f passes an uncompressed assembly through unchanged
    gzip -cdf ${assembly} > assembly.fa

    # Only the contigs Flye marks circular (assembly_info.txt column 4, circ. = Y)
    # go to Berokka: on a linear contig a short end match (an unresolved repeat)
    # would be trimmed and tagged circular
    # (side files are read with getline: with NR == FNR an empty first file
    # would make awk treat the assembly as that file)
    awk -F '\\t' '!/^#/ && \$4 == "Y" { print \$1 }' ${assembly_info} > circular.txt
    awk '
        BEGIN { while ((getline name < "circular.txt") > 0) circular[name] }
        /^>/  { keep = (substr(\$1, 2) in circular) }
        keep
    ' assembly.fa > circular.fa

    if [ -s circular.fa ]; then
        berokka \\
            $args \\
            --outdir berokka \\
            circular.fa
        mv berokka/02.trimmed.fa trimmed.fa
        mv berokka/03.results.tab ${prefix}.berokka.tsv
    else
        : > trimmed.fa
        printf '#sequence\\tstatus\\told_len\\tnew_len\\ttrimmed\\n' > ${prefix}.berokka.tsv
    fi

    # Flye's contig order: circular contigs as Berokka wrote them (a contig it
    # removed as a control sequence stays out), linear contigs unchanged
    # (Berokka's lines are kept as an array with each contig's first and last
    # line: appending them to one string per contig is quadratic, 30 min for a
    # 5.3 Mb chromosome)
    awk '
        BEGIN {
            while ((getline name < "circular.txt") > 0) circular[name]
            while ((getline line < "trimmed.fa") > 0) {
                trimmed[++n] = line
                if (line ~ /^>/) { split(substr(line, 2), word, /[ \\t]/); id = word[1]; first[id] = n }
                last[id] = n
            }
        }
        /^>/ {
            id = substr(\$1, 2)
            linear = !(id in circular)
            if (!linear && id in first) for (i = first[id]; i <= last[id]; i++) print trimmed[i]
        }
        linear
    ' assembly.fa > ${prefix}.berokka.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        berokka: \$( berokka --version 2>&1 | sed 's/berokka //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.berokka.fasta
    printf '#sequence\\tstatus\\told_len\\tnew_len\\ttrimmed\\n' > ${prefix}.berokka.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        berokka: 0.2.3
    END_VERSIONS
    """
}
