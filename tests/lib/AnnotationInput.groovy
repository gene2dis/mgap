import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Checks which FASTA each annotation and typing task staged, from a finished
 * nf-test pipeline run (its trace.csv and work directory).
 *
 * Every task of these processes must stage the uncompressed SEQKIT_FILTER
 * output of its own sample (GTDB-Tk, one batch task: of every sample) and no
 * gzipped file.
 */
class AnnotationInput {

    static final List<String> TOOLS = [
        'CHECKM2', 'MLST', 'BAKTA', 'GTDBTK', 'GENOMAD', 'MOBSUITE_RECON',
        'RGI_MAIN', 'KLEBORATE', 'SCCMEC', 'SISTR'
    ]

    /** Returns [ tools: processes that ran, errors: what was wrong ] */
    static Map check(metaDir) {
        def workDir = Paths.get(metaDir.toString()).resolveSibling('work')
        def rows = Paths.get(metaDir.toString(), "trace.csv").readLines()
        def header = rows[0].split('\t') as List
        def tasks = rows.drop(1).collect { line ->
            def cols = line.split('\t')
            def name = cols[header.indexOf('name')]
            def matcher = name =~ /:([A-Z0-9_]+)(?: \((.+)\))?$/
            [
                process: matcher ? matcher[0][1] : name,
                sample : matcher ? matcher[0][2] : null,
                dir    : taskDir(workDir, cols[header.indexOf('hash')])
            ]
        }

        // sample -> the SEQKIT_FILTER output it must reach the tools as
        def filtered = tasks
            .findAll { it.process == 'SEQKIT_FILTER' }
            .collectEntries { [ (it.sample): it.dir.resolve("${it.sample}.fasta".toString()).toRealPath() ] }

        def errors = []
        def ran = [] as Set
        tasks.findAll { it.process in TOOLS }.each { task ->
            ran << task.process
            def staged = stagedLinks(task.dir)
            staged.findAll { it.fileName.toString().endsWith('.gz') }.each { link ->
                errors << "${task.process} (${task.sample}): staged gzipped ${link.fileName}"
            }
            def targets = staged.collect { it.toRealPath() }
            // GTDB-Tk is one batch task (tag gtdbtk_batch) over every sample
            def expected = task.process == 'GTDBTK' ? filtered.values() : [ filtered[task.sample] ]
            expected.each { fasta ->
                if (!fasta || !(fasta in targets)) {
                    errors << "${task.process} (${task.sample}): did not stage ${fasta}"
                } else if (!fasta.text.startsWith('>')) {
                    errors << "${task.process} (${task.sample}): ${fasta} is not plain FASTA"
                }
            }
        }
        return [ tools: ran, errors: errors ]
    }

    private static Path taskDir(Path workDir, String hash) {
        def (prefix, rest) = hash.tokenize('/')
        return Files.list(workDir.resolve(prefix)).find { it.fileName.toString().startsWith(rest) }
    }

    // symlinks Nextflow staged into the task directory, in any subdirectory
    private static List<Path> stagedLinks(Path dir) {
        def links = []
        Files.walk(dir).withCloseable { stream ->
            stream.filter { Files.isSymbolicLink(it) }.forEach { links << it }
        }
        return links
    }
}
