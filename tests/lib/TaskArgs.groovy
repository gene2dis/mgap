import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

/**
 * The ext.args each task of a process resolved to, from a finished nf-test
 * pipeline run (its trace.csv and work directory). Stubs do not run the tool
 * command, so tests/config/ext_args.config writes the resolved ext.args as a
 * comment into each task's .command.run, and this reads it back.
 */
class TaskArgs {

    /** Returns [ sample (task tag): resolved ext.args ] for every task of the process */
    static Map<String, String> of(metaDir, String process) {
        def workDir = Paths.get(metaDir.toString()).resolveSibling('work')
        def rows = Paths.get(metaDir.toString(), "trace.csv").readLines()
        def header = rows[0].split('\t') as List
        def args = [:]
        rows.drop(1).each { line ->
            def cols = line.split('\t')
            def matcher = cols[header.indexOf('name')] =~ /:([A-Z0-9_]+)(?: \((.+)\))?$/
            if (matcher && matcher[0][1] == process) {
                def run = taskDir(workDir, cols[header.indexOf('hash')]).resolve('.command.run').text
                def recorded = run.readLines().find { it.trim().startsWith('# ext.args:') }
                assert recorded != null : "${process} (${matcher[0][2]}): no ext.args line; is tests/config/ext_args.config loaded?"
                args[matcher[0][2]] = recorded.trim().substring('# ext.args:'.length()).trim()
            }
        }
        return args
    }

    /**
     * Returns [ sample (task tag): sorted names of the folders the task's inputs were staged in ]
     * for every task of the process (inputs declared with stageAs: '<folder>/*')
     */
    static Map<String, List<String>> inputFolders(metaDir, String process) {
        def workDir = Paths.get(metaDir.toString()).resolveSibling('work')
        def rows = Paths.get(metaDir.toString(), "trace.csv").readLines()
        def header = rows[0].split('\t') as List
        def folders = [:]
        rows.drop(1).each { line ->
            def cols = line.split('\t')
            def matcher = cols[header.indexOf('name')] =~ /:([A-Z0-9_]+)(?: \((.+)\))?$/
            if (matcher && matcher[0][1] == process) {
                def dir = taskDir(workDir, cols[header.indexOf('hash')])
                def names = [] as SortedSet
                Files.list(dir).withCloseable { stream ->
                    stream.filter { Files.isDirectory(it) }.forEach { sub ->
                        Files.list(sub).withCloseable { inner ->
                            if (inner.anyMatch { Files.isSymbolicLink(it) }) {
                                names << sub.fileName.toString()
                            }
                        }
                    }
                }
                folders[matcher[0][2]] = names as List
            }
        }
        return folders
    }

    /** Returns the sorted names of the files staged (as symlinks, in any folder) into the tasks of the process */
    static List<String> stagedNames(metaDir, String process) {
        def workDir = Paths.get(metaDir.toString()).resolveSibling('work')
        def rows = Paths.get(metaDir.toString(), "trace.csv").readLines()
        def header = rows[0].split('\t') as List
        def names = [] as SortedSet
        rows.drop(1).each { line ->
            def cols = line.split('\t')
            def matcher = cols[header.indexOf('name')] =~ /:([A-Z0-9_]+)(?: \((.+)\))?$/
            if (matcher && matcher[0][1] == process) {
                Files.walk(taskDir(workDir, cols[header.indexOf('hash')])).withCloseable { stream ->
                    stream.filter { Files.isSymbolicLink(it) }.forEach { names << it.fileName.toString() }
                }
            }
        }
        return names as List
    }

    private static Path taskDir(Path workDir, String hash) {
        def (prefix, rest) = hash.tokenize('/')
        def dir = Files.list(workDir.resolve(prefix)).withCloseable { stream ->
            stream.filter { it.fileName.toString().startsWith(rest) }.findFirst().orElse(null)
        }
        assert dir != null : "no work directory for task hash ${hash} under ${workDir}"
        return dir
    }
}
