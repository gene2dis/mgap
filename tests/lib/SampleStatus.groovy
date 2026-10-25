import java.nio.file.Paths

/**
 * Helpers for the sample QC tests (tests/sample_qc.nf.test): write a samplesheet
 * of local test files before the run, and read pipeline_info/sample_status.tsv
 * after it.
 */
class SampleStatus {

    /** Rows of sample_status.tsv as [ sample: [ column: value ] ] */
    static Map<String, Map<String, String>> of(outputDir) {
        def lines = Paths.get(outputDir.toString(), 'pipeline_info', 'sample_status.tsv').readLines()
        def header = lines[0].split('\t', -1) as List
        return lines.drop(1).collectEntries { line ->
            def row = [ header, line.split('\t', -1) as List ].transpose().collectEntries()
            [ row.sample, row ]
        }
    }

    /** Writes a samplesheet (header line, then one line per sample) and returns its path */
    static String samplesheet(outputDir, String name, String header, List<String> lines) {
        def sheet = Paths.get(outputDir.toString(), name).toFile()
        sheet.parentFile.mkdirs()
        sheet.text = ([ header ] + lines).join('\n') + '\n'
        return sheet.toString()
    }
}
