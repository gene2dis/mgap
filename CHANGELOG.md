# gene2dis/mgap: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [Unreleased]

Major refactoring to align with nf-core best practices and modern Nextflow patterns.

### `Added`

- Migrated to `nf-schema@2.2.0` plugin for parameter validation (replacing deprecated `nf-validation`)
- Added `PIPELINE_INITIALISATION` and `PIPELINE_COMPLETION` subworkflows for standardized initialization
- Added cloud configuration profiles: `awsbatch`, `googlebatch`, `tower`
- Added Wave container and Fusion filesystem support for cloud execution
- Added nf-test framework with pipeline and module tests
- Added a single GitHub CI workflow (`test.yml`: stub test-profile run + nf-test suite)
- Added `meta.yml` and `environment.yml` to all local modules
- Added stub sections to key modules for fast testing
- Added pre-commit hooks configuration
- Added issue templates for bug reports and feature requests
- Added `mobsuite/recon` (nf-core module) for optional plasmid detection and reconstruction, gated by `--run_mobsuite`
- Added local `sistr` module (v1.1.3) for automatic Salmonella serotype prediction, triggered when MLST scheme is `salmonella`
- Added `SALMONELLA` local subworkflow wrapping SISTR, following the existing `KLEBSIELLA` / `SCCMEC` taxa-specific pattern
- Working test profile with tiny public test data (`-profile test`), plus stub pipeline tests for paired/single-end Illumina, ONT (both assembler modes) and contig inputs
- Single-end Illumina support (detected from the samplesheet)
- Run-level MultiQC report aggregating fastp/fastplong, Kraken2, QUAST and Bakta results (`--multiqc_config`, `--multiqc_title`, `--multiqc_logo`)
- `--medaka_model` parameter for basecaller-matched Medaka polishing
- `--brackendb` is now actually used by Bracken (falls back to the Kraken2 DB)
- Per-sample cross-validation of samplesheet columns against `--seq_type` with clear error messages
- Launch-time warnings for annotation steps skipped because their database was not provided
- Stub blocks in every local module (fast `-stub` runs cover all three input modes)

### `Changed`

- Modernized `nextflow.config` with `resourceLimits`, updated profiles, and cloud parameters
- Refactored `main.nf` to use subworkflow-based initialization pattern
- Updated `workflows/mgap.nf` with explicit closure parameters and version collection
- Modernized `subworkflows/local/illumina.nf` and `ont.nf` with proper channel handling
- Updated `nextflow_schema.json` to JSON Schema draft 2020-12 with organized sections
- Updated documentation in `docs/usage.md` with cloud execution examples
- All local modules now use `${moduleDir}/environment.yml` syntax
- **Updated nf-core modules to latest versions:**
  - `flye`: Updated to v2.9.5, outputs gzipped FASTA
  - `spades`: Updated to v4.1.0, outputs gzipped FASTA
  - `fastp`: Updated module with new input signature (adapter in tuple)
  - `quast`: Updated module with 3-input signature (consensus, reference, gff)
  - `bakta`: Updated module with 6-input signature (+regions, +hmms)
  - `checkm2/predict`: Updated module with database tuple input
  - `amrfinderplus/run`: Updated to v4.2.5 (customized local copy)
  - `genomad/endtoend`: Updated module
  - `gtdbtk/classifywf`: Updated module with 3-input signature
  - `kraken2/kraken2`: Updated module with enhanced test coverage
  - `mlst`: Updated module with environment.yml
  - `multiqc`: Updated module with custom prefix support
  - `seqtk/sample`: Updated module with standard config
- CheckM2, Bakta, AMRFinderPlus and geNomad now run only when their database parameter is set (previously crashed on unset databases)
- GTDB-Tk receives genomes via a batchfile, so it works for any assembly file extension (the `gtdbtk_extension` parameter is gone)
- The Mash/Seqtk coverage-adjustment logic is shared by the Illumina and ONT paths (`COVERAGE_ADJUST` subworkflow) and tolerates unparseable coverage estimates
- The Autocycler container image is pinned by digest; assembler-failure tolerance moved to configuration with loud reporting when all assemblers fail for a sample
- `software_versions.yml` is valid YAML listing each tool once, plus pipeline/Nextflow versions
- **Improved samplesheet handling for contig mode:**
  - `CreateSampleSheet.py` now uses `fasta` column instead of `fastq_1` for contig mode
  - CSV output format now matches data type (sample,fasta vs sample,fastq_1,fastq_2)
  - Schema validation now uses `oneOf` constraint for mutually exclusive fastq/fasta inputs
  - Enabled lenient mode for more flexible input validation
  - Added support for gzipped FASTA files (.fasta.gz, .fa.gz, .fna.gz)
- **Removed deprecated local modules** (replaced by nf-core versions): checkm2, flye, genomad, nanoq, spades (amrfinderplus and kleborate keep customized local copies)

### `Fixed`

- Fixed implicit closure parameter warnings in map/branch operations
- Fixed deprecated `Channel` factory usage (now uses lowercase `channel`)
- Fixed dnaapler version command (was incorrectly calling medaka)
- Fixed AMRFinderPlus null organism error when MLST scheme not in taxa_map
- Fixed the misspelled _Burkholderia pseudomallei_ entry in the MLST-to-AMRFinderPlus organism map, and added verified mappings for seven more species
- Mash coverage estimation no longer crashes on paired-end input (per-file estimates are summed); MLST output is parsed defensively
- Samples no longer vanish silently when dnaapler produces no reoriented assembly (falls back to the unoriented assembly with a warning)
- SISTR's conditionally-produced outputs are optional, so runs no longer fail over missing novel-allele files
- `ConsolidateResults.py`: fixed silently-dropped kleborate/sccmec tables, added ONT support and defensive parsing

### `Removed`

- Azure configuration and the nf-core community scaffolding (issue/PR templates, gitpod/devcontainer, lint workflows)
- Legacy `lib/` Groovy classes and the bundled jar; legacy samplesheet-check chain (`input_check.nf`, `check_samplesheet.py`)
- The no-op email/webhook notification stubs and their parameters/assets
- Dead parameters (`gtdbtk_mash_db`, `gtdbtk_min_perc_aa`, `gtdbtk_min_af`, `unicycler_*`, `antismash_*`, and friends); the schema now matches `nextflow.config` exactly
- Unused nf-core modules (antismash, macrel, and the stale kleborate/amrfinderplus duplicates)
- `accesory_scripts/` renamed to `accessory_scripts/`

### `Dependencies`

- Nextflow: `>=23.04.0` (minimum version)
- nf-schema: `2.2.0` (replacing nf-validation)

## v1.0dev - [date]

Initial release of gene2dis/mgap, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Support for ONT data (December 2023)

### `Fixed`

### `Dependencies`

- Bakta 1.8.2, MLST 2.23, MultiQC 1.18, AMRFinderPlus 3.11.18, geNomad 1.5.2 (December 2023)

### `Deprecated`
