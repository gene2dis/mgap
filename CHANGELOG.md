# gene2dis/mgap: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [Unreleased]

Major refactoring to align with nf-core best practices and modern Nextflow patterns.

### `Added`

- Migrated to `nf-schema@2.2.0` plugin for parameter validation (replacing deprecated `nf-validation`)
- Added `PIPELINE_INITIALISATION` and `PIPELINE_COMPLETION` subworkflows for standardized initialization
- Added cloud configuration profiles: `mgap_awsbatch`, `mgap_googlebatch`, `tower` (prefixed so they do not merge with nf-core's institutional `awsbatch`/`googlebatch` profiles)
- Added Wave container and Fusion filesystem support for cloud execution
- Added nf-test framework with pipeline and module tests
- Added a single GitHub CI workflow (`test.yml`: stub test-profile run + nf-test suite)
- CI lint job (`pre-commit run --all-files`, `nextflow lint`); the stub run also runs on the minimum Nextflow version (26.04.0)
- Stub tests for species routing (MLST scheme → Kleborate, SCCmec, AMRFinderPlus organism), coverage adjustment, Kraken2/Bracken, every Autocycler module, `MASH_SKETCH`, `RGI_MAIN` and the `COVERAGE_ADJUST` subworkflow
- Added `meta.yml` and `environment.yml` to all local modules
- Added stub sections to key modules for fast testing
- Added pre-commit hooks configuration
- Added issue templates for bug reports and feature requests
- Added `mobsuite/recon` (nf-core module) for optional plasmid detection and reconstruction, gated by `--run_mobsuite`
- Added local `sistr` module (v1.1.3) for automatic Salmonella serotype prediction, triggered when MLST scheme is `senterica_achtman_2`
- Added the `SPECIES_TYPING` local subworkflow: one place that routes samples by MLST-derived organism to Kleborate, SCCmec and SISTR (replaces the one-module `KLEBSIELLA` and `SALMONELLA` subworkflows and the bare `SCCMEC` call)
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
- `software_versions.yml` is valid YAML listing each tool once (MultiQC included), plus pipeline/Nextflow versions; all but MultiQC's own appear as a Software Versions section in the MultiQC report (grouped by process name, so tools MultiQC also detects are listed once)
- `AMRFINDERPLUS_RUN` reads the tool and database versions once, no longer emits the unused `VER`/`DBVER` env outputs and no longer writes a file named `stdout`
- **Improved samplesheet handling for contig mode:**
  - `CreateSampleSheet.py` now uses `fasta` column instead of `fastq_1` for contig mode
  - CSV output format now matches data type (sample,fasta vs sample,fastq_1,fastq_2)
  - Schema validation now uses `oneOf` constraint for mutually exclusive fastq/fasta inputs
  - Enabled lenient mode for more flexible input validation
  - Added support for gzipped FASTA files (.fasta.gz, .fa.gz, .fna.gz)
- **Removed deprecated local modules** (replaced by nf-core versions): checkm2, flye, genomad, nanoq, spades (amrfinderplus and kleborate keep customized local copies)
- `conf/ont.config` is merged into `conf/modules.config`, which now holds the options of every process
- Bakta no longer publishes its `versions.yml` into `<sample>/annotation/bakta/` (versions are collected in `pipeline_info/software_versions.yml`)
- RGI loads CARD once per run: a new `RGI_LOAD` step runs `rgi load` and builds the BLAST/DIAMOND databases, and `RGI_MAIN` works on a copy of that local database instead of loading CARD for every sample
- Kleborate (`-p kpsc --trim_headers`) and SISTR (`--qc`) options moved from the module scripts to `ext.args` in `conf/modules.config` (same command lines)
- Every annotation and typing tool reads one prepared FASTA per sample (new local `SEQKIT_FILTER` step, SeqKit 2.13.0): contigs shorter than `--min_contig_length` (default 1000 bp) removed, uncompressed, contigs renamed `contig_1`, `contig_2`, ... with the original header kept as the description, published with a before/after `seqkit stats` table in `<sample>/assemblies/annotation_input/`. Before, Bakta alone filtered by length and renamed contigs, geNomad, MOB-suite and the species tools read Bakta's `.fna` (or the raw assembly without `--bakta_db`), and CheckM2, MLST, RGI and GTDB-Tk read the raw assembly. QUAST still reports the whole assembly. Outputs that change: CheckM2, MLST, RGI and GTDB-Tk no longer see contigs below 1000 bp; without `--bakta_db`, geNomad, MOB-suite and the species tools no longer see them either; Bakta keeps these names (`--keep-contig-headers`), so every tool reports the same `contig_N` IDs; tools that read the raw assembly before (RGI always, and without `--bakta_db` also geNomad, MOB-suite and the species tools) now report `contig_N` instead of the assembler's names (e.g. `NODE_1_length_…`); GTDB-Tk's `user_genome` is the sample name (was e.g. `<sample>.scaffolds` for Illumina). Kleborate, SCCmec and SISTR results no longer depend on whether `--bakta_db` is set: without it, short low-coverage contigs used to add SCCmec BLAST hits and change SISTR's closest cgMLST genome; the results now match the runs with Bakta
- `--min_contig_length` applies to every annotation tool, not only Bakta, and must be at least 1 (`0` used to be rejected only by Bakta). A sample with no contig that long is left out of annotation and typing with a warning instead of reaching the tools as an empty file

### `Fixed`

- RGI's version is recorded in `software_versions.yml` (it was empty: the container's busybox `grep` has no `-P`)
- Annotation tools no longer receive gzipped assemblies: without `--bakta_db`, Kleborate, sccmec, MOB-suite and geNomad got SPAdes' `.fa.gz` (or gzipped samplesheet contigs) directly. The assembly is decompressed once by `SEQKIT_FILTER`; SISTR's own gunzip step is removed
- `AMRFINDERPLUS_RUN` passes `--organism` only when the given AMRFinderPlus database lists it, and otherwise runs without it with a warning, instead of failing the sample (e.g. an older `--amrfinder_db`)
- `RGI_LOAD` fails when the protein, BLAST or DIAMOND database was not built (RGI ignores those build errors), instead of every `RGI_MAIN` rebuilding them
- _Salmonella_ samples now reach SISTR and get the AMRFinderPlus `Salmonella` organism: mlst 2.25.0 reports the scheme as `senterica_achtman_2`, which the MLST-to-organism map did not list. Their AMRFinderPlus results now include point mutations (`<sample>-mutations.tsv`, e.g. `gyrA_D87Y`) and apply the organism's gene suppression (e.g. `fieF` is no longer reported). The map was checked against `mlst --longlist` (2.25.0) and `amrfinder --list_organisms` (4.2.5): `vcholerae_2`, `diphtheria_3` and `helicobacter` added; `sepidermidis` removed, because AMRFinderPlus has no `Staphylococcus_epidermidis` organism and exited with an error for those samples
- Fixed implicit closure parameter warnings in map/branch operations
- Fixed deprecated `Channel` factory usage (now uses lowercase `channel`)
- Fixed dnaapler version command (was incorrectly calling medaka)
- Fixed AMRFinderPlus null organism error when MLST scheme not in taxa_map
- Fixed the misspelled _Burkholderia pseudomallei_ entry in the MLST-to-AMRFinderPlus organism map, and added verified mappings for seven more species
- Mash coverage estimation no longer crashes on paired-end input (per-file estimates are summed); MLST output is parsed defensively
- Samples no longer vanish silently when dnaapler produces no reoriented assembly (falls back to the unoriented assembly with a warning)
- SISTR's conditionally-produced outputs are optional, so runs no longer fail over missing novel-allele files
- `ConsolidateResults.py`: fixed silently-dropped kleborate/sccmec tables, added ONT support and defensive parsing
- Duplicate sample names in the samplesheet are rejected at launch (`uniqueEntries`); before, joins keyed on the sample could mis-join or drop samples
- Database paths are checked at launch (`exists` in the schema, `checkIfExists` for Kraken2, Bracken, RGI and MOB-suite)
- Cloud profiles: resources come from `conf/base.config` with retry scaling (Google Batch machine types no longer smaller than the requested memory), `google.location` replaces the ignored `google.region`, and the dead `google.lifeSciences` block and the AMI-specific AWS CLI default are gone
- MultiQC shows QUAST: it was given the `<sample>.tsv` link, which its QUAST module does not read (it reads `report.tsv`), so the report never had a QUAST section, and `--seq_type contig` runs without `--bakta_db` failed (MultiQC found no results and wrote no report). QUAST now labels each assembly with the sample ID (`--labels`), so the `Assembly` row of `<sample>/qc/quast/<sample>.tsv` is the sample ID instead of the FASTA file name, and samples whose FASTA files share a name no longer merge into one MultiQC row

### `Removed`

- Azure configuration and the nf-core community scaffolding (issue/PR templates, gitpod/devcontainer, lint workflows)
- Legacy `lib/` Groovy classes and the bundled jar; legacy samplesheet-check chain (`input_check.nf`, `check_samplesheet.py`)
- The no-op email/webhook notification stubs and their parameters/assets
- Dead parameters (`gtdbtk_mash_db`, `gtdbtk_min_perc_aa`, `gtdbtk_min_af`, `unicycler_*`, `antismash_*`, and friends); the schema now matches `nextflow.config` exactly
- Unused nf-core modules (antismash, macrel, and the stale kleborate/amrfinderplus duplicates)
- `accesory_scripts/` renamed to `accessory_scripts/`

### `Dependencies`

- Nextflow: `>=26.04.0` (minimum version; the version the pipeline is developed and CI-tested on)
- nf-schema: `2.2.0` (replacing nf-validation)

## v1.0dev - [date]

Initial release of gene2dis/mgap, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Support for ONT data (December 2023)

### `Fixed`

### `Dependencies`

- Bakta 1.8.2, MLST 2.23, MultiQC 1.18, AMRFinderPlus 3.11.18, geNomad 1.5.2 (December 2023)

### `Deprecated`
