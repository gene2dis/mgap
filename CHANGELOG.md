# gene2dis/mgap: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [Unreleased]

Major refactoring to align with nf-core best practices and modern Nextflow patterns.

### `Added`

- Migrated to `nf-schema@2.2.0` plugin for parameter validation (replacing deprecated `nf-validation`)
- Added `PIPELINE_INITIALISATION` and `PIPELINE_COMPLETION` subworkflows for standardized initialization
- Added cloud configuration profiles: `awsbatch`, `googlebatch`, `azurebatch`, `tower`
- Added Wave container and Fusion filesystem support for cloud execution
- Added nf-test framework with pipeline and module tests
- Added GitHub CI/CD workflows for linting, testing, and releases
- Added `meta.yml` and `environment.yml` to all local modules
- Added stub sections to key modules for fast testing
- Added pre-commit hooks configuration
- Added issue templates for bug reports and feature requests

### `Changed`

- Modernized `nextflow.config` with `resourceLimits`, updated profiles, and cloud parameters
- Refactored `main.nf` to use subworkflow-based initialization pattern
- Updated `workflows/mgap.nf` with explicit closure parameters and version collection
- Modernized `subworkflows/local/illumina.nf` and `ont.nf` with proper channel handling
- Updated `nextflow_schema.json` to JSON Schema draft 2020-12 with organized sections
- Updated documentation in `docs/usage.md` with cloud execution examples
- All local modules now use `${moduleDir}/environment.yml` syntax
- **Updated nf-core modules to latest versions:**
  - `nanoq`: Updated module with output_format parameter
  - `porechop/abi`: Updated module with custom_adapters parameter
  - `flye`: Updated to v2.9.5, outputs gzipped FASTA
  - `spades`: Updated to v4.1.0, outputs gzipped FASTA
  - `fastp`: Updated module with new input signature (adapter in tuple)
  - `quast`: Updated module with 3-input signature (consensus, reference, gff)
  - `bakta`: Updated module with 6-input signature (+regions, +hmms)
  - `checkm2/predict`: Updated module with database tuple input
  - `amrfinderplus/run`: Updated to v4.2.5
  - `genomad/endtoend`: Updated module
  - `kleborate`: Updated module (no species parameter)
  - `staphopiasccmec`: Updated module
  - `gtdbtk/classifywf`: Updated module with 3-input signature
- Changed default `gtdbtk_extension` from `fa` to `gz` to match gzipped SPAdes output

### `Fixed`

- Fixed implicit closure parameter warnings in map/branch operations
- Fixed deprecated `Channel` factory usage (now uses lowercase `channel`)
- Fixed dnaapler version command (was incorrectly calling medaka)
- Fixed AMRFinderPlus null organism error when MLST scheme not in taxa_map
- Fixed GTDB-Tk extension parameter not being passed to the module

### `Dependencies`

- Nextflow: `>=23.04.0` (minimum version)
- nf-schema: `2.2.0` (replacing nf-validation)

### `Deprecated`

- Legacy `lib/` Groovy files (NfcoreSchema, NfcoreTemplate, WorkflowMain, WorkflowMgap, Utils)
  - Functionality replaced by nf-schema plugin and subworkflows
  - Will be removed in a future release

## v1.0dev - [date]

Initial release of gene2dis/mgap, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
