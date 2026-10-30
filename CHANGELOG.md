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
- `--flye_asm_coverage` (off by default) for very deep ONT samples: Flye builds its first assembly step from the longest reads up to that coverage (`--asm-coverage`), with the genome size from `--genome_size` or, if unset, a per-sample estimate by the new local `LRGE` module (lrge 1.0.0, `quay.io/biocontainers/lrge:1.0.0--hf1b6044_0`; published in `<sample>/assemblies/lrge/`). In the benchmark (`docs/benchmarks/coverage_2026-10-24.md`) `--asm-coverage 40` made Flye 14–26 % faster but lost a sample's 3.9–5.2 kb plasmids, so it stays opt-in; a zero `--genome_size` fails at launch
- `--brackendb` is now actually used by Bracken (falls back to the Kraken2 DB)
- Per-sample cross-validation of samplesheet columns against `--seq_type` with clear error messages
- Launch-time warnings for annotation steps skipped because their database was not provided
- Stub blocks in every local module (fast `-stub` runs cover all three input modes)
- Sample QC gate (issues #14, #3, #2; D-7): samples with fewer than `--min_bases` (default 50,000,000) bases after fastp/fastplong skip Kraken2/Bracken, assembly and annotation; SPAdes, Flye and Medaka retry an out-of-memory failure once and then ignore it, so one failed sample no longer stops the run; QUAST failures are a warning. `pipeline_info/sample_status.tsv` lists every sample as pass, skipped or failed with the stage, reason and read and assembly numbers (reads, bases, mean and, for ONT, approximate N50 read length from the fastp/fastplong reports; contigs, length and N50 from `seqkit stats`), and the completion log prints a summary line (with how many samples with an assembly QUAST failed for). The run exits 0 unless every sample that went into assembly failed
- Run report (D-8): the local `MGAP_REPORT` module (`bin/mgap_report.py`, Seqera Containers image with Python 3.12, pandas 3.0.6, Plotly 7.1.0, PyYAML 6.0.3 and openpyxl 3.1.5) writes `report/mgap_report.html` (one offline page: run overview with database versions, read QC, species in the reads, assembly, CheckM2, typing, AMR heatmaps, mobile elements, per-sample cards, versions and citations), `report/mgap_summary.tsv` and `report/mgap_summary.xlsx` (the summary and one sheet per tool). Above 100 samples the per-sample bar charts become distributions with the outliers named; sample cards start collapsed, with a search box and a list that puts samples with warnings first. `--skip_report` turns it off; `accessory_scripts/build_report_manifest.py` builds the report's inputs from an existing results folder. `pipeline_info/software_versions.yml` now includes the report's own libraries
- CI job running the `bin/` script tests (`pytest tests/bin`) with the report's package versions
- Plasmid consensus (issue #17, D-22): PlasmidFinder 2.1.6 (nf-core `plasmidfinder`, `quay.io/biocontainers/plasmidfinder:2.1.6--py310hdfd78af_1`, patched with an optional database input) runs by default on the annotation FASTA (`--run_plasmidfinder`, `--plasmidfinder_db`; without it the database bundled in the container, from 2019); with `--run_mobsuite`, MOB-typer 3.1.9 (nf-core `mobsuite/typer`, `--multi`) types each circular contig of ONT assemblies (Flye `assembly_info.txt`, Autocycler `circular=true`) and contig-mode input (header tags, else every contig), where `mob_recon` reconstructs no plasmids; the local `PLASMID_CONSENSUS` module (`bin/plasmid_consensus.py`) writes `<sample>/annotation/plasmid_consensus/<sample>.plasmid_consensus.tsv`, one row per contig with the geNomad, MOB-suite and PlasmidFinder evidence, circularity, mobility and a call (plasmid with 2 or more sources, candidate with 1, check on contigs of 1 Mb or more, chromosome, no call; rules in docs/output.md). The run report gets a Plasmid consensus section, summary columns and an XLSX sheet
- Berokka (issue #17, D-23): with Flye (default ONT assembler), the new local `BEROKKA` module (Berokka 0.2.3, `quay.io/biocontainers/berokka:0.2.3--hdfd78af_3`) trims the end overlap of the contigs Flye marks circular before Medaka (`--run_berokka`, default true; linear contigs unchanged). Flye can assemble a small plasmid as tandem copies: on the SCL10571 fixture a 5,216 bp plasmid came out as a 10,432 bp contig and a 4,163 bp plasmid as 12,428 bp, so their genes were annotated two and three times; with Berokka they are single copies, as in the Autocycler consensus. Flye-mode assemblies change only where a circular contig's start repeats at its end; `<sample>/assemblies/berokka/<sample>.berokka.tsv` lists what was trimmed. Berokka's control-sequence filter is off (`--filter ''`), and a sample whose Berokka task fails is polished from the untrimmed Flye assembly with a warning
- ISEScan (issue #17, D-9): the new local `ISESCAN` module (ISEScan 1.7.3, `quay.io/biocontainers/isescan:1.7.3--hab16a5f_1`) finds insertion sequences in the annotation FASTA when `--run_isescan` is set (default false; no database). Results in `<sample>/annotation/isescan/` (`<sample>.isescan.tsv` one row per element with family and complete/partial type, header only when none is found; `.sum`, `.gff`, `.is.fna`, `.orf.fna`, `.orf.faa`). The run report gets an "Insertion sequences" section (elements per family and sample; complete, partial, IS length, elements on plasmid and candidate plasmid contigs), a card line, `isescan_elements`, `isescan_complete` and `isescan_families` summary columns and an ISEScan workbook sheet. About 20 min per genome with 6 CPUs and under 2 GB (SCL2535 Illumina: 31 elements, 13 complete; SCL10571 ONT: 84, 66 complete). A failed ISEScan task is retried once after an out-of-memory exit (with more memory, the same 6 CPUs), then ignored: that sample has no ISEScan results, `ISEScan failed` is noted in `sample_status.tsv`, the report and the completion log, and the run goes on. An ISEScan run that writes no table without reaching its last step fails rather than report no IS element. Other outputs unchanged
- Species check for routing (B11, S1–S3; D-13, D-15): with the new `--species_refs`, each annotation FASTA is compared with skani 0.4.0 (nf-core `skani/search`, patched) against a curated set of 53 NCBI type or reference genomes (one per AMRFinderPlus organism plus close relatives such as _K. variicola_, _E. fergusonii_, _S. argenteus_, _S. bongori_, _E. hormaechei_, _C. lari_, _B. bronchiseptica_; `assets/species_refs.tsv`), built once by `accessory_scripts/build_species_refs.py` (NCBI datasets, never by the pipeline; `docs/usage.md`). The new local `SPECIES_CHECK` module (`bin/species_check.py`) combines the MLST scheme, the reference call (ANI ≥ 95 %, ≥ 50 % of the assembly aligned), the top species in the reads (Bracken, Kraken2 on ONT) and the samplesheet species into one decision per sample, which sets the AMRFinderPlus organism and the species tools. When the reference species is not typed by MLST's scheme, the reference call wins, mlst runs again with its scheme (`<sample>.<scheme>.tsv`) and the log warns. `pipeline_info/species_check.tsv` lists every sample's sources, organism, decision and reason; skani hits go to `<sample>/annotation/skani/`. Without `--species_refs` routing follows MLST (launch warning). Two references within 0.5 ANI points are `ambiguous` when their AMRFinderPlus organisms differ, or when neither has one (_B. bronchiseptica_ and _B. parapertussis_). mlst is not run again when `--mlst_datadir`, or without it the container's database (`assets/mlst_schemes.txt`), lacks the scheme (warning and reason in `species_check.tsv`; the report keeps the automatic result), and a _K. aerogenes_ call keeps the `Klebsiella_pneumoniae` AMRFinderPlus organism but does not go to Kleborate (KpSC preset)
- Optional samplesheet `species` column (S4, D-15b): recorded in `species_check.tsv`, flagged when the reference call or the MLST genus differs (genus and epithet compared; subspecies, strain and `sp.` ignored), and used for routing only when neither gives an AMRFinderPlus organism and MLST's scheme does not type another genus. The run report's typing table shows the species used, the closest reference with its ANI and the species check decision (flagged ones first), sample cards and warnings show flagged checks, and `mgap_summary.tsv` gains `species_check_decision`, `skani_species`, `skani_ani` and `amrfinder_organism` (`mlst_organism` is now the organism of the MLST scheme alone)

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
- The Mash/Seqtk coverage-adjustment logic is in the `COVERAGE_ADJUST` subworkflow (Illumina only) and tolerates unparseable coverage estimates
- Illumina coverage is now estimated as fastp's output bases ÷ Mash's genome size, with `mash sketch -k 32 -m 3` (was Mash's own coverage line with its defaults, `-k 21` and no k-mer copy filter, which counted sequencing-error k-mers as genome: 5–21 Mb genome sizes and 0.17–0.57× the real coverage, so subsampling never triggered, not even at 244x). The new estimate is within 0.90–1.09× of the coverage over the assembly (13 samples, `docs/benchmarks/coverage_2026-10-24.md`). `--max_coverage` is real base coverage and its default rises from 110 to 150 (in a depth sweep assemblies stopped improving between 110x and 150x). Outputs change for samples above 150x: their reads are subsampled to ~150x before SPAdes. `MASH_SKETCH` runs as `process_single` (Mash ignores `-p` with `-r`), publishes `<sample>.genome_size` instead of `<sample>.mash_coverage`, and no longer re-emits the reads, so `-resume` reuses SPAdes when only the coverage step changes
- ONT reads are no longer subsampled before Flye: the Flye path does not run `COVERAGE_ADJUST`, so `--adjust_coverage` and `--max_coverage` apply to Illumina only. Mash cannot estimate ONT coverage (genome sizes of 30–115 Mb and 12–15x for samples at 81x and 251x in the benchmark), so ONT reads were in practice never reduced; outputs do not change. A `--max_coverage` other than the default with `--seq_type ont` prints a launch warning
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
- Autocycler: compress, cluster, trim/resolve and combine run as one task per sample (`AUTOCYCLER_CONSENSUS`, replacing `AUTOCYCLER_COMPRESS`, `AUTOCYCLER_CLUSTER`, `AUTOCYCLER_TRIM_RESOLVE` and `AUTOCYCLER_COMBINE`), so the `autocycler_out` directory is built in one place and published once to `<sample>/assemblies/autocycler/autocycler_out/` (it now also holds `consensus_assembly.{fasta,gfa,yaml}`); `trim` and `combine` get the task's threads
- Autocycler image `microds/autocycler:0.8.0` (was `0.6.0`, which despite its tag contained Autocycler 0.5.2), built from `docker/autocycler-suite/Dockerfile` with every package and the base image pinned; `modules/local/autocycler/environment.yml` pins the same versions, and the image is defined once, in `conf/modules.config`. Version changes: Autocycler 0.5.2 → 0.8.0, metaMDBG 1.3 → 1.4, minimap2 2.28 → 2.31, minipolish 0.2.0 → 0.2.1, myloasm 0.4.0 → 0.7.0, Plassembler 1.8.2 → 1.8.5 (Dnaapler 1.3.0 → 1.4.0 inside the image); Canu 2.3, Flye 2.9.6, miniasm 0.3, NECAT, NextDenovo 2.5.2, NextPolish 1.4.1, Racon 1.5.0 and Raven 1.8.3 unchanged
- Autocycler resources from a real run (5.9 Mb genome, ~100x ONT): `AUTOCYCLER_ASSEMBLY` reserves 6 CPUs and 36 GB (72 GB on the out-of-memory retry, same CPUs; 16 h per attempt as before) instead of `process_high` (12 CPUs, 72 GB) (peak 11.6 GB, NECAT; 2–4 cores busy), so twice as many assemblers run at once; `AUTOCYCLER_GENOME_SIZE`, which runs Raven on all reads, keeps 6 CPUs and retries an out-of-memory exit twice (36, 72, then 108 GB; 12 h per attempt) (#17)
- `--autocycler_read_type` reaches the assembler helper through `ext.args` of `AUTOCYCLER_ASSEMBLY` (no `params` read inside the module)
- Bracken gets a read length (`-r`) per sample: the longest `database<N>mers.kmer_distrib` in the Bracken database at or below the sample's sequenced read length (fastp `read1_before_filtering.total_cycles`), or the shortest with a warning when none fits; new `--bracken_read_length` sets it for all samples. Before, no `-r` was passed, so Bracken used its default of 100 for 150 and 300 bp reads. Bracken fractions change slightly (at most 0.49 points in the benchmark; the top species did not change)
- Kraken2 runs with `--memory-mapping` (identical reports; per-sample time 147–216 s → 46–148 s with the database in the page cache, 438 s → 178 s from disk; `docs/benchmarks/coverage_2026-10-24.md`). The memory request stays `process_high` (72 GB). On the AWS Batch and Google Batch executors (`mgap_awsbatch`, `mgap_googlebatch`), where each task starts on a fresh machine, it stays off
- Medaka now polishes with the model named in the reads (dorado `RG:Z:` tag or MinKNOW `basecall_model_version_id=`) when `--medaka_model` is unset; Medaka 1.11.3 could not read either and used `r1041_e82_400bps_sup_v4.3.0` for every sample. When the reads name no model it knows, Medaka 2.2.2 uses `r1041_e82_400bps_sup_v5.2.0` (set `--medaka_model` for such reads). Polished ONT assemblies change where the new model differs: on SCL10571 (`sup@v5.0.0` reads) at 35 sites, 29 of them in one 22 kb repeat contig, with the same MLST and AMRFinderPlus calls. New `--medaka_bacteria` (default false) polishes with Medaka's bacterial methylation model (R10.4.1 400 bps hac/sup v4.2+ reads; for other reads the sample fails at polishing) (`docs/benchmarks/tool_updates_2026-10-29.md`)
- `bordetella_3` and `campylobacter_nonjejuni*` no longer give an AMRFinderPlus organism from the MLST scheme alone (D-14): they type species the organisms are not curated for (_B. parapertussis_, _B. bronchiseptica_; _C. lari_, _C. upsaliensis_ …). With `--species_refs`, a reference call of _B. pertussis_ or _C. jejuni_/_C. coli_ gives them `Bordetella_pertussis` or `Campylobacter`. Samples typed by these schemes lose point-mutation calls unless the reference call confirms the species
- Species tools and the AMRFinderPlus organism follow the species check instead of the MLST scheme alone; with `--species_refs`, a _K. pneumoniae_ ST11 that mlst typed `ecoli_achtman_4` ST14464 (SCL2461) gets `Klebsiella_pneumoniae` and Kleborate (ST11); the mlst re-run with the `klebsiella` scheme reports no ST for it (two phoE hits, `phoE(1,1)`), so the report's MLST column shows `klebsiella` `-` instead of `ecoli_achtman_4` ST14464

### `Fixed`

- A run no longer stops when one sample has too few reads to assemble (SPAdes exit 21, "Invalid kmer coverage histogram", issue #3) or when QUAST finds no contig of at least 500 bp in a highly fragmented assembly (exit 4, issue #2); with QUAST failed for every sample in contig mode without Bakta, MultiQC (which then has no report to read and fails) is skipped with a warning
- Medaka's out-of-memory kills are retried: `medaka_consensus` exits 1 when a step is killed, so the module now exits 137 when the task's cgroup counted an out-of-memory kill
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
- MLST output is parsed defensively
- Samples no longer vanish silently when dnaapler produces no reoriented assembly (falls back to the unoriented assembly with a warning)
- SISTR's conditionally-produced outputs are optional, so runs no longer fail over missing novel-allele files
- `ConsolidateResults.py`: fixed silently-dropped kleborate/sccmec tables, added ONT support and defensive parsing
- `ConsolidateResults.py` reads Kleborate 3's `klebsiella_pneumo_complex_output.txt` (Kleborate 3 names its output after the `-p kpsc` preset); it looked for `<sample>.results.txt`, which Kleborate 3 does not write, and stopped with `FileNotFoundError` on any run with a Kleborate sample. It also reads CheckM2's `<sample>_checkm2_report.tsv` and geNomad's summaries from `<sample>_summary/`; the old paths matched nothing, so the CheckM2 and geNomad tables were silently left out
- Duplicate sample names in the samplesheet are rejected at launch (`uniqueEntries`); before, joins keyed on the sample could mis-join or drop samples
- Database paths are checked at launch (`exists` in the schema, `checkIfExists` for Kraken2, Bracken, RGI and MOB-suite)
- Cloud profiles: resources come from `conf/base.config` with retry scaling (Google Batch machine types no longer smaller than the requested memory), `google.location` replaces the ignored `google.region`, and the dead `google.lifeSciences` block and the AMI-specific AWS CLI default are gone
- MultiQC shows QUAST: it was given the `<sample>.tsv` link, which its QUAST module does not read (it reads `report.tsv`), so the report never had a QUAST section, and `--seq_type contig` runs without `--bakta_db` failed (MultiQC found no results and wrote no report). QUAST now labels each assembly with the sample ID (`--labels`), so the `Assembly` row of `<sample>/qc/quast/<sample>.tsv` is the sample ID instead of the FASTA file name, and samples whose FASTA files share a name no longer merge into one MultiQC row
- Autocycler assemblers that run out of memory are retried once with more memory: `autocycler helper` reports a killed assembler only in its log (exit 0, no assembly; `genome_size` exits 1), so the out-of-memory retry never triggered and the assembler was silently missing from the consensus. The modules now exit 137 on that message; `AUTOCYCLER_ASSEMBLY` retries an out-of-memory exit once and then ignores the task, and `AUTOCYCLER_GENOME_SIZE` retries like any other process
- NextDenovo contributes to the Autocycler consensus: in the old image NextPolish could not load `calgs.so`, so every NextDenovo task failed, and `autocycler helper` exits 0 when an assembler fails, so it went unnoticed. An assembler task that ends without an assembly now fails (exit 1), which Nextflow reports as an ignored error; Plassembler on a genome without plasmids is not a failure (no assembly, exit 0). The helper's stderr reaches `.command.err` as it is written, so it is kept when a task is killed
- Singularity/Apptainer can pull the Autocycler image: it is pinned by digest alone (`docker.io/microds/autocycler@sha256:…`), as Apptainer rejects a reference with both a tag and a digest
- `--plassembler_db` is staged into the Plassembler tasks; before, only its host path was exported, which a container cannot see
- A sample whose Autocycler clustering leaves no QC-pass cluster gets a warning and is left out of downstream analysis; before, the trim/resolve step looped over the unmatched glob and failed the run
- Autocycler's `versions.yml` entry is the version number (`0.8.0`); it was `Autocycler 0.5.2`
- Flye's assembly and graph are published (`<sample>/assemblies/flye/*.assembly.fasta.gz`, `*.assembly_graph.gfa.gz`): nf-core flye gzips them, and the publish pattern had no `.gz`, so only the log, assembly info and parameters were published
- On/off parameters can be switched off on the command line: Nextflow 26.04 passes `--adjust_coverage false` (and `--run_kraken2 false`, `--run_dnaapler false`, `--validate_params false`, `--rgi_include_loose false`, ...) as the string `"false"`, which counted as true, so the step still ran. Boolean parameters are now read as booleans wherever they are used

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
- Autocycler suite image: `microds/autocycler:0.8.0` (digest-pinned; see `docs/usage.md` for the full version list)
- fastplong: `0.3.0` → `0.7.1` (`biocontainers/fastplong:0.7.1`, bioconda build `h43da1c4_0`; biocontainers publishes 0.7.1 only under this tag). nf-core module updated (0.4.1 upstream, now reporting its version through a topic tuple) and patched to 0.7.1 (`modules/nf-core/fastplong/fastplong.diff`). Filter thresholds unchanged; in the benchmark 0.7.1 kept the same bases within 0.1 point on current R10.4.1 sup data and 0.3 points fewer on an older R10.4 sample
- AMRFinderPlus: `4.2.5` → `4.2.7` (`biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0`), a bug-fix release (a crash when an HMM hit has no BLAST alignment, ncbi/amr#180). Outputs unchanged: all result and mutation tables of the 10-sample Illumina benchmark set were byte-identical to 4.2.5 with database 2025-12-03.1 (`docs/benchmarks/tool_updates_2026-10-29.md`)
- Dnaapler: `1.3.0` → `1.4.0` (`biocontainers/dnaapler:1.4.0--pyhdfd78af_0`, the version the Autocycler image already has): reads compressed input and fixes the start coordinate when the `largest`/`mystery`/`nearest` fallbacks pick a CDS on the negative strand. The reoriented SCL10571 assembly was byte-identical to 1.3.0's
- Kleborate container: `3.2.4--pyhdfd78af_0` → `3.2.4--pyhdfd78af_1`, the same Kleborate rebuilt with Kaptive 3.2.2 (was 3.1.0), minimap2 2.31 and Python 3.14. With `-p kpsc` every output field was identical on three K. pneumoniae samples
- RGI: `6.0.3` → `6.0.4` (`biocontainers/rgi:6.0.4--pyh05cac1d_1`). 6.0.3 cannot read some CARD 4.0.1 mutation models (`invalid literal for int()` in its log) and skips them, so it reports no point mutations; 6.0.4 keeps every 6.0.3 hit and adds them: on the 10-sample benchmark set 25 new `protein variant model` rows, most also called by AMRFinderPlus (gyrA, parC, murA, uhpT). Not the newest 6.0.8: from 6.0.5 on RGI drops Perfect hits such as CTX-M-65 and AAC(6')-Ib-cr5 on short contigs (`docs/benchmarks/tool_updates_2026-10-29.md`)
- GTDB-Tk: `2.5.2` → `2.8.0` (`biocontainers/gtdbtk:2.8.0--pyhdfd78af_0`; nf-core module patched, `gtdbtk-classifywf.diff`). 2.5.2 reads only GTDB r220/r226; r232 needs GTDB-Tk 2.7.0 or later, and 2.8.0 no longer reads r220/r226: `--gtdbtk_db` must be the r232 package. 2.8.0 classified the 10-sample benchmark set on r232 with the module's arguments unchanged (1 h 20 min, 177 GB peak memory of the 200 GB requested). Classifications can differ from older GTDB-Tk/GTDB releases (new taxonomy; 2.8.0 assigns the closest reference with AF ≥ `--min_af` before checking its ANI radius)
- Medaka: `1.11.3` → `2.2.2` in the new image `docker.io/microds/medaka` (tag 2.2.2, pinned by digest; built from `docker/medaka/Dockerfile`: the `medaka:2.2.2--py312h3050eb1_0` biocontainer plus the model files Medaka can choose automatically and its bacterial model, which the biocontainer would otherwise download at run time and fail to store inside a task). About 2.5× faster (15.7 → 6.2 min on the ~100x SCL10571 fixture) with less memory; with the same model the consensus is byte-identical to 1.11.3's
- skani: `0.4.0` (new; `quay.io/biocontainers/skani:0.4.0--hec9b1f2_0`, nf-core `skani/search` patched from 0.2.2, `skani-search.diff`; 0.4 databases cannot be read by older skani)

## v1.0dev - [date]

Initial release of gene2dis/mgap, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Support for ONT data (December 2023)

### `Fixed`

### `Dependencies`

- Bakta 1.8.2, MLST 2.23, MultiQC 1.18, AMRFinderPlus 3.11.18, geNomad 1.5.2 (December 2023)

### `Deprecated`
