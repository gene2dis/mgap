# gene2dis/mgap: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

gene2dis/mgap is a Nextflow pipeline for microbial genome analysis. It supports:

- **Illumina short-read assembly** (FastP → SPAdes)
- **Oxford Nanopore long-read assembly** (fastplong → Flye → Medaka → Dnaapler, or Autocycler consensus → Dnaapler)
- **Pre-assembled contigs** (direct annotation)

The pipeline performs quality control, assembly, annotation (Bakta), taxonomic classification (GTDB-Tk, optional), AMR detection (AMRFinderPlus), resistome analysis (RGI, optional), mobile element detection (geNomad), plasmid reconstruction (MOB-suite, optional), and species-specific analyses (Kleborate for _Klebsiella_, SCCmec typing for _S. aureus_, SISTR serotyping for _Salmonella_).

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1 or ONT reads. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". Leave empty for ONT single-end reads.                      |
| `fasta`   | Full path to FASTA file for pre-assembled contigs. Supports `.fasta`, `.fa`, `.fna` extensions (with optional `.gz` compression). Used only for `--seq_type contig` mode.              |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Automatic samplesheet generation

For convenience, a helper script is provided to automatically generate samplesheets from a directory of sequencing files:

```bash
python accesory_scripts/CreateSampleSheet.py <input_directory> <output_samplesheet.csv>
```

The script supports three types of sequencing data and can auto-detect the type. It automatically generates the correct column format based on the detected data type:

- **Illumina paired-end reads**: Automatically pairs R1/R2 files → outputs `sample,fastq_1,fastq_2` columns
- **Oxford Nanopore reads**: Single FASTQ files → outputs `sample,fastq_1` columns  
- **Pre-assembled contigs**: FASTA files → outputs `sample,fasta` columns

The script intelligently detects the data type by examining file extensions and naming patterns, ensuring the samplesheet format matches the pipeline's requirements for each sequencing mode.

#### Usage examples

```bash
# Auto-detect data type (recommended)
python accesory_scripts/CreateSampleSheet.py /path/to/fastq_files samplesheet.csv

# Explicitly specify Illumina paired-end data
python accesory_scripts/CreateSampleSheet.py /path/to/fastq_files samplesheet.csv --type illumina

# Oxford Nanopore long reads
python accesory_scripts/CreateSampleSheet.py /path/to/ont_reads samplesheet.csv --type ont

# Pre-assembled contigs
python accesory_scripts/CreateSampleSheet.py /path/to/contigs samplesheet.csv --type contig
```

#### Sample name extraction

The script intelligently extracts sample names from filenames:

- For Illumina data with standard naming (e.g., `Sample_S123_R1_001.fastq.gz`), it extracts `Sample`
- For simple R1/R2 naming (e.g., `Sample_R1.fastq.gz`), it extracts `Sample`
- For contigs and ONT data, it uses the base filename without extensions

#### Supported file formats

- **FASTQ files**: `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`
- **FASTA files**: `.fasta`, `.fa`, `.fna` (all optionally gzipped: `.fasta.gz`, `.fa.gz`, `.fna.gz`)

**Note:** The pipeline fully supports gzipped FASTA files for contig mode, allowing you to work with compressed assemblies directly without manual decompression.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
# Illumina short-read data
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --bakta_db /path/to/bakta_db \
    --checkm2_db /path/to/checkm2_db \
    -profile docker

# Oxford Nanopore long-read data
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type ont \
    --flye_mode nano-hq \
    -profile docker

# Pre-assembled contigs
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type contig \
    -profile docker

# With GTDB-Tk taxonomic classification (optional)
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --run_gtdbtk \
    --gtdbtk_db /path/to/gtdbtk_db \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

### Using a parameters file

Instead of passing every option on the command line, you can collect them in a YAML (or JSON) file and pass it to Nextflow with `-params-file`:

```bash
nextflow run gene2dis/mgap -params-file params.yaml -profile docker
```

A ready-to-edit template is provided at [`docs/params_template.yaml`](params_template.yaml). Copy it, fill in the required fields (`input`, `outdir`, `seq_type`) and any database paths you need, and leave everything else at its default:

```bash
cp docs/params_template.yaml my_run_params.yaml
# edit my_run_params.yaml to set input, outdir, seq_type, databases, etc.
nextflow run gene2dis/mgap -params-file my_run_params.yaml -profile docker
```

Any parameters you leave out of the YAML file fall back to the pipeline defaults defined in `nextflow.config`. Values passed on the command line (e.g. `--outdir results_v2`) override values in the params file.

### Autocycler Consensus Assembly (ONT)

For ONT data, the pipeline supports an alternative assembler: [Autocycler](https://github.com/rrwick/Autocycler) (v0.6.0+). Autocycler generates consensus long-read assemblies by running multiple assemblers on subsampled reads and combining the results, producing higher-quality assemblies than any single assembler alone.

**Basic usage:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type ont \
    --ont_assembler autocycler \
    -profile docker
```

**With custom assembler set (including slow assemblers):**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type ont \
    --ont_assembler autocycler \
    --autocycler_assemblers 'raven,miniasm,flye,metamdbg,necat,nextdenovo,canu,myloasm' \
    -profile docker
```

**With Plassembler (requires database):**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type ont \
    --ont_assembler autocycler \
    --autocycler_assemblers 'raven,miniasm,flye,metamdbg,necat,nextdenovo,plassembler' \
    --plassembler_db /path/to/plassembler_db \
    -profile docker
```

**Autocycler parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ont_assembler` | `flye` | ONT assembler: `flye` (Flye + Medaka) or `autocycler` (consensus multi-assembler) |
| `--autocycler_assemblers` | `raven,miniasm,flye,metamdbg,necat,nextdenovo` | Comma-separated list of assemblers for Autocycler |
| `--autocycler_read_type` | `ont_r10` | Read type: `ont_r9`, `ont_r10`, `pacbio_clr`, `pacbio_hifi` |
| `--autocycler_max_contigs` | `50` | Maximum number of contigs to retain per sample during compress and cluster steps |
| `--plassembler_db` | `null` | Path to Plassembler database. Plassembler is skipped if not provided |

**Supported assemblers:** `raven`, `miniasm`, `flye`, `metamdbg`, `necat`, `nextdenovo`, `canu`, `myloasm`, `plassembler`

> **Note:** The default assembler set excludes slow assemblers (Canu) and those requiring external databases (Plassembler). Add them explicitly via `--autocycler_assemblers` if needed.

> **Note:** When using `--ont_assembler autocycler`, the coverage adjustment step (Mash + Seqtk) is automatically skipped, as Autocycler handles its own read subsampling internally.

#### Autocycler Docker image

The Autocycler modules require a Docker image containing Autocycler and all supported assemblers. A pre-built image is available on Docker Hub and will be pulled automatically:

```
microds/autocycler:0.6.0
```

A Dockerfile is also provided in the repository if you need to build a custom image:

```bash
cd docker/autocycler-suite
docker build -t microds/autocycler:0.6.0 .
```

For Singularity users, the image is pulled automatically via `docker://microds/autocycler:0.6.0`. To build manually:

```bash
singularity build autocycler.sif docker://microds/autocycler:0.6.0
```

### Read QC Parameters

#### Illumina (FastP)

[FastP](https://github.com/OpenGene/fastp) performs quality trimming, adapter removal, and filtering of Illumina paired-end reads. The following parameters control its behaviour:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastp_qualified_quality` | `15` | Minimum base quality score (phred). Bases below this are considered unqualified |
| `--fastp_unqualified_percent_limit` | `40` | Maximum percentage of unqualified bases allowed per read (0–100) |
| `--fastp_cut_front_window_size` | `4` | Sliding window size for 5′ end trimming |
| `--fastp_cut_front_mean_quality` | `20` | Mean quality threshold for 5′ end sliding window trimming |
| `--fastp_cut_right_window_size` | `4` | Sliding window size for 3′ end trimming |
| `--fastp_cut_right_mean_quality` | `20` | Mean quality threshold for 3′ end sliding window trimming |
| `--fastp_reads_minlength` | `50` | Minimum read length after trimming. Shorter reads are discarded |
| `--fastp_n_base_limit` | `5` | Maximum number of N bases allowed per read |

#### ONT (fastplong)

[fastplong](https://github.com/OpenGene/fastplong) performs quality filtering and adapter trimming of Oxford Nanopore long reads.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fastplong_qualified_quality` | `15` | Minimum base quality score (phred). Bases below this are considered unqualified |
| `--fastplong_unqualified_percent_limit` | `40` | Maximum percentage of unqualified bases allowed per read (0–100) |
| `--fastplong_min_read_length` | `1000` | Minimum read length after filtering. Shorter reads are discarded |

### Annotation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_contig_length` | `1000` | Minimum contig length (bp) for Bakta annotation. Contigs shorter than this are skipped |

### Dnaapler Contig Reorientation (ONT)

[Dnaapler](https://github.com/gbouras13/dnaapler) (v1.3.0) reorients assembled microbial sequences so that each contig starts at a consistent location (e.g., at the *dnaA*, *repA*, or *terL* gene). This is an optional step enabled by default (`--run_dnaapler true`).

The behavior differs depending on the assembler:

- **Flye mode:** Dnaapler runs on the polished FASTA output from Medaka. All contigs are candidates for reorientation.
- **Autocycler mode:** Dnaapler runs on the consensus GFA from `autocycler combine`. Only circular contigs are reoriented (GFA-aware mode). The reoriented GFA is then converted back to FASTA via `autocycler gfa2fasta`. This follows the [recommended Autocycler workflow](https://github.com/rrwick/Autocycler/wiki/Reorienting-with-Dnaapler).

**Dnaapler parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_dnaapler` | `true` | Enable contig reorientation with Dnaapler |

**To disable Dnaapler:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type ont \
    --run_dnaapler false \
    -profile docker
```

### GTDB-Tk Taxonomic Classification

GTDB-Tk provides taxonomic classification using the Genome Taxonomy Database (GTDB). This is an optional step that can be enabled with the `--run_gtdbtk` flag.

**Batch Processing:** The pipeline processes all genomes together in a single GTDB-Tk run, which significantly reduces runtime compared to processing genomes individually. Results are consolidated into a single set of output files in the `gtdbtk/` directory.

**Basic usage:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --run_gtdbtk \
    --gtdbtk_db /path/to/gtdbtk_db \
    -profile docker
```

**Advanced GTDB-Tk parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_gtdbtk` | `false` | Enable GTDB-Tk taxonomic classification |
| `--gtdbtk_db` | `null` | Path to GTDB-Tk reference database (required if enabled) |
| `--gtdbtk_mash_db` | `null` | Path to Mash database for ANI screening (optional, skips ANI if not provided) |
| `--gtdbtk_extension` | `gz` | File extension for genome files (see table below) |
| `--gtdbtk_min_perc_aa` | `10` | Minimum percentage of amino acids in MSA |
| `--gtdbtk_min_af` | `0.65` | Minimum alignment fraction |
| `--gtdbtk_pplacer_scratch` | `true` | Use scratch directory for pplacer to reduce memory usage |

**Important: File extension by sequencing type**

The `--gtdbtk_extension` parameter must match the file extension of your assembled genomes:

| `--seq_type` | Assembler | Output extension | `--gtdbtk_extension` |
|--------------|-----------|------------------|----------------------|
| `illumina` | SPAdes | `.scaffolds.fa.gz` | `gz` (default) |
| `ont` | Medaka | `.fasta` | `fasta` |
| `contig` | N/A | varies | match your input files |

**Example with ONT data (requires extension override):**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type ont \
    --run_gtdbtk \
    --gtdbtk_db /path/to/gtdbtk_db \
    --gtdbtk_extension fasta \
    -profile docker
```

**Example with Illumina data (uses default extension):**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --run_gtdbtk \
    --gtdbtk_db /path/to/gtdbtk_db \
    -profile docker
```

**Note:** The GTDB-Tk database is large (~70GB) and must be downloaded separately. See the [GTDB-Tk documentation](https://ecogenomics.github.io/GTDBTk/installing/index.html) for database installation instructions.

### RGI Antimicrobial Resistance Gene Prediction

RGI (Resistance Gene Identifier) predicts resistomes from genome assemblies using the Comprehensive Antibiotic Resistance Database (CARD). This is an optional step that can be enabled with the `--run_rgi` flag.

**Basic usage:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --run_rgi \
    --rgi_db /path/to/card_db \
    -profile docker
```

**RGI parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_rgi` | `false` | Enable RGI antimicrobial resistance gene prediction |
| `--rgi_db` | `null` | Path to RGI/CARD database directory (required if enabled) |
| `--rgi_include_loose` | `false` | Include loose hits (below detection model cut-off) |
| `--rgi_include_nudge` | `false` | Nudge loose hits to strict for partial gene sequences |
| `--rgi_low_quality` | `false` | Use low quality mode for short contigs to predict partial genes |
| `--rgi_alignment_tool` | `DIAMOND` | Alignment tool: `BLAST` or `DIAMOND` |

**Example with loose hits and low quality mode:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --run_rgi \
    --rgi_db /path/to/card_db \
    --rgi_include_loose \
    --rgi_low_quality \
    -profile docker
```

**Example using BLAST instead of DIAMOND:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --run_rgi \
    --rgi_db /path/to/card_db \
    --rgi_alignment_tool BLAST \
    -profile docker
```

**Database setup:**

The CARD database must be downloaded and set up before running RGI. See the [RGI documentation](https://github.com/arpcard/rgi#card-database) for database installation instructions:

```bash
# Download CARD database
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json

# Load database (creates localDB directory)
rgi load --card_json card.json --local
```

Then provide the path to the database directory (containing the loaded database files) via `--rgi_db`.

### MOB-suite Plasmid Detection

[MOB-suite](https://github.com/phac-nml/mob-suite) performs plasmid identification, typing, and reconstruction from bacterial assemblies. This is an optional step that can be enabled with the `--run_mobsuite` flag. It consumes the Bakta-annotated nucleotide sequences (`.fna`).

**Basic usage:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --run_mobsuite \
    -profile docker
```

**MOB-suite parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_mobsuite` | `false` | Enable MOB-suite plasmid detection and reconstruction |
| `--mobsuite_db` | `null` | Path to a pre-built MOB-suite database directory (optional) |

**Database handling:**

MOB-suite ships its reference databases with the container/conda environment and can also download them on first execution. For air-gapped or repeated runs, pre-download the database and pass it via `--mobsuite_db` — the pipeline forwards it to `mob_recon --database_directory` through `ext.args`:

```bash
# One-time database initialisation (inside the mob_suite conda env or container)
mob_init -d /path/to/mobsuite_db

# Use the pre-built database
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --run_mobsuite \
    --mobsuite_db /path/to/mobsuite_db \
    -profile docker
```

### MLST Typing

[MLST](https://github.com/tseemann/mlst) (Torsten Seemann's classic MLST) scans assemblies against PubMLST typing schemes and is run automatically on every sample. The scheme it reports also drives downstream logic in the pipeline (for example, a `senterica_achtman_2` call triggers SISTR).

The tool ships with its own bundled BLAST database and PubMLST data directory inside the container, so no databases are required by default. If you need to run against an updated or custom set of schemes, you can override either one:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mlst_blastdb` | `null` | Path to an alternative MLST BLAST database fasta file. Forwarded to `mlst --blastdb`. |
| `--mlst_datadir` | `null` | Path to an alternative MLST PubMLST data directory. Forwarded to `mlst --datadir`. |

Both parameters are optional and independent — set either, both, or neither. When unset, `mlst` uses the database and data directory bundled with its container.

**Example — using a custom MLST database:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --mlst_blastdb /path/to/mlst/blast/mlst.fa \
    --mlst_datadir /path/to/mlst/pubmlst \
    -profile docker
```

See the [mlst documentation](https://github.com/tseemann/mlst#updating-the-database) for instructions on building or updating these databases.

### SISTR Salmonella Serotyping

[SISTR](https://github.com/phac-nml/sistr_cmd) (Salmonella In Silico Typing Resource) predicts serotypes and cgMLST subtypes directly from assemblies. The pipeline runs SISTR **v1.1.3** automatically whenever MLST detects a _Salmonella enterica_ assembly (MLST scheme `senterica_achtman_2`). No additional flags or databases are required — SISTR bundles its own reference data in the container.

SISTR is skipped for any other species and adds no runtime cost to non-Salmonella samples.

### Kraken2 Contamination Detection

[Kraken2](https://ccb.jhu.edu/software/kraken2/) classifies reads against a reference database to detect contamination. For Illumina data, [Bracken](https://ccb.jhu.edu/software/bracken/) is also run to refine abundance estimates. For ONT data, Kraken2 runs alone (no Bracken step).

This step is **enabled by default** but requires a database (`--kraken2db`). If `--kraken2db` is not provided, the step is silently skipped even if `--run_kraken2` is `true`.

**Kraken2/Bracken parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_kraken2` | `true` | Enable Kraken2 contamination detection |
| `--kraken2db` | `null` | Path to Kraken2 database directory (required to run) |
| `--brackendb` | `null` | Path to Bracken database directory (Illumina only; uses same DB as Kraken2 if built with Bracken) |

**Basic usage with contamination detection:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --kraken2db /path/to/kraken2_db \
    --brackendb /path/to/kraken2_db \
    -profile docker
```

**Disable contamination detection:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --run_kraken2 false \
    -profile docker
```

> **Note:** Bracken runs only for `--seq_type illumina`. ONT mode runs Kraken2 only (no Bracken re-estimation).

### Coverage Adjustment

The pipeline estimates genome coverage using [Mash](https://mash.readthedocs.io/en/latest/) and subsamples reads with [Seqtk](https://github.com/lh3/seqtk) if coverage exceeds the target threshold. This prevents over-assembly and reduces runtime for high-coverage datasets.

This step is **enabled by default** and applies to both Illumina (after FastP) and ONT Flye mode (after fastplong). It is automatically **skipped** for `--ont_assembler autocycler`, as Autocycler handles its own read subsampling internally.

**Coverage adjustment parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--adjust_coverage` | `true` | Enable coverage estimation and subsampling |
| `--max_coverage` | `110` | Target maximum coverage (x). Reads are subsampled if coverage exceeds this value |

**Disable coverage adjustment:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    --adjust_coverage false \
    -profile docker
```

**Use a lower coverage target:**

```bash
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type ont \
    --max_coverage 60 \
    -profile docker
```

### Cloud execution

The pipeline supports execution on AWS, Google Cloud, and Azure:

```bash
# AWS Batch
nextflow run gene2dis/mgap \
    --input s3://bucket/samplesheet.csv \
    --outdir s3://bucket/results \
    --seq_type illumina \
    --awsqueue my-batch-queue \
    -profile awsbatch,docker

# Google Cloud Batch
nextflow run gene2dis/mgap \
    --input gs://bucket/samplesheet.csv \
    --outdir gs://bucket/results \
    --seq_type illumina \
    --gcp_project my-project \
    -profile googlebatch,docker

# Seqera Platform (Tower)
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    -profile tower,docker
```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull gene2dis/mgap
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [gene2dis/mgap releases page](https://github.com/gene2dis/mgap/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of rnaseq](https://nf-co.re/rnaseq/3.9/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Legacy / Inactive Parameters

The following parameters exist in `nextflow.config` but correspond to tools that are currently **not active** in the pipeline. They are reserved for future use and have no effect on pipeline runs.

| Parameter group | Tools | Status |
|-----------------|-------|--------|
| `unicycler_*` (`unicycler_min_fasta_length`, `unicycler_mode`) | Unicycler | Not used; Unicycler is not part of the current assembly workflow |
| `antismash_*` (11 params: `antismash_db`, `antismash_install`, `antismash_cbgeneral`, etc.) | antiSMASH | Module present but commented out; not executed |

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
