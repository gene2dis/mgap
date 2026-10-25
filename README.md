## Introduction

**Microbial Genome Analysis Pipeline**

**MGAP** is a bioinformatics best-practice analysis pipeline for Microbial Genome Analysis.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Documentation

The gene2dis/mgap documentation is available in the [docs](docs/) folder:

- [Usage](docs/usage.md) - An overview of how the pipeline works, how to run it, and a description of all command-line flags
- [Output](docs/output.md) - An overview of the different results produced by the pipeline and how to interpret them
- [Workflow Diagram](docs/workflow_diagram.md) - A detailed Mermaid diagram of the full pipeline architecture

## Pipeline Overview

This pipeline supports three input modes:

- **Illumina short-read assembly** (FastP → SPAdes)
- **Oxford Nanopore long-read assembly** (fastplong → Flye → Medaka → Dnaapler, or Autocycler consensus → Dnaapler)
- **Pre-assembled contigs** (direct annotation)

### Genome Assembly

The genome assembly workflow processes Illumina and Oxford Nanopore data using technology-specific tools:

#### Illumina Assembly

1. Read QC, cleaning, and filtering ([`FastP`](https://github.com/OpenGene/fastp))
2. Optional coverage estimation ([`Mash`](https://mash.readthedocs.io/en/latest/)) and reduction ([`Seqtk`](https://github.com/lh3/seqtk))
3. Optional contamination check ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
4. Genome assembly ([`SPAdes`](https://github.com/ablab/spades))

#### ONT Assembly (Flye mode, default)

1. Read QC, quality filtering, and adapter trimming ([`fastplong`](https://github.com/OpenGene/fastplong))
2. Genome assembly from all reads ([`Flye`](https://github.com/fenderglass/Flye)); optionally from the longest reads up to `--flye_asm_coverage`, with the genome size from `--genome_size` or [`lrge`](https://github.com/mbhall88/lrge)
3. Genome polishing ([`Medaka`](https://github.com/nanoporetech/medaka))
4. Optional contig reorientation ([`Dnaapler`](https://github.com/gbouras13/dnaapler))

#### ONT Assembly (Autocycler mode)

1. Read QC, quality filtering, and adapter trimming ([`fastplong`](https://github.com/OpenGene/fastplong))
2. Consensus multi-assembler assembly ([`Autocycler`](https://github.com/rrwick/Autocycler))
3. Optional contig reorientation on GFA ([`Dnaapler`](https://github.com/gbouras13/dnaapler)) + GFA-to-FASTA conversion ([`Autocycler gfa2fasta`](https://github.com/rrwick/Autocycler/wiki/Autocycler-gfa2fasta))

#### Pre-assembled Contigs

- Direct annotation workflow (skips assembly steps)

### Genome Annotation

With the assembled genome (or provided contigs), contigs shorter than `--min_contig_length` (default 1000 bp) are removed
and the rest written uncompressed and renamed `contig_1`, `contig_2`, ... to one FASTA per sample
([`SeqKit`](https://github.com/shenwei356/seqkit)). Every annotation step reads that file; QUAST reports the whole
assembly. The annotation steps include:

1. Quality assessment ([`CheckM2`](https://github.com/chklovski/CheckM2))
2. MLST analysis ([`MLST`](https://github.com/tseemann/mlst))
3. Annotation ([`Bakta`](https://github.com/oschwengers/bakta))
4. Taxonomic classification ([`GTDB-Tk`](https://github.com/Ecogenomics/GTDBTk)) - optional, uses GTDB reference database
5. Antibiotic resistance prediction ([`AMRFinderPlus`](https://github.com/ncbi/amr)) - includes point mutation analysis for supported organisms
6. Antibiotic resistance gene prediction ([`RGI`](https://github.com/arpcard/rgi)) - optional, uses CARD database for resistome analysis
7. Mobile element detection ([`geNomad`](https://github.com/apcamargo/genomad)) - prophages and plasmids
8. Plasmid detection, typing, and reconstruction ([`MOB-suite`](https://github.com/phac-nml/mob-suite)) - optional
9. Species-specific analyses:
   - _Klebsiella_: [`Kleborate`](https://github.com/klebgenomics/Kleborate)
   - _S. aureus_: SCCmec classification using [`sccmec`](https://github.com/rpetit3/sccmec)
   - _Salmonella_: Serotype prediction using [`SISTR`](https://github.com/phac-nml/sistr_cmd)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=26.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility. Conda is also supported as a last resort.

3. Download the pipeline by cloning the repository or downloading the zip file

## Running the Pipeline

For detailed usage instructions, see the [Usage documentation](docs/usage.md).

### Samplesheet

Create a CSV samplesheet with your samples. Examples:

**Illumina reads:**

```csv
sample,fastq_1,fastq_2
SAMPLE1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
```

**ONT reads:**

```csv
sample,fastq_1
SAMPLE1,/path/to/sample1.fastq.gz
```

**Pre-assembled contigs:**

```csv
sample,fasta
SAMPLE1,/path/to/sample1.fasta
SAMPLE2,/path/to/sample2.fa.gz
```

An [example samplesheet](assets/samplesheet.csv) is provided with the pipeline.

**Automatic samplesheet generation:**

A helper script is available to automatically generate samplesheets from a directory of sequencing files:

```bash
python accessory_scripts/CreateSampleSheet.py /path/to/data samplesheet.csv
```

The script auto-detects the data type (Illumina/ONT/contigs) and intelligently extracts sample names from filenames. It automatically generates the correct column format based on the detected data type:

- Illumina: `sample,fastq_1,fastq_2`
- ONT: `sample,fastq_1`
- Contigs: `sample,fasta`

See the [Usage documentation](docs/usage.md#automatic-samplesheet-generation) for detailed examples.

### Basic Usage

```bash
# Illumina short-read data
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type illumina \
    -profile docker

# Oxford Nanopore long-read data
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type ont \
    -profile docker

# Pre-assembled contigs
nextflow run gene2dis/mgap \
    --input samplesheet.csv \
    --outdir results \
    --seq_type contig \
    -profile docker
```

### Parameters

Every user-facing parameter, grouped as in `nextflow_schema.json`. Database parameters default to `null`; the step
that needs them is skipped when they are unset (CheckM2, Bakta, AMRFinderPlus and geNomad print a warning at launch).

#### Input and output

| Parameter         | Default  | Description                                             |
| ----------------- | -------- | ------------------------------------------------------- |
| `--input`         | required | Path to samplesheet CSV file                            |
| `--outdir`        | required | Output directory for results                            |
| `--seq_type`      | required | Sequencing type: `illumina`, `ont`, or `contig`         |
| `--multiqc_title` | `null`   | MultiQC report title (page header and file name)        |
| `-profile`        | —        | Configuration profile: `docker`, `singularity`, `conda` |

#### Pipeline steps

| Parameter           | Default | Description                                                  |
| ------------------- | ------- | ------------------------------------------------------------ |
| `--run_kraken2`     | `true`  | Enable Kraken2 contamination detection (needs `--kraken2db`) |
| `--adjust_coverage` | `true`  | Enable coverage estimation and read subsampling              |
| `--max_coverage`    | `110`   | Target maximum coverage (x) for subsampling                  |
| `--run_gtdbtk`      | `false` | Enable GTDB-Tk taxonomic classification                      |
| `--run_rgi`         | `false` | Enable RGI antimicrobial resistance gene prediction          |
| `--run_mobsuite`    | `false` | Enable MOB-suite plasmid detection and reconstruction        |
| `--run_dnaapler`    | `true`  | Enable contig reorientation with Dnaapler (ONT only)         |

#### Databases

| Parameter                  | Default | Description                                                                                   |
| -------------------------- | ------- | --------------------------------------------------------------------------------------------- |
| `--bakta_db`               | `null`  | Path to Bakta database (Bakta and AMRFinderPlus are skipped if unset)                         |
| `--checkm2_db`             | `null`  | Path to CheckM2 database file (CheckM2 is skipped if unset)                                   |
| `--amrfinder_db`           | `null`  | Path to AMRFinderPlus database (AMRFinderPlus is skipped if unset; also needs `--bakta_db`)   |
| `--genomad_db`             | `null`  | Path to geNomad database (geNomad is skipped if unset)                                        |
| `--kraken2db`              | `null`  | Path to Kraken2 database (required to run contamination detection)                            |
| `--brackendb`              | `null`  | Path to Bracken database (Illumina only; falls back to `--kraken2db`)                         |
| `--gtdbtk_db`              | `null`  | Path to GTDB-Tk database (required if `--run_gtdbtk` is enabled)                              |
| `--gtdbtk_pplacer_scratch` | `true`  | Use a scratch directory for pplacer to reduce GTDB-Tk memory use                              |
| `--rgi_db`                 | `null`  | Path to RGI/CARD database (required if `--run_rgi` is enabled)                                |
| `--mobsuite_db`            | `null`  | Path to pre-built MOB-suite database (optional; uses bundled DB if unset)                     |
| `--mlst_blastdb`           | `null`  | Alternative MLST BLAST database directory containing `mlst.fa` and its index files (optional) |
| `--mlst_datadir`           | `null`  | Alternative MLST PubMLST data directory (optional; uses the one bundled with mlst if unset)   |
| `--plassembler_db`         | `null`  | Path to Plassembler database (only needed when `plassembler` is in `--autocycler_assemblers`) |

#### Illumina read processing (fastp)

| Parameter                           | Default | Description                                      |
| ----------------------------------- | ------- | ------------------------------------------------ |
| `--fastp_qualified_quality`         | `15`    | Minimum base quality for a base to be qualified  |
| `--fastp_unqualified_percent_limit` | `40`    | Maximum percentage of unqualified bases per read |
| `--fastp_cut_front_window_size`     | `4`     | Window size for front trimming                   |
| `--fastp_cut_front_mean_quality`    | `20`    | Mean quality threshold for front trimming        |
| `--fastp_cut_right_window_size`     | `4`     | Window size for right trimming                   |
| `--fastp_cut_right_mean_quality`    | `20`    | Mean quality threshold for right trimming        |
| `--fastp_reads_minlength`           | `50`    | Minimum read length after trimming               |
| `--fastp_n_base_limit`              | `5`     | Maximum number of N bases per read               |

#### ONT read processing and assembly

| Parameter                               | Default                                        | Description                                                                                                                           |
| --------------------------------------- | ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `--ont_assembler`                       | `flye`                                         | `flye` (Flye + Medaka) or `autocycler` (consensus of several assemblers)                                                              |
| `--fastplong_qualified_quality`         | `15`                                           | Minimum base quality for a base to be qualified (fastplong)                                                                           |
| `--fastplong_unqualified_percent_limit` | `40`                                           | Maximum percentage of unqualified bases per read (fastplong)                                                                          |
| `--fastplong_min_read_length`           | `1000`                                         | Minimum read length (fastplong)                                                                                                       |
| `--flye_mode`                           | `nano-hq`                                      | Flye read mode: `nano-hq`, `nano-raw`, `nano-corr`, `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`                                        |
| `--flye_asm_coverage`                   | `null`                                         | Flye `--asm-coverage` for very deep samples (e.g. `40`): first assembly step from the longest reads up to this coverage; off if unset |
| `--genome_size`                         | `null`                                         | Genome size for `--flye_asm_coverage`, e.g. `5.5m`; estimated per sample with lrge if unset                                           |
| `--medaka_model`                        | `null`                                         | Medaka model matching the basecaller model (`medaka tools list_models`); medaka's default if unset                                    |
| `--autocycler_assemblers`               | `raven,miniasm,flye,metamdbg,necat,nextdenovo` | Assemblers for Autocycler; also supported: `canu`, `myloasm`, `plassembler`                                                           |
| `--autocycler_read_type`                | `ont_r10`                                      | Read type for the Autocycler assembler helper: `ont_r9`, `ont_r10`, `pacbio_clr`, `pacbio_hifi`                                       |
| `--autocycler_max_contigs`              | `50`                                           | Maximum contigs per assembly in the Autocycler compress and cluster steps                                                             |

#### RGI

| Parameter              | Default   | Description                                                           |
| ---------------------- | --------- | --------------------------------------------------------------------- |
| `--rgi_include_loose`  | `false`   | Include loose hits (below the detection model cut-off)                |
| `--rgi_include_nudge`  | `false`   | Include hits nudged from loose to strict for partial gene sequences   |
| `--rgi_low_quality`    | `false`   | Low quality/coverage mode for short contigs, to predict partial genes |
| `--rgi_alignment_tool` | `DIAMOND` | Alignment tool for RGI searches: `DIAMOND` or `BLAST`                 |

#### Annotation

| Parameter             | Default | Description                                                                                        |
| --------------------- | ------- | -------------------------------------------------------------------------------------------------- |
| `--min_contig_length` | `1000`  | Minimum contig length (bp) for every annotation and typing tool (QUAST reports the whole assembly) |

> **Note:** Pipeline parameters use double dashes (`--`), while Nextflow parameters use a single dash (`-`).

For a complete list of parameters and advanced configuration options, see the [Usage documentation](docs/usage.md).

## Output

For detailed information about the pipeline outputs, see the [Output documentation](docs/output.md).

## Credits

gene2dis/mgap was originally written by the Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello. Its development was led by Juan A. Ugalde.

## Development

`main` is the production pipeline; development happens on `dev`. Work on short-lived feature branches off `dev` and
open pull requests into `dev`; `dev` is merged into `main` only for a release. Before merging, run `pre-commit run --all-files` and the test profile (`nextflow run . -profile test,docker --outdir results`).

Test-profile runs clean up their `work/` task directories automatically. To remove all remaining run/test scratch (all gitignored): `rm -rf work .nextflow .nextflow.log* .nf-test .nf-test.log`.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
