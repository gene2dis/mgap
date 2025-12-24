## Introduction

**Microbial Genome Analysis Pipeline**

**MGAP** is a bioinformatics best-practice analysis pipeline for Microbial Genome Analysis.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Documentation

The gene2dis/mgap documentation is available in the [docs](docs/) folder:

- [Usage](docs/usage.md) - An overview of how the pipeline works, how to run it, and a description of all command-line flags
- [Output](docs/output.md) - An overview of the different results produced by the pipeline and how to interpret them

## Pipeline Overview

This pipeline supports three input modes:

- **Illumina short-read assembly** (FastP → SPAdes)
- **Oxford Nanopore long-read assembly** (Nanoq → Porechop → Flye → Medaka)
- **Pre-assembled contigs** (direct annotation) 

### Genome Assembly

The genome assembly workflow processes Illumina and Oxford Nanopore data using technology-specific tools:

#### Illumina Assembly
1. Read QC, cleaning, and filtering ([`FastP`](https://github.com/OpenGene/fastp))
2. Optional coverage estimation ([`Mash`](https://mash.readthedocs.io/en/latest/)) and reduction ([`Seqtk`](https://github.com/lh3/seqtk))
3. Optional contamination check ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
4. Genome assembly ([`SPAdes`](https://github.com/ablab/spades))

#### ONT Assembly
1. Read QC ([`Nanoq`](https://github.com/esteinig/nanoq))
2. Adapter removal ([`Porechop_ABI`](https://github.com/bonsai-team/Porechop_ABI))
3. Optional coverage estimation ([`Mash`](https://mash.readthedocs.io/en/latest/)) and reduction ([`Seqtk`](https://github.com/lh3/seqtk))
4. Genome assembly ([`Flye`](https://github.com/fenderglass/Flye))
5. Genome polishing ([`Medaka`](https://github.com/nanoporetech/medaka))

#### Pre-assembled Contigs
- Direct annotation workflow (skips assembly steps)

### Genome Annotation

With the assembled genome (or provided contigs), the annotation steps include:

1. Quality assessment ([`CheckM2`](https://github.com/chklovski/CheckM2))
2. MLST analysis ([`MLST`](https://github.com/tseemann/mlst))
3. Annotation ([`Bakta`](https://github.com/oschwengers/bakta))
4. Antibiotic resistance prediction ([`AMRFinderPlus`](https://github.com/ncbi/amr)) - includes point mutation analysis for supported organisms
5. Mobile element detection ([`geNomad`](https://github.com/apcamargo/genomad)) - prophages and plasmids
6. Species-specific analyses:
   - _Klebsiella_: [`Kleborate`](https://github.com/klebgenomics/Kleborate)
   - _S. aureus_: SCCmec classification using [`staphopia-sccmec`](https://github.com/staphopia/staphopia-sccmec)


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

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
```

An [example samplesheet](assets/samplesheet.csv) is provided with the pipeline.

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

### Key Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV file |
| `--outdir` | Output directory for results |
| `--seq_type` | Sequencing type: `illumina`, `ont`, or `contig` |
| `-profile` | Configuration profile: `docker`, `singularity`, `conda` |
| `--bakta_db` | Path to Bakta database |
| `--checkm2_db` | Path to CheckM2 database |

> **Note:** Pipeline parameters use double dashes (`--`), while Nextflow parameters use a single dash (`-`).

For a complete list of parameters and advanced configuration options, see the [Usage documentation](docs/usage.md).


## Output

For detailed information about the pipeline outputs, see the [Output documentation](docs/output.md).

## Credits

gene2dis/mgap was originally written by the Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello. Its development was led by Juan A. Ugalde.


## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  gene2dis/mgap for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
