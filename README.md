## Introduction

**Microbial Genome Analysis Pipeline**

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**MGAP** is a bioinformatics best-practice analysis pipeline for Microbial Genome Analysis Pipeline.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC, clean, and filter reads. ([`FastP`](https://github.com/OpenGene/fastp))
2. If requested, calculate coverage of genome ([`Mash`](https://mash.readthedocs.io/en/latest/)) and reduce coverage ([`Seqtk`](https://github.com/lh3/seqtk).
3. If requested, check the reads with Kraken2 against a standard (small) database to evaluate possible contaminations ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
4. Genome assembly ([`Unicycler`](https://github.com/rrwick/Unicycler))
5. MLST analysis ([`MLST`](https://github.com/tseemann/mlst))
6. Annotation ([`BAKTA`](https://github.com/oschwengers/bakta))
7. Antibiotic resistance prediction using AMRFinderPlus. If the organism is on the list provided by AMRFinderPlus, also point mutations are evaluated ([`AMRFinderPLus`](https://github.com/ncbi/amr))
8. Prophage and plasmid search using [`Genomad`](https://github.com/apcamargo/genomad)
9. For _Klebsiella_, evaluation of the genome using [`Kleborate`](https://github.com/klebgenomics/Kleborate)
10. For _S. aureus_, SCCmec classification using [`staphophia-sccmec`](https://github.com/staphopia/staphopia-sccmec)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

Currently the pipeline has been tested with Docker, but should work without issues with Singularity or Conda.

3. Download the pipeline, either cloning the repository or downloading the zip file

## Running the pipeline

### Sample sheet file
First you need to create a Samplesheet file, that contain the name of the samples and the location of the reads. This is an example of this file:

```
sample,fastq_1,fastq_2
SCL10095,/home/jugalde/germlab/data/S190_R1.fastq.gz,/home/jugalde/germlab/data/S190_R2.fastq.gz
```

The file **always** has to include the header `sample,fastq_1,fastq_2`

### Pipeline parameters

The file `nextflow.config` contains all the parameteres used by the pipeline, including path to database files. Currently the path work in our server (_Arrakis_), but if you are running elsewhere, these need to be updated. *Later a section on how to download and prepare the DBs needed for the analysis, will be added to this Readme file*.


#### Starting the pipeline

Once you have everything ready, you can run the pipeline by calling the name of the repository, and adding the appropiate parameters. For example:

```bash
nextflow run 
```

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run gene2dis-mgap --input <SAMPLESHEET> --outdir <OUTDIR> -profile docker 
   ```

   Here the parameters are:
   - gene2dis-mgap: The name of the repository. This could be the path to a location on the computer (e.g. /home/jugalde/pipelines/gene2dis-mgap).
   - --input: The samplesheet file
   - --outdir: The output directory for the output of the pipeline
   - -profile: either docker or conda (not tested yet)

   *Notice that the pipeline parameters have two dashes (--), while parameters that are for nextflow only have one (-).


## Credits

gene2dis/mgap was originally written by the Microbial Data Science Lab, Center for Bioinformatics and Integrative Biology, Universidad Andres Bello. Its development was led by Juan A. Ugalde

We thank the following people for their extensive assistance in the development of this pipeline:

- Juan A. Ugalde


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
