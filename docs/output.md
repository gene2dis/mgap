# gene2dis/mgap: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [GTDB-Tk](#gtdb-tk) - Taxonomic classification (optional)
- [RGI](#rgi) - Antimicrobial resistance gene prediction (optional)
- [Kleborate](#kleborate) - Klebsiella-specific analysis (when detected)
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### GTDB-Tk

<details markdown="1">
<summary>Output files</summary>

- `gtdbtk/`
  - `gtdbtk.batch.*.summary.tsv`: Taxonomic classification summary for all bacterial and archaeal genomes in the batch
  - `gtdbtk.batch.*.classify.tree.gz`: Phylogenetic tree with genome placement (optional)
  - `gtdbtk.batch.*.markers_summary.tsv`: Summary of identified marker genes (optional)
  - `gtdbtk.batch.*.msa.fasta.gz`: Multiple sequence alignment of marker genes (optional)
  - `gtdbtk.batch.*.filtered.tsv`: Genomes filtered during classification (optional)
  - `gtdbtk.batch.failed_genomes.tsv`: Genomes that failed classification (optional)
  - `gtdbtk.batch.log`: GTDB-Tk execution log
  - `gtdbtk.batch.warnings.log`: Warnings generated during classification

</details>

[GTDB-Tk](https://ecogenomics.github.io/GTDBTk/) is a toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy (GTDB). It uses average nucleotide identity (ANI) to assign genomes to species clusters and phylogenetic placement for higher-level taxonomic assignments.

**Batch Processing:** GTDB-Tk processes all genomes in a single batch run, which significantly reduces runtime compared to processing genomes individually. All samples are analyzed together and results are consolidated into a single set of output files.

The main output file is the `summary.tsv` which contains:
- Taxonomic classification (domain to species level) for all genomes
- Classification method used
- Closest reference genome
- ANI to reference genome
- Alignment fraction
- Warnings and notes

This step is optional and only runs when `--run_gtdbtk` is enabled. The GTDB-Tk database must be provided via `--gtdbtk_db`.

### RGI

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/rgi/`
  - `<sample_id>.txt`: Tab-separated file containing RGI predictions with detailed information about detected resistance genes
  - `<sample_id>.json`: JSON file containing comprehensive RGI results including gene sequences and additional metadata

</details>

[RGI](https://github.com/arpcard/rgi) (Resistance Gene Identifier) predicts resistomes from protein or nucleotide sequences using the Comprehensive Antibiotic Resistance Database (CARD). The tool identifies antimicrobial resistance genes based on homology and SNP models.

The main output file (`*.txt`) contains:
- **ORF_ID**: Open reading frame identifier
- **Contig**: Source contig name
- **Cut_Off**: Detection paradigm (Strict, Perfect, or Loose)
- **Best_Hit_ARO**: Best matching ARO (Antibiotic Resistance Ontology) term
- **Drug Class**: Antibiotic drug class(es) the gene confers resistance to
- **Resistance Mechanism**: Mechanism of resistance (e.g., antibiotic efflux, antibiotic inactivation)
- **AMR Gene Family**: Gene family classification
- **% Identity**: Sequence identity to reference
- **% Coverage**: Sequence coverage of reference

The JSON output provides additional details including:
- Full nucleotide and protein sequences
- CARD database annotations
- Model information and detection parameters

This step is optional and only runs when `--run_rgi` is enabled. The CARD database must be provided via `--rgi_db`.

### Dnaapler

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/assemblies/dnaapler/`
  - `<sample_id>_reoriented.fasta`: Reoriented genome assembly in FASTA format (Flye mode)
  - `<sample_id>_reoriented.gfa`: Reoriented genome assembly in GFA format (Autocycler mode)
- `<sample_id>/assemblies/autocycler/`
  - `<sample_id>.fasta`: Final FASTA assembly converted from reoriented GFA (Autocycler mode, via `autocycler gfa2fasta`)

</details>

[Dnaapler](https://github.com/gbouras13/dnaapler) (v1.3.0) reorients complete circular microbial genome assemblies so that each sequence starts at a consistent location, typically at a gene like *dnaA* (chromosomes), *repA* (plasmids), or *terL* (phages).

- **Flye mode:** Operates on FASTA input from Medaka, reorienting all contigs.
- **Autocycler mode:** Operates on GFA input from `autocycler combine`, reorienting only circular contigs. The reoriented GFA is then converted to FASTA via `autocycler gfa2fasta`.

This step is optional and enabled by default (`--run_dnaapler true`). It only applies to ONT assemblies.

### Kleborate

<details markdown="1">
<summary>Output files</summary>

- `kleborate/`
  - `*.results.txt`: Kleborate analysis results for each Klebsiella sample

</details>

[Kleborate](https://github.com/klebgenomics/Kleborate) is a tool for screening genome assemblies of _Klebsiella pneumoniae_ and the _Klebsiella pneumoniae_ species complex (KpSC) for:

- MLST sequence type
- Species identification (within the KpSC)
- ICEKp-associated virulence loci: yersiniabactin (*ybt*), colibactin (*clb*), salmochelin (*iro*), hypermucoidy (*rmpA*)
- Virulence plasmid associated loci: salmochelin (*iro*), aerobactin (*iuc*), hypermucoidy (*rmpA*, *rmpA2*)
- Antimicrobial resistance determinants: acquired genes, SNPs, gene truncations and intrinsic β-lactamases
- K (capsule) and O antigen (LPS) serotype prediction

**Automatic detection:** Kleborate analysis is automatically triggered when _Klebsiella_ species are detected based on MLST results. The analysis runs on a per-sample basis, generating detailed reports for each identified Klebsiella genome.

The output file contains comprehensive typing and screening information that helps characterize the virulence and resistance profiles of Klebsiella isolates.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
