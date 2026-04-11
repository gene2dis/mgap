# gene2dis/mgap: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Output Folder Structure

All results are written to the directory specified with `--outdir`. The layout below shows the full directory tree; conditional directories are annotated with the flag that controls them.

```
<outdir>/
├── <sample_id>/                         # One subdirectory per sample
│   ├── read_processing/
│   │   ├── fastp/                       # Illumina QC: *.html, *.json, *.fastq.gz
│   │   ├── fastplong/                   # ONT QC: *.html, *.json, *.fastq.gz
│   │   ├── kraken2/                     # Contamination detection: *.report.txt, *.fastq.gz
│   │   │                                #   (requires --run_kraken2 true and --kraken2db)
│   │   ├── bracken/                     # Abundance re-estimation: *.tsv
│   │   │                                #   (Illumina only; requires --brackendb)
│   │   └── subsampled/                  # Subsampled reads: *.fastq.gz
│   │                                    #   (only when coverage > --max_coverage)
│   ├── assemblies/
│   │   ├── *.scaffolds.fa.gz            # SPAdes output (--seq_type illumina)
│   │   ├── flye/                        # Flye assembly: *.fasta, *.gfa, *.txt, *.log
│   │   │                                #   (ONT, --ont_assembler flye)
│   │   ├── medaka/                      # Medaka-polished assembly: *.fasta
│   │   │                                #   (ONT, --ont_assembler flye)
│   │   ├── autocycler/                  # Autocycler consensus assembly
│   │   │   ├── genome_size/             #   Genome size estimate: *_genome_size.txt
│   │   │   ├── autocycler_out/          #   Cluster/trim/resolve working directory
│   │   │   ├── *.fasta                  #   Final consensus assembly
│   │   │   └── *.gfa                    #   Final consensus assembly graph
│   │   │                                #   (ONT, --ont_assembler autocycler)
│   │   └── dnaapler/                    # Reoriented assembly: *_reoriented.fasta / *.gfa / *.fasta
│   │                                    #   (ONT only, --run_dnaapler true [default])
│   ├── qc/
│   │   └── quast/                       # Assembly quality metrics: *.tsv
│   └── annotation/
│       ├── checkm2/                     # Genome completeness/contamination: *.tsv
│       ├── mlst/                        # Multi-locus sequence typing: *.tsv
│       ├── bakta/                       # Genome annotation: *.gff3, *.gbff, *.fna, *.faa, *.tsv, ...
│       ├── amrfinder/                   # AMR detection: *.tsv
│       ├── genomad/                     # Mobile genetic elements: *_summary/, *_annotate/, *_find_proviruses/
│       ├── mobsuite/                    # Plasmid detection: chromosome.fasta, contig_report.txt, plasmid_*.fasta, mobtyper_results.txt
│       │                                #   (requires --run_mobsuite true)
│       ├── rgi/                         # Resistance gene prediction: *.txt, *.json
│       │                                #   (requires --run_rgi true)
│       ├── kleborate/                   # Klebsiella virulence/resistance: *.txt
│       │                                #   (auto-triggered when Klebsiella is detected by MLST)
│       ├── sccmec/                      # S. aureus SCCmec typing: *.tsv
│       │                                #   (auto-triggered when S. aureus is detected by MLST)
│       └── sistr/                       # Salmonella serotype prediction: *.tab, *-allele.{json,fasta}, *-cgmlst.csv
│                                        #   (auto-triggered when Salmonella is detected by MLST)
├── gtdbtk/                              # GTDB-Tk taxonomic classification (batch, all samples)
│   ├── gtdbtk.batch.*.summary.tsv       #   (requires --run_gtdbtk true and --gtdbtk_db)
│   └── gtdbtk.batch.log
└── pipeline_info/                       # Nextflow execution reports and software versions
    ├── execution_report_*.html
    ├── execution_timeline_*.html
    ├── execution_trace_*.txt
    ├── pipeline_dag_*.html
    └── software_versions.yml
```

> **Notes:**
> - `<sample_id>` directories are created for each sample. `gtdbtk/` is a single shared directory processed in batch mode across all samples.
> - `bracken/` only appears for `--seq_type illumina`; ONT mode runs Kraken2 only.
> - `subsampled/` only appears when reads exceed `--max_coverage` and `--adjust_coverage true` (default).
> - `dnaapler/` also receives `autocycler gfa2fasta` output when using Autocycler mode.
> - `rgi/`, `mobsuite/`, `kleborate/`, `sccmec/`, and `sistr/` are conditional and only appear when triggered.

---

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

### Read Processing

- [FastP](#fastp) - Illumina read QC and trimming
- [fastplong](#fastplong) - ONT read QC, filtering, and adapter trimming
- [Kraken2 + Bracken](#kraken2--bracken) - Contamination detection *(optional)*
- [Mash + Seqtk](#mash--seqtk) - Coverage estimation and read subsampling *(optional)*

### Genome Assembly

- [SPAdes](#spades) - Illumina *de novo* assembly
- [Flye](#flye) - ONT *de novo* assembly (default ONT assembler)
- [Medaka](#medaka) - ONT assembly polishing (Flye mode)
- [Autocycler](#autocycler) - ONT consensus multi-assembler assembly *(optional)*
- [Dnaapler](#dnaapler) - Contig reorientation *(optional, ONT only)*

### Genome Annotation and Analysis

- [QUAST](#quast) - Assembly quality metrics
- [CheckM2](#checkm2) - Genome completeness and contamination assessment
- [MLST](#mlst) - Multi-locus sequence typing
- [Bakta](#bakta) - Genome annotation
- [GTDB-Tk](#gtdb-tk) - Taxonomic classification *(optional)*
- [AMRFinderPlus](#amrfinderplus) - Antimicrobial resistance detection
- [geNomad](#genomad) - Mobile genetic element identification
- [RGI](#rgi) - Resistance gene prediction *(optional)*
- [MOB-suite](#mob-suite) - Plasmid detection and reconstruction *(optional)*

### Species-Specific Analysis

- [Kleborate](#kleborate) - *Klebsiella*-specific virulence and resistance scoring
- [sccmec](#sccmec) - *S. aureus* SCCmec cassette typing
- [SISTR](#sistr) - *Salmonella* serotype prediction

### Reporting

- [Pipeline information](#pipeline-information) - Report metrics generated during workflow execution

---

## Read Processing

### FastP

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/read_processing/fastp/`
  - `*.fastp.html`: Interactive HTML report with quality metrics before and after trimming.
  - `*.fastp.json`: JSON report with detailed trimming and quality statistics.
  - `*.fastq.gz`: Trimmed and filtered reads.

</details>

[FastP](https://github.com/OpenGene/fastp) performs quality trimming, adapter removal, and filtering of Illumina paired-end reads. The HTML report provides a visual summary of read quality before and after processing.

This step runs for `--seq_type illumina` only.

### fastplong

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/read_processing/fastplong/`
  - `*.html`: Interactive HTML report with quality metrics.
  - `*.json`: JSON report with detailed filtering statistics.
  - `*.fastq.gz`: Filtered and trimmed long reads.

</details>

[fastplong](https://github.com/OpenGene/fastplong) performs quality filtering and adapter trimming of Oxford Nanopore long reads. It replaces the need for separate tools for quality filtering and adapter removal.

This step runs for `--seq_type ont` only.

### Kraken2 + Bracken

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/read_processing/kraken2/`
  - `*.report.txt`: Kraken2 classification report with taxonomic abundance.
  - `*.classifiedreads.fastq.gz`: Reads classified by Kraken2.
  - `*.unclassifiedreads.fastq.gz`: Reads not classified by Kraken2.
- `<sample_id>/read_processing/bracken/`
  - `*.tsv`: Bracken abundance re-estimation at a specified taxonomic level.

</details>

[Kraken2](https://ccb.jhu.edu/software/kraken2/) classifies reads against a reference database for contamination detection. [Bracken](https://ccb.jhu.edu/software/bracken/) refines Kraken2 abundance estimates at a specified taxonomic level.

This step is optional and runs when `--run_kraken2` is enabled and `--kraken2db` is provided.

### Mash + Seqtk

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/read_processing/mash/`
  - `*.msh`: Mash sketch file.
  - `*.mash_stats`: Mash sketch statistics.
  - `*.mash_coverage`: Estimated genome coverage.
- `<sample_id>/read_processing/subsampled/`
  - `*.fastq.gz`: Subsampled reads (only present when coverage exceeded `--max_coverage`).

</details>

[Mash](https://mash.readthedocs.io/en/latest/) estimates genome coverage from read sketches. If coverage exceeds the `--max_coverage` threshold (default: 110x), [Seqtk](https://github.com/lh3/seqtk) subsamples reads to the target coverage.

This step is optional and runs when `--adjust_coverage` is enabled (default: true). It is skipped when using `--ont_assembler autocycler`, as Autocycler handles its own read subsampling.

---

## Genome Assembly

### SPAdes

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/assemblies/`
  - `*.scaffolds.fa.gz`: Gzipped scaffold sequences (primary assembly output).
  - `*.contigs.fa.gz`: Gzipped contig sequences.
  - `*.assembly.gfa.gz`: Gzipped assembly graph in GFA format.
  - `*.spades.log`: SPAdes log file.

</details>

[SPAdes](https://github.com/ablab/spades) performs *de novo* genome assembly from Illumina short reads using the `--isolate` mode optimized for bacterial isolate genomes.

This step runs for `--seq_type illumina` only.

### Flye

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/assemblies/flye/`
  - `*.fasta`: Assembled genome contigs.
  - `*.gfa`: Assembly graph in GFA format.
  - `*.log`: Flye log file.
  - `*.txt`: Assembly info with contig statistics.

</details>

[Flye](https://github.com/fenderglass/Flye) performs *de novo* assembly of long reads. The assembly mode is controlled by `--flye_mode` (default: `nano-hq`).

This step runs for `--seq_type ont` with `--ont_assembler flye` (default).

### Medaka

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/assemblies/medaka/`
  - `*.fasta`: Polished genome assembly.

</details>

[Medaka](https://github.com/nanoporetech/medaka) polishes the Flye assembly using the original long reads to improve consensus accuracy.

This step runs for `--seq_type ont` with `--ont_assembler flye` (default).

### Autocycler

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/assemblies/autocycler/`
  - `genome_size/*_genome_size.txt`: Estimated genome size.
  - `autocycler_out/`: Autocycler working directory with clustering, trimming, and resolution results (excludes intermediate `input_assemblies/` and `subsampled_reads/` subdirectories).
  - `*.fasta`: Final consensus assembly in FASTA format.
  - `*.gfa`: Final consensus assembly in GFA format.

</details>

[Autocycler](https://github.com/rrwick/Autocycler) generates consensus long-read assemblies by running multiple assemblers on subsampled reads and combining the results. The pipeline runs the full Autocycler workflow: genome size estimation, subsampling, multi-assembler assembly, compression, clustering, trim/resolve, and combine.

This step runs for `--seq_type ont` with `--ont_assembler autocycler`.

### Dnaapler

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/assemblies/dnaapler/`
  - `*_reoriented.fasta`: Reoriented genome assembly in FASTA format (Flye mode).
  - `*_reoriented.gfa`: Reoriented genome assembly in GFA format (Autocycler mode).
  - `*.fasta`: Final FASTA assembly converted from reoriented GFA (Autocycler mode, via `autocycler gfa2fasta`).

</details>

[Dnaapler](https://github.com/gbouras13/dnaapler) reorients complete circular microbial genome assemblies so that each sequence starts at a consistent location, typically at a gene like *dnaA* (chromosomes), *repA* (plasmids), or *terL* (phages).

- **Flye mode:** Operates on FASTA input from Medaka, reorienting all contigs.
- **Autocycler mode:** Operates on GFA input from `autocycler combine`, reorienting only circular contigs. The reoriented GFA is then converted to FASTA via `autocycler gfa2fasta`.

This step is optional and enabled by default (`--run_dnaapler true`). It only applies to ONT assemblies.

---

## Genome Annotation and Analysis

### QUAST

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/qc/quast/`
  - `*.tsv`: Tab-separated assembly quality metrics.

</details>

[QUAST](https://github.com/ablab/quast) evaluates genome assembly quality by computing metrics such as total length, number of contigs, N50, GC content, and largest contig size.

### CheckM2

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/checkm2/`
  - `*.tsv`: Tab-separated completeness and contamination assessment.

</details>

[CheckM2](https://github.com/chklovski/CheckM2) assesses genome quality by estimating completeness and contamination using machine learning models. The output reports completeness percentage, contamination percentage, and overall quality assessment for each genome.

### MLST

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/mlst/`
  - `*.tsv`: Tab-separated MLST typing results including scheme, sequence type, and allele profiles.

</details>

[MLST](https://github.com/tseemann/mlst) scans genome assemblies against PubMLST schemes to determine sequence types. The MLST result is also used internally to identify species for AMRFinderPlus organism-specific analysis and to trigger species-specific tools (Kleborate, sccmec, SISTR).

### Bakta

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/bakta/`
  - `*.gff3`: Annotation in GFF3 format.
  - `*.gbff`: Annotation in GenBank format.
  - `*.fna`: Annotated nucleotide sequences.
  - `*.faa`: Annotated protein sequences.
  - `*.tsv`: Tab-separated annotation summary.
  - `*.txt`: Human-readable annotation summary.
  - `*.embl`: Annotation in EMBL format.
  - `*.hypotheticals.tsv`: Hypothetical protein details.
  - `*.hypotheticals.faa`: Hypothetical protein sequences.

</details>

[Bakta](https://github.com/oschwengers/bakta) provides rapid and standardized annotation of bacterial genomes. It identifies CDSs, tRNAs, rRNAs, ncRNAs, CRISPR arrays, and other genomic features. The Bakta outputs (nucleotide FASTA, protein FASTA, GFF3) are used as input for downstream AMRFinderPlus and geNomad analyses.

### GTDB-Tk

<details markdown="1">
<summary>Output files</summary>

- `gtdbtk/`
  - `gtdbtk.batch.*.summary.tsv`: Taxonomic classification summary for all bacterial and archaeal genomes in the batch.
  - `gtdbtk.batch.*.classify.tree.gz`: Phylogenetic tree with genome placement *(optional)*.
  - `gtdbtk.batch.*.markers_summary.tsv`: Summary of identified marker genes *(optional)*.
  - `gtdbtk.batch.*.msa.fasta.gz`: Multiple sequence alignment of marker genes *(optional)*.
  - `gtdbtk.batch.*.filtered.tsv`: Genomes filtered during classification *(optional)*.
  - `gtdbtk.batch.failed_genomes.tsv`: Genomes that failed classification *(optional)*.
  - `gtdbtk.batch.log`: GTDB-Tk execution log.
  - `gtdbtk.batch.warnings.log`: Warnings generated during classification.

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

### AMRFinderPlus

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/amrfinder/`
  - `*.tsv`: Tab-separated AMR detection results.

</details>

[AMRFinderPlus](https://github.com/ncbi/amr) identifies antimicrobial resistance genes, stress response genes, and virulence factors in assembled genomes. It uses the Bakta-annotated nucleotide sequences, protein sequences, and GFF3 annotations alongside species information derived from MLST to perform organism-specific point mutation analysis for supported species.

The output TSV contains:
- Gene symbol and name
- Protein/nucleotide accession
- Sequence identity and coverage
- AMR gene family and subclass
- Element type (AMR, STRESS, VIRULENCE)
- Method of detection (BLAST, HMM, or point mutation)

### geNomad

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/genomad/`
  - `*_summary/`: Summary directory with classification results.
  - `*_annotate/`: Annotation details for identified mobile elements.
  - `*_find_proviruses/`: Provirus detection results.

</details>

[geNomad](https://github.com/apcamargo/genomad) identifies mobile genetic elements (prophages and plasmids) in genome assemblies. It uses the Bakta-annotated nucleotide sequences as input and classifies contigs as chromosome, plasmid, or virus based on gene content and genomic signatures.

### RGI

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/rgi/`
  - `<sample_id>.txt`: Tab-separated file containing RGI predictions with detailed information about detected resistance genes.
  - `<sample_id>.json`: JSON file containing comprehensive RGI results including gene sequences and additional metadata.

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

### MOB-suite

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/mobsuite/`
  - `chromosome.fasta`: All contigs assigned to the chromosome.
  - `contig_report.txt`: Per-contig chromosome/plasmid assignment table.
  - `plasmid_*.fasta`: Reconstructed plasmid sequences, one FASTA per plasmid group.
  - `mobtyper_results.txt`: Plasmid typing results (replicon, MOB-typer predictions).

</details>

[MOB-suite](https://github.com/phac-nml/mob-suite) reconstructs and types plasmids in draft bacterial assemblies using the `mob_recon` workflow. It clusters contigs into plasmid groups, assigns replicon and relaxase types, and produces one FASTA per reconstructed plasmid together with a contig-level chromosome/plasmid assignment report.

This step is optional and only runs when `--run_mobsuite` is enabled. A pre-built database can be supplied via `--mobsuite_db`; otherwise MOB-suite uses the database bundled with the container.

---

## Species-Specific Analysis

### Kleborate

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/kleborate/`
  - `*.txt`: Kleborate analysis results.

</details>

[Kleborate](https://github.com/klebgenomics/Kleborate) is a tool for screening genome assemblies of _Klebsiella pneumoniae_ and the _Klebsiella pneumoniae_ species complex (KpSC) for:

- MLST sequence type
- Species identification (within the KpSC)
- ICEKp-associated virulence loci: yersiniabactin (*ybt*), colibactin (*clb*), salmochelin (*iro*), hypermucoidy (*rmpA*)
- Virulence plasmid associated loci: salmochelin (*iro*), aerobactin (*iuc*), hypermucoidy (*rmpA*, *rmpA2*)
- Antimicrobial resistance determinants: acquired genes, SNPs, gene truncations and intrinsic β-lactamases
- K (capsule) and O antigen (LPS) serotype prediction

**Automatic detection:** Kleborate analysis is automatically triggered when _Klebsiella_ species are detected based on MLST results.

### sccmec

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/sccmec/`
  - `*.tsv`: SCCmec cassette typing results.

</details>

[sccmec](https://github.com/rpetit3/sccmec) classifies the SCCmec cassette type in *Staphylococcus aureus* genome assemblies. SCCmec (Staphylococcal Cassette Chromosome *mec*) is a mobile genetic element that carries the *mecA* gene responsible for methicillin resistance (MRSA).

**Automatic detection:** This analysis is automatically triggered when *S. aureus* is detected based on MLST results.

### SISTR

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/sistr/`
  - `<sample_id>.tab`: Tab-separated serovar prediction including serogroup, H1/H2 antigens, and cgMLST-based subtype.
  - `<sample_id>-allele.json`: Novel allele annotations in JSON format.
  - `<sample_id>-allele.fasta`: Novel allele sequences.
  - `<sample_id>-cgmlst.csv`: cgMLST profile.

</details>

[SISTR](https://github.com/phac-nml/sistr_cmd) (Salmonella In Silico Typing Resource) predicts serovars and subtypes from _Salmonella enterica_ assemblies using antigen gene detection and cgMLST. The pipeline pins SISTR to **v1.1.3**.

**Automatic detection:** SISTR analysis is automatically triggered when MLST assigns the `senterica_achtman_2` scheme to a sample.

---

## Reporting

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.html`.
  - `software_versions.yml`: Software versions used in the pipeline run.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
