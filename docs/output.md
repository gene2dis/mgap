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
│   │   ├── mash/                        # Coverage estimate: *.msh, *.mash_stats, *.genome_size
│   │   │                                #   (Illumina, --adjust_coverage true [default])
│   │   └── subsampled/                  # Subsampled reads: *.fastq.gz
│   │                                    #   (only when coverage > --max_coverage)
│   ├── assemblies/
│   │   ├── *.scaffolds.fa.gz            # SPAdes output (--seq_type illumina)
│   │   ├── flye/                        # Flye assembly: *.fasta.gz, *.gfa.gz, *.txt, *.log
│   │   │                                #   (ONT, --ont_assembler flye)
│   │   ├── lrge/                        # Genome size estimate for Flye: *.genome_size.txt
│   │   │                                #   (ONT, --flye_asm_coverage without --genome_size)
│   │   ├── medaka/                      # Medaka-polished assembly: *.fasta
│   │   │                                #   (ONT, --ont_assembler flye)
│   │   ├── autocycler/                  # Autocycler consensus assembly
│   │   │   ├── genome_size/             #   Genome size estimate: *_genome_size.txt
│   │   │   ├── autocycler_out/          #   Compress/cluster/trim/resolve/combine directory
│   │   │   ├── *.fasta                  #   Final consensus assembly
│   │   │   └── *.gfa                    #   Final consensus assembly graph
│   │   │                                #   (ONT, --ont_assembler autocycler)
│   │   ├── dnaapler/                    # Reoriented assembly: *_reoriented.fasta / *.gfa / *.fasta
│   │   │                                #   (ONT only, --run_dnaapler true [default])
│   │   └── annotation_input/            # FASTA read by every annotation tool: *.fasta (contigs >= --min_contig_length,
│   │                                    #   uncompressed), *.filter_stats.tsv (contigs and length before/after)
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
├── report/                              # Run report (skipped with --skip_report)
│   ├── mgap_report.html                 #   Offline HTML report of every sample and tool
│   ├── mgap_summary.tsv                 #   One row per sample: key numbers and warnings
│   └── mgap_summary.xlsx                #   The summary, then one sheet per tool
└── pipeline_info/                       # Nextflow execution reports and software versions
    ├── execution_report_*.html
    ├── execution_timeline_*.html
    ├── execution_trace_*.txt
    ├── pipeline_dag_*.html
    ├── sample_status.tsv                 # pass, skipped or failed per sample, with reason and key numbers
    └── software_versions.yml
```

> **Notes:**
>
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
- [Kraken2 + Bracken](#kraken2--bracken) - Contamination detection _(optional)_
- [Mash + Seqtk](#mash--seqtk) - Coverage estimation and read subsampling _(optional)_

### Genome Assembly

- [SPAdes](#spades) - Illumina _de novo_ assembly
- [Flye](#flye) - ONT _de novo_ assembly (default ONT assembler)
- [Medaka](#medaka) - ONT assembly polishing (Flye mode)
- [Autocycler](#autocycler) - ONT consensus multi-assembler assembly _(optional)_
- [Dnaapler](#dnaapler) - Contig reorientation _(optional, ONT only)_
- [Annotation input](#annotation-input) - Length filter and decompression before annotation

### Genome Annotation and Analysis

- [QUAST](#quast) - Assembly quality metrics
- [CheckM2](#checkm2) - Genome completeness and contamination assessment
- [MLST](#mlst) - Multi-locus sequence typing
- [Bakta](#bakta) - Genome annotation
- [GTDB-Tk](#gtdb-tk) - Taxonomic classification _(optional)_
- [AMRFinderPlus](#amrfinderplus) - Antimicrobial resistance detection
- [geNomad](#genomad) - Mobile genetic element identification
- [RGI](#rgi) - Resistance gene prediction _(optional)_
- [MOB-suite](#mob-suite) - Plasmid detection and reconstruction _(optional)_

### Species-Specific Analysis

- [Kleborate](#kleborate) - _Klebsiella_-specific virulence and resistance scoring
- [sccmec](#sccmec) - _S. aureus_ SCCmec cassette typing
- [SISTR](#sistr) - _Salmonella_ serotype prediction

### Reporting

- [Run report](#run-report) - One offline HTML page, a per-sample TSV and an XLSX workbook for the whole run
- [MultiQC](#multiqc) - QC reports of fastp/fastplong, Kraken2, QUAST and Bakta in one page
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

[Kraken2](https://ccb.jhu.edu/software/kraken2/) classifies reads against a reference database for contamination detection. [Bracken](https://ccb.jhu.edu/software/bracken/) refines Kraken2 abundance estimates at a specified taxonomic level, with the database read length (`-r`) matched to each sample's reads (or `--bracken_read_length`).

This step is optional and runs when `--run_kraken2` is enabled and `--kraken2db` is provided.

### Mash + Seqtk

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/read_processing/mash/`
  - `*.msh`: Mash sketch file.
  - `*.mash_stats`: Mash sketch statistics (estimated genome size and k-mer coverage).
  - `*.genome_size`: Mash's genome size estimate in bases, used for the coverage.
- `<sample_id>/read_processing/subsampled/`
  - `*.fastq.gz`: Subsampled reads (only present when coverage exceeded `--max_coverage`).

</details>

[Mash](https://mash.readthedocs.io/en/latest/) estimates the genome size from the reads' k-mers; the read coverage is the number of bases fastp kept divided by that size. If coverage exceeds `--max_coverage` (default: 150x), [Seqtk](https://github.com/lh3/seqtk) subsamples the reads to that coverage.

This step is optional and runs when `--adjust_coverage` is enabled (default: true), for `--seq_type illumina` only: ONT reads are not subsampled.

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

[SPAdes](https://github.com/ablab/spades) performs _de novo_ genome assembly from Illumina short reads using the `--isolate` mode optimized for bacterial isolate genomes.

This step runs for `--seq_type illumina` only.

### Flye

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/assemblies/flye/`
  - `*.assembly.fasta.gz`: Assembled genome contigs (before Medaka polishing).
  - `*.assembly_graph.gfa.gz`: Assembly graph in GFA format.
  - `*.flye.log`: Flye log file.
  - `*.assembly_info.txt`: Assembly info with contig statistics (length, coverage, circularity).
  - `*.params.json`: Flye run parameters.

</details>

[Flye](https://github.com/fenderglass/Flye) performs _de novo_ assembly of long reads. The assembly mode is controlled by `--flye_mode` (default: `nano-hq`). Flye uses all reads that pass fastplong; with `--flye_asm_coverage` its first step (disjointigs) uses only the longest reads up to that coverage.

This step runs for `--seq_type ont` with `--ont_assembler flye` (default).

- `<sample_id>/assemblies/lrge/`
  - `*.genome_size.txt`: Genome size in bases estimated by [lrge](https://github.com/mbhall88/lrge) from the fastplong reads, passed to Flye as `--genome-size`. Only with `--flye_asm_coverage` and no `--genome_size`.

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
  - `autocycler_out/`: Autocycler directory written by compress, cluster, trim, resolve and combine (one task per sample): `input_assemblies.gfa`, `clustering/` (with `clustering.tsv` and one directory per cluster), and `consensus_assembly.{fasta,gfa,yaml}`. A sample with no cluster passing QC has no `qc_pass/` clusters and no consensus files.
  - `*.fasta`: Final consensus assembly in FASTA format.
  - `*.gfa`: Final consensus assembly in GFA format.

</details>

[Autocycler](https://github.com/rrwick/Autocycler) generates consensus long-read assemblies by running multiple assemblers on subsampled reads and combining the results. The pipeline runs the full Autocycler workflow: genome size estimation, subsampling, multi-assembler assembly, then compression, clustering, trim/resolve and combine in one task per sample.

This step runs for `--seq_type ont` with `--ont_assembler autocycler`.

### Dnaapler

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/assemblies/dnaapler/`
  - `*_reoriented.fasta`: Reoriented genome assembly in FASTA format (Flye mode).
  - `*_reoriented.gfa`: Reoriented genome assembly in GFA format (Autocycler mode).
  - `*.fasta`: Final FASTA assembly converted from reoriented GFA (Autocycler mode, via `autocycler gfa2fasta`).

</details>

[Dnaapler](https://github.com/gbouras13/dnaapler) reorients complete circular microbial genome assemblies so that each sequence starts at a consistent location, typically at a gene like _dnaA_ (chromosomes), _repA_ (plasmids), or _terL_ (phages).

- **Flye mode:** Operates on FASTA input from Medaka, reorienting all contigs.
- **Autocycler mode:** Operates on GFA input from `autocycler combine`, reorienting only circular contigs. The reoriented GFA is then converted to FASTA via `autocycler gfa2fasta`.

This step is optional and enabled by default (`--run_dnaapler true`). It only applies to ONT assemblies.

### Annotation input

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/assemblies/annotation_input/`
  - `<sample_id>.fasta`: The assembly (SPAdes scaffolds, ONT assembly or the samplesheet contigs) without contigs shorter than `--min_contig_length` (default 1000 bp), uncompressed, contigs renamed `contig_1`, `contig_2`, ... in assembly order, with the original header kept as the description (e.g. `>contig_1 NODE_1_length_1682762_cov_30.909424`).
  - `<sample_id>.filter_stats.tsv`: `seqkit stats` table with one row for the assembly before the filter and one for the file above: number of contigs, total length, N50 and other length statistics.

</details>

[SeqKit](https://bioinf.shenwei.me/seqkit/) prepares one FASTA per sample that every annotation and typing tool reads: CheckM2, MLST, Bakta, GTDB-Tk, geNomad, MOB-suite, RGI, Kleborate, sccmec and SISTR (AMRFinderPlus reads Bakta's outputs, which keep the same contig names). All of them therefore see the same contigs under the same names, and the original assembler name of any `contig_N` can be looked up in this file. QUAST runs on the unfiltered assembly, so its statistics describe the assembly as produced. A sample with no contig of at least `--min_contig_length` gets an empty file here and is left out of annotation and typing, with a warning in the log.

This step runs for every `--seq_type`.

---

## Genome Annotation and Analysis

### QUAST

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/qc/quast/`
  - `<sample_id>.tsv`: Tab-separated assembly quality metrics; the `Assembly` row is the sample ID.

</details>

[QUAST](https://github.com/ablab/quast) evaluates genome assembly quality by computing metrics such as total length, number of contigs, N50, GC content, and largest contig size. It runs on the unfiltered assembly; the other tools in this section read the [annotation input](#annotation-input) (contigs >= `--min_contig_length`).

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

[Bakta](https://github.com/oschwengers/bakta) provides rapid and standardized annotation of bacterial genomes. It identifies CDSs, tRNAs, rRNAs, ncRNAs, CRISPR arrays, and other genomic features. Bakta runs with `--keep-contig-headers`, so its outputs use the `contig_N` names of the [annotation input](#annotation-input). The Bakta outputs (nucleotide FASTA, protein FASTA, GFF3) are the input of AMRFinderPlus.

### GTDB-Tk

<details markdown="1">
<summary>Output files</summary>

- `gtdbtk/`
  - `gtdbtk.batch.*.summary.tsv`: Taxonomic classification summary for all bacterial and archaeal genomes in the batch.
  - `gtdbtk.batch.*.classify.tree.gz`: Phylogenetic tree with genome placement _(optional)_.
  - `gtdbtk.batch.*.markers_summary.tsv`: Summary of identified marker genes _(optional)_.
  - `gtdbtk.batch.*.msa.fasta.gz`: Multiple sequence alignment of marker genes _(optional)_.
  - `gtdbtk.batch.*.filtered.tsv`: Genomes filtered during classification _(optional)_.
  - `gtdbtk.batch.failed_genomes.tsv`: Genomes that failed classification _(optional)_.
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

[geNomad](https://github.com/apcamargo/genomad) identifies mobile genetic elements (prophages and plasmids) in genome assemblies. It reads the [annotation input](#annotation-input) and classifies contigs as chromosome, plasmid, or virus based on gene content and genomic signatures.

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
- ICEKp-associated virulence loci: yersiniabactin (_ybt_), colibactin (_clb_), salmochelin (_iro_), hypermucoidy (_rmpA_)
- Virulence plasmid associated loci: salmochelin (_iro_), aerobactin (_iuc_), hypermucoidy (_rmpA_, _rmpA2_)
- Antimicrobial resistance determinants: acquired genes, SNPs, gene truncations and intrinsic β-lactamases
- K (capsule) and O antigen (LPS) serotype prediction

**Automatic detection:** Kleborate analysis is automatically triggered when _Klebsiella_ species are detected based on MLST results.

### sccmec

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/annotation/sccmec/`
  - `*.tsv`: SCCmec cassette typing results.

</details>

[sccmec](https://github.com/rpetit3/sccmec) classifies the SCCmec cassette type in _Staphylococcus aureus_ genome assemblies. SCCmec (Staphylococcal Cassette Chromosome _mec_) is a mobile genetic element that carries the _mecA_ gene responsible for methicillin resistance (MRSA).

**Automatic detection:** This analysis is automatically triggered when _S. aureus_ is detected based on MLST results.

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

**Automatic detection:** SISTR analysis is automatically triggered when MLST assigns the `senterica_achtman_2` scheme to a sample (`salmonella` with mlst databases older than 2.23).

---

## Reporting

### Run report

<details markdown="1">
<summary>Output files</summary>

- `report/`
  - `mgap_report.html`: one self-contained page (Plotly is inlined, about 5 MB): it opens offline, from a shared drive or as an e-mail attachment.
  - `mgap_summary.tsv`: one row per sample with its status, read, assembly, quality and typing numbers, AMR genes and warnings; columns no sample has a value for are left out.
  - `mgap_summary.xlsx`: the same summary as the first sheet, then one sheet per tool that ran (sample status, reads, species in the reads, QUAST, length filter, CheckM2, MLST, GTDB-Tk, Kleborate, SISTR, SCCmec, AMRFinderPlus, RGI without its sequence columns, geNomad plasmids and viruses, MOB-suite contigs and plasmids, software versions), each with a `sample` column.

</details>

`MGAP_REPORT` (`bin/mgap_report.py`) reads every tool's output for the run and writes a report for lab staff: each section says what the tool does and how to read its result. Sections appear only for the tools that ran; skipped and failed samples are listed from `pipeline_info/sample_status.tsv`.

| Section                            | What it shows                                                                                                                                                                                                                                                                                                                                 |
| ---------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Run overview                       | Pipeline and Nextflow versions, run date and sequencing type; pass, skipped and failed counts with the samples that did not pass and why; key thresholds; the databases used, with their version where the database records one (Bakta, AMRFinderPlus, geNomad, GTDB-Tk, CARD; Kraken2, CheckM2 and MOB-suite databases show "not recorded"). |
| Read quality and species           | Reads and bases before and after fastp/fastplong with the `--min_bases` limit, Q30, GC, duplication and (Illumina) coverage from Mash's genome size; the species in the reads from Bracken (or Kraken2), with a warning when a second species has more than 5 % of the reads.                                                                 |
| Assembly                           | QUAST total length (with the usual genome length for the MLST species), contigs and N50, and the contigs kept for annotation by `--min_contig_length`.                                                                                                                                                                                        |
| Genome quality (CheckM2)           | Completeness against contamination, with the 90 % / 5 % high-quality limits; genomes with low completeness, high contamination or both have their own colour and marker and are listed first.                                                                                                                                                 |
| Species and typing                 | MLST scheme and ST, the GTDB-Tk species and ANI and the top species in the reads; rows where the genus from MLST differs from GTDB-Tk or the reads are flagged and listed first. Kleborate, SISTR and SCCmec tables for the samples routed to them; SCCmec calls not made from the cassette's target genes are flagged and listed first.      |
| Antimicrobial resistance           | Heatmaps of the AMRFinderPlus acquired genes and the RGI Perfect and Strict hits (gene by sample, colour = identity), the point mutations and every AMRFinderPlus hit.                                                                                                                                                                        |
| Plasmids and other mobile elements | geNomad plasmid and virus sequences per sample and the MOB-suite reconstructed plasmids with replicon, relaxase, mobility and host range.                                                                                                                                                                                                     |
| Samples                            | One card per sample with its key numbers and warnings. Cards start collapsed; a search box and a sample list open them, and samples with warnings are marked (warning sign and red outline) and listed first.                                                                                                                                 |
| Methods, versions and citations    | Every tool version from `pipeline_info/software_versions.yml` and the citations from `CITATIONS.md`, with any tool that has no citation entry.                                                                                                                                                                                                |

With more than 100 samples the one-bar-per-sample charts (read QC, species, assembly, geNomad) become distributions: one box per metric with a point per sample, the outliers (outside 1.5 times the interquartile range or a QC limit) marked, and the most extreme named; long tables scroll inside their box. On copies of the 10 test genomes, a 300-sample report is about 9 MB and opens in about 13 s; real runs with more varied resistance genes have larger heatmaps.

`--skip_report` turns the report off. `accessory_scripts/build_report_manifest.py` builds the report's inputs from an existing results folder, so the report can be rebuilt for a run made before it existed.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: interactive report aggregating fastp/fastplong, Kraken2, QUAST and Bakta results across all samples, with a Software Versions section (every tool except MultiQC itself, whose version is only known after the report is built; it is in `pipeline_info/software_versions.yml`).
  - `multiqc_data/`, `multiqc_plots/`: parsed data and exported plots behind the report.

</details>

[MultiQC](https://multiqc.info) aggregates the per-sample QC results of the run into a single interactive HTML report. Customise it with `--multiqc_config`, `--multiqc_title` and `--multiqc_logo`.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report_<timestamp>.html`, `execution_timeline_<timestamp>.html`, `execution_trace_<timestamp>.txt` and `pipeline_dag_<timestamp>.html`.
  - `software_versions.yml`: Software versions used in the pipeline run (one entry per process, MultiQC included, plus the pipeline and Nextflow versions).
  - `sample_status.tsv`: one row per samplesheet sample, sorted by sample: whether it passed, was skipped or failed, and why (see below).

</details>

`sample_status.tsv` columns (an empty cell is a number the run did not produce, for example read numbers in contig mode):

| Column               | Description                                                                                                                                                                                            |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sample`             | Sample ID from the samplesheet                                                                                                                                                                         |
| `seq_type`           | `illumina`, `ont` or `contig`                                                                                                                                                                          |
| `status`             | `pass` (reached annotation), `skipped` (below `--min_bases`) or `failed`                                                                                                                               |
| `stage`              | Last stage the sample reached; for skipped and failed samples, the stage where it stopped: `read_qc`, `assembly` (SPAdes, Flye, Autocycler), `polishing` (Medaka), `annotation_input`, `annotation`    |
| `reason`             | Why the sample was skipped or failed, and warnings for samples that went on (`QUAST failed`), separated by `;`                                                                                         |
| `input_reads`        | Reads before filtering (fastp/fastplong `before_filtering`; for paired-end, R1 and R2 together)                                                                                                        |
| `input_bases`        | Bases before filtering                                                                                                                                                                                 |
| `input_mean_length`  | Mean read length before filtering (bases ÷ reads)                                                                                                                                                      |
| `input_n50`          | ONT only: read N50 before filtering, approximate (interpolated within the bin of fastplong's read length histogram that holds it; bins are about 12 % wide). fastp reports no read length distribution |
| `filtered_reads`     | Reads after fastp/fastplong                                                                                                                                                                            |
| `filtered_bases`     | Bases after fastp/fastplong: the number `--min_bases` is compared with                                                                                                                                 |
| `assembly_contigs`   | Contigs in the assembly (or the samplesheet FASTA in contig mode), from `SEQKIT_FILTER`'s `seqkit stats`                                                                                               |
| `assembly_length`    | Total assembly length (bp)                                                                                                                                                                             |
| `assembly_n50`       | Assembly N50 (bp)                                                                                                                                                                                      |
| `annotation_contigs` | Contigs of at least `--min_contig_length` in the annotation FASTA                                                                                                                                      |
| `annotation_length`  | Total length of the annotation FASTA (bp)                                                                                                                                                              |

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
