# MGAP Pipeline Workflow Diagram

## Overview

The **gene2dis/mgap** pipeline supports three input modes based on the `--seq_type` parameter: `illumina`, `ont`, or `contig`. Each mode feeds into a shared downstream analysis that includes assembly quality assessment, taxonomic classification, gene annotation, and antimicrobial resistance detection.

## Full Pipeline Diagram

```mermaid
flowchart TD
    %% ── Input ──
    INPUT["Samplesheet\n(--input)"] --> SEQ_TYPE{seq_type?}

    %% ── Illumina subworkflow ──
    SEQ_TYPE -- "illumina" --> FASTP["Fastp\n(quality trimming)"]
    FASTP --> KR2_ILL{"run_kraken2?"}
    KR2_ILL -- yes --> KRAKEN2_ILL["Kraken2"] --> BRACKEN["Bracken"]
    KR2_ILL -- no --> COV_ILL
    BRACKEN --> COV_ILL{"adjust_coverage?"}
    FASTP --> COV_ILL
    COV_ILL -- yes --> MASH_ILL["Mash Sketch\n(estimate coverage)"]
    MASH_ILL --> SEQTK_ILL["Seqtk Sample\n(subsample if needed)"]
    SEQTK_ILL --> SPADES
    COV_ILL -- no --> SPADES["SPAdes\n(assembly)"]

    %% ── ONT subworkflow ──
    SEQ_TYPE -- "ont" --> FASTPLONG["fastplong\n(QC + adapter trimming)"]
    FASTPLONG --> KR2_ONT{"run_kraken2?"}
    KR2_ONT -- yes --> KRAKEN2_ONT["Kraken2"]
    KR2_ONT -- no --> ASM_ONT
    KRAKEN2_ONT --> ASM_ONT
    FASTPLONG --> ASM_ONT{"ont_assembler?"}

    %% ── Flye branch ──
    ASM_ONT -- "flye" --> COV_ONT{"adjust_coverage?"}
    COV_ONT -- yes --> MASH_ONT["Mash Sketch\n(estimate coverage)"]
    MASH_ONT --> SEQTK_ONT["Seqtk Sample\n(subsample if needed)"]
    SEQTK_ONT --> FLYE
    COV_ONT -- no --> FLYE["Flye\n(assembly)"]
    FLYE --> MEDAKA["Medaka\n(polishing)"]
    MEDAKA --> DNAAPLER_FLYE{"run_dnaapler?"}
    DNAAPLER_FLYE -- yes --> DNAAPLER_F["Dnaapler\n(reorientation, FASTA)"]
    DNAAPLER_FLYE -- no --> ASSEMBLY
    DNAAPLER_F --> ASSEMBLY

    %% ── Autocycler branch ──
    ASM_ONT -- "autocycler" --> AC_GSIZE["Autocycler\n(genome size)"]
    AC_GSIZE --> AC_SUB["Autocycler\n(subsample)"]
    AC_SUB --> AC_ASM["Autocycler\n(multi-assembler)"]
    AC_ASM --> AC_COMP["Autocycler\n(compress)"]
    AC_COMP --> AC_CLUST["Autocycler\n(cluster)"]
    AC_CLUST --> AC_TR["Autocycler\n(trim + resolve)"]
    AC_TR --> AC_COMB["Autocycler\n(combine)"]
    AC_COMB --> DNAAPLER_AC{"run_dnaapler?"}
    DNAAPLER_AC -- yes --> DNAAPLER_G["Dnaapler\n(reorientation, GFA)"]
    DNAAPLER_G --> GFA2FASTA["Autocycler\n(gfa2fasta)"]
    GFA2FASTA --> ASSEMBLY
    DNAAPLER_AC -- no --> ASSEMBLY

    %% ── Contig input ──
    SEQ_TYPE -- "contig" --> ASSEMBLY

    %% ── Merge into shared assembly channel ──
    SPADES --> ASSEMBLY["Genome Assembly"]

    %% ── Shared downstream analysis ──
    ASSEMBLY --> QUAST["QUAST\n(assembly QC)"]
    ASSEMBLY --> CHECKM2["CheckM2\n(completeness / contamination)"]
    ASSEMBLY --> MLST["MLST\n(sequence typing)"]
    ASSEMBLY --> BAKTA["Bakta\n(annotation)"]

    %% ── Taxonomy ──
    ASSEMBLY --> GTDBTK_CHK{"run_gtdbtk?"}
    GTDBTK_CHK -- yes --> GTDBTK["GTDB-Tk\n(taxonomic classification)"]

    %% ── AMR detection ──
    MLST --> SPECIES["Species Code\n(MLST → taxa map)"]
    BAKTA --> AMR_JOIN["Join Bakta outputs\n+ species code"]
    SPECIES --> AMR_JOIN
    AMR_JOIN --> AMRFINDER["AMRFinderPlus\n(AMR detection)"]

    %% ── Genomad ──
    BAKTA --> GENOMAD["geNomad\n(mobile genetic elements)"]

    %% ── RGI (optional) ──
    ASSEMBLY --> RGI_CHK{"run_rgi?"}
    RGI_CHK -- yes --> RGI["RGI\n(AMR gene prediction)"]

    %% ── Taxa-specific tools ──
    SPECIES --> TAXA_BRANCH{"Species?"}
    BAKTA --> TAXA_BRANCH
    TAXA_BRANCH -- "K. pneumoniae" --> KLEBORATE["Kleborate"]
    TAXA_BRANCH -- "S. aureus" --> STAPHOPIASCCMEC["staphopia-sccmec"]

    %% ── Software versions ──
    QUAST & CHECKM2 & MLST & BAKTA & AMRFINDER & GENOMAD --> VERSIONS["CUSTOM_DUMPSOFTWAREVERSIONS\n(collect versions)"]
    KLEBORATE & STAPHOPIASCCMEC --> VERSIONS

    %% ── Styling ──
    classDef decision fill:#ffd966,stroke:#333,color:#000
    classDef optional fill:#d5e8d4,stroke:#82b366,color:#000
    classDef process fill:#dae8fc,stroke:#6c8ebf,color:#000
    classDef input fill:#f8cecc,stroke:#b85450,color:#000
    classDef merge fill:#e1d5e7,stroke:#9673a6,color:#000

    class SEQ_TYPE,KR2_ILL,KR2_ONT,COV_ILL,COV_ONT,ASM_ONT,DNAAPLER_FLYE,DNAAPLER_AC,GTDBTK_CHK,RGI_CHK,TAXA_BRANCH decision
    class KRAKEN2_ILL,BRACKEN,KRAKEN2_ONT,MASH_ILL,SEQTK_ILL,MASH_ONT,SEQTK_ONT,DNAAPLER_F,DNAAPLER_G,GFA2FASTA,GTDBTK,RGI optional
    class FASTP,FASTPLONG,SPADES,FLYE,MEDAKA,AC_GSIZE,AC_SUB,AC_ASM,AC_COMP,AC_CLUST,AC_TR,AC_COMB,QUAST,CHECKM2,MLST,BAKTA,AMRFINDER,GENOMAD,KLEBORATE,STAPHOPIASCCMEC,VERSIONS,SPECIES,AMR_JOIN process
    class INPUT input
    class ASSEMBLY merge
```

## Subworkflow Details

### Illumina Subworkflow

| Step | Tool | Description |
|------|------|-------------|
| 1 | **Fastp** | Quality trimming of paired-end reads |
| 2 | **Kraken2 + Bracken** | Contamination detection *(optional)* |
| 3 | **Mash Sketch** | Coverage estimation *(optional, if `adjust_coverage`)* |
| 4 | **Seqtk Sample** | Read subsampling when coverage exceeds `max_coverage` |
| 5 | **SPAdes** | *De novo* genome assembly |

### ONT Subworkflow (Flye mode, default)

| Step | Tool | Description |
|------|------|-------------|
| 1 | **fastplong** | Quality filtering and adapter trimming of long reads |
| 2 | **Kraken2** | Contamination detection *(optional)* |
| 3 | **Mash Sketch** | Coverage estimation *(optional, if `adjust_coverage`)* |
| 4 | **Seqtk Sample** | Read subsampling when coverage exceeds `max_coverage` |
| 5 | **Flye** | *De novo* long-read assembly |
| 6 | **Medaka** | Assembly polishing |
| 7 | **Dnaapler** | Contig reorientation using FASTA input *(optional, `run_dnaapler`)* |

### ONT Subworkflow (Autocycler mode)

| Step | Tool | Description |
|------|------|-------------|
| 1 | **fastplong** | Quality filtering and adapter trimming of long reads |
| 2 | **Kraken2** | Contamination detection *(optional)* |
| 3 | **Autocycler genome_size** | Genome size estimation via Raven |
| 4 | **Autocycler subsample** | Subsample reads into independent subsets |
| 5 | **Autocycler assembly** | Fan-out assembly across multiple assemblers × subsamples |
| 6 | **Autocycler compress** | Compress assemblies into unitig graph |
| 7 | **Autocycler cluster** | Cluster unitigs into putative genomic sequences |
| 8 | **Autocycler trim + resolve** | Trim and resolve each QC-pass cluster |
| 9 | **Autocycler combine** | Combine resolved clusters into consensus assembly (FASTA + GFA) |
| 10 | **Dnaapler** | Reorient circular contigs using GFA input *(optional, `run_dnaapler`)* |
| 11 | **Autocycler gfa2fasta** | Convert reoriented GFA back to FASTA *(only when Dnaapler is enabled)* |

### Shared Downstream Analysis

| Step | Tool | Description |
|------|------|-------------|
| 1 | **QUAST** | Assembly quality metrics |
| 2 | **CheckM2** | Genome completeness and contamination assessment |
| 3 | **MLST** | Multi-locus sequence typing |
| 4 | **Bakta** | Genome annotation |
| 5 | **GTDB-Tk** | Taxonomic classification *(optional)* |
| 6 | **AMRFinderPlus** | Antimicrobial resistance gene detection (uses Bakta + MLST outputs) |
| 7 | **geNomad** | Mobile genetic element identification |
| 8 | **RGI** | Resistance gene prediction *(optional)* |

### Taxa-Specific Analysis

| Species | Tool | Description |
|---------|------|-------------|
| *Klebsiella pneumoniae* | **Kleborate** | Virulence and resistance scoring |
| *Staphylococcus aureus* | **staphopia-sccmec** | SCCmec cassette typing |
