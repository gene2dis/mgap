# MGAP Pipeline Workflow Diagram

## Overview

The **gene2dis/mgap** pipeline supports three input modes based on the `--seq_type` parameter: `illumina`, `ont`, or `contig`. Each mode feeds into a shared downstream analysis that includes assembly quality assessment, taxonomic classification, gene annotation, and antimicrobial resistance detection.

## Full Pipeline Diagram

```mermaid
flowchart TD
    %% ── Input ──
    INPUT["Samplesheet<br>(--input)"] --> SEQ_TYPE{seq_type?}

    %% ── Illumina subworkflow ──
    SEQ_TYPE -- "illumina" --> FASTP["Fastp<br>(quality trimming)"]
    FASTP --> KR2_ILL{"run_kraken2?"}
    KR2_ILL -- yes --> KRAKEN2_ILL["Kraken2"] --> BRACKEN["Bracken<br>(Illumina only)"]
    KR2_ILL -- no --> COV_ILL
    BRACKEN --> COV_ILL{"adjust_coverage?"}
    FASTP --> COV_ILL
    COV_ILL -- yes --> MASH_ILL["Mash Sketch<br>(genome size;<br>coverage = fastp bases / size)"]
    MASH_ILL --> SEQTK_ILL["Seqtk Sample<br>(subsample if > max_coverage)"]
    SEQTK_ILL --> SPADES
    COV_ILL -- no --> SPADES["SPAdes<br>(assembly)"]

    %% ── ONT subworkflow ──
    SEQ_TYPE -- "ont" --> FASTPLONG["fastplong<br>(QC + adapter trimming)"]
    FASTPLONG --> KR2_ONT{"run_kraken2?"}
    KR2_ONT -- yes --> KRAKEN2_ONT["Kraken2<br>(ONT: no Bracken)"]
    KR2_ONT -- no --> ASM_ONT
    KRAKEN2_ONT --> ASM_ONT
    FASTPLONG --> ASM_ONT{"ont_assembler?"}

    %% ── Flye branch ──
    ASM_ONT -- "flye" --> ASMCOV_ONT{"flye_asm_coverage?"}
    ASMCOV_ONT -- "yes, no genome_size" --> LRGE["lrge<br>(genome size)"]
    LRGE --> FLYE
    ASMCOV_ONT -- "no, or genome_size given" --> FLYE["Flye<br>(assembly, all reads)"]
    FLYE --> MEDAKA["Medaka<br>(polishing)"]
    MEDAKA --> DNAAPLER_FLYE{"run_dnaapler?"}
    DNAAPLER_FLYE -- yes --> DNAAPLER_F["Dnaapler<br>(reorientation, FASTA)"]
    DNAAPLER_FLYE -- no --> ASSEMBLY
    DNAAPLER_F --> ASSEMBLY

    %% ── Autocycler branch ──
    ASM_ONT -- "autocycler" --> AC_GSIZE["Autocycler<br>(genome size)"]
    AC_GSIZE --> AC_SUB["Autocycler<br>(subsample)"]
    AC_SUB --> AC_ASM["Autocycler<br>(multi-assembler)"]
    AC_ASM --> AC_CONS["Autocycler consensus<br>(compress, cluster,<br>trim + resolve, combine)"]
    AC_CONS --> DNAAPLER_AC{"run_dnaapler?"}
    DNAAPLER_AC -- yes --> DNAAPLER_G["Dnaapler<br>(reorientation, GFA)"]
    DNAAPLER_G --> GFA2FASTA["Autocycler<br>(gfa2fasta)"]
    GFA2FASTA --> ASSEMBLY
    DNAAPLER_AC -- no --> ASSEMBLY

    %% ── Contig input ──
    SEQ_TYPE -- "contig" --> ASSEMBLY

    %% ── Merge into shared assembly channel ──
    SPADES --> ASSEMBLY["Genome Assembly"]

    %% ── Shared downstream analysis ──
    ASSEMBLY --> QUAST["QUAST<br>(assembly QC)"]
    ASSEMBLY --> ANNOT_IN["SeqKit<br>(contigs >= min_contig_length,<br>uncompressed, renamed contig_N)"]
    ANNOT_IN --> CHECKM2["CheckM2<br>(completeness / contamination)"]
    ANNOT_IN --> MLST["MLST<br>(sequence typing)"]
    ANNOT_IN --> BAKTA["Bakta<br>(annotation)"]

    %% ── Taxonomy ──
    ANNOT_IN --> GTDBTK_CHK{"run_gtdbtk?"}
    GTDBTK_CHK -- yes --> GTDBTK["GTDB-Tk<br>(taxonomic classification)"]

    %% ── AMR detection ──
    MLST --> SPECIES["Species Code<br>(MLST → taxa map)"]
    BAKTA --> AMR_JOIN["Join Bakta outputs<br>+ species code"]
    SPECIES --> AMR_JOIN
    AMR_JOIN --> AMRFINDER["AMRFinderPlus<br>(AMR detection)"]

    %% ── Genomad ──
    ANNOT_IN --> GENOMAD["geNomad<br>(mobile genetic elements)"]

    %% ── RGI (optional) ──
    ANNOT_IN --> RGI_CHK{"run_rgi?"}
    RGI_CHK -- yes --> RGI_LOAD["RGI load<br>(CARD, once per run)"]
    RGI_LOAD --> RGI["RGI<br>(AMR gene prediction)"]

    %% ── MOB-suite (optional) ──
    ANNOT_IN --> MOBSUITE_CHK{"run_mobsuite?"}
    MOBSUITE_CHK -- yes --> MOBSUITE["MOB-suite<br>(plasmid reconstruction)"]

    %% ── Taxa-specific tools ──
    SPECIES --> TAXA_BRANCH{"SPECIES_TYPING<br>Species?"}
    ANNOT_IN --> TAXA_BRANCH
    TAXA_BRANCH -- "K. pneumoniae" --> KLEBORATE["Kleborate"]
    TAXA_BRANCH -- "S. aureus" --> SCCMEC["sccmec"]
    TAXA_BRANCH -- "Salmonella" --> SISTR["SISTR<br>(serotype prediction)"]

    %% ── Software versions ──
    QUAST & ANNOT_IN & CHECKM2 & MLST & BAKTA & AMRFINDER & GENOMAD & RGI --> TOOL_VERSIONS["Tool versions"]
    KLEBORATE & SCCMEC & SISTR & MOBSUITE --> TOOL_VERSIONS
    QUAST & BAKTA --> MULTIQC["MULTIQC<br>(aggregate QC report)"]
    TOOL_VERSIONS -- "Software Versions section" --> MULTIQC
    TOOL_VERSIONS & MULTIQC --> VERSIONS["software_versions.yml<br>(tools + MultiQC)"]

    %% ── Styling ──
    classDef decision fill:#ffd966,stroke:#333,color:#000
    classDef optional fill:#d5e8d4,stroke:#82b366,color:#000
    classDef process fill:#dae8fc,stroke:#6c8ebf,color:#000
    classDef input fill:#f8cecc,stroke:#b85450,color:#000
    classDef merge fill:#e1d5e7,stroke:#9673a6,color:#000

    class SEQ_TYPE,KR2_ILL,KR2_ONT,COV_ILL,ASMCOV_ONT,ASM_ONT,DNAAPLER_FLYE,DNAAPLER_AC,GTDBTK_CHK,RGI_CHK,MOBSUITE_CHK,TAXA_BRANCH decision
    class KRAKEN2_ILL,BRACKEN,KRAKEN2_ONT,MASH_ILL,SEQTK_ILL,LRGE,DNAAPLER_F,DNAAPLER_G,GFA2FASTA,GTDBTK,RGI_LOAD,RGI,MOBSUITE optional
    class FASTP,FASTPLONG,SPADES,FLYE,MEDAKA,AC_GSIZE,AC_SUB,AC_ASM,AC_COMP,AC_CLUST,AC_TR,AC_COMB,QUAST,ANNOT_IN,CHECKM2,MLST,BAKTA,AMRFINDER,GENOMAD,KLEBORATE,SCCMEC,SISTR,TOOL_VERSIONS,VERSIONS,SPECIES,AMR_JOIN process
    class INPUT input
    class ASSEMBLY merge
```

## Subworkflow Details

### Illumina Subworkflow

| Step | Tool             | Description                                                                                      |
| ---- | ---------------- | ------------------------------------------------------------------------------------------------ |
| 1    | **Fastp**        | Quality trimming of paired-end reads                                                             |
| 2    | **Kraken2**      | Contamination detection _(optional, requires `--kraken2db`)_                                     |
| 2b   | **Bracken**      | Abundance re-estimation, read length (`-r`) matched to the reads _(requires a Bracken database)_ |
| 3    | **Mash Sketch**  | Genome size for the coverage estimate (fastp bases ÷ size) _(optional, if `--adjust_coverage`)_  |
| 4    | **Seqtk Sample** | Read subsampling when coverage exceeds `--max_coverage` (default 150x)                           |
| 5    | **SPAdes**       | _De novo_ genome assembly                                                                        |

### ONT Subworkflow (Flye mode, default)

| Step | Tool          | Description                                                                               |
| ---- | ------------- | ----------------------------------------------------------------------------------------- |
| 1    | **fastplong** | Quality filtering and adapter trimming of long reads                                      |
| 2    | **Kraken2**   | Contamination detection _(optional, requires `--kraken2db`; Bracken not run for ONT)_     |
| 3    | **lrge**      | Genome size for Flye _(only with `--flye_asm_coverage` and no `--genome_size`)_           |
| 4    | **Flye**      | _De novo_ long-read assembly from all reads (`--asm-coverage` with `--flye_asm_coverage`) |
| 5    | **Medaka**    | Assembly polishing                                                                        |
| 6    | **Dnaapler**  | Contig reorientation using FASTA input _(optional, `--run_dnaapler`)_                     |

### ONT Subworkflow (Autocycler mode)

| Step | Tool                       | Description                                                                           |
| ---- | -------------------------- | ------------------------------------------------------------------------------------- |
| 1    | **fastplong**              | Quality filtering and adapter trimming of long reads                                  |
| 2    | **Kraken2**                | Contamination detection _(optional, requires `--kraken2db`; Bracken not run for ONT)_ |
| 3    | **Autocycler genome_size** | Genome size estimation via Raven                                                      |
| 4    | **Autocycler subsample**   | Subsample reads into independent subsets                                              |
| 5    | **Autocycler assembly**    | Fan-out assembly across multiple assemblers × subsamples                              |
| 6    | **Autocycler consensus**   | Compress, cluster, trim/resolve, combine (one task, FASTA + GFA)                      |
| 7    | **Dnaapler**               | Reorient circular contigs using GFA input _(optional, `run_dnaapler`)_                |
| 8    | **Autocycler gfa2fasta**   | Convert reoriented GFA back to FASTA _(only when Dnaapler is enabled)_                |

### Shared Downstream Analysis

| Step | Tool              | Description                                                                                                  |
| ---- | ----------------- | ------------------------------------------------------------------------------------------------------------ |
| 1    | **QUAST**         | Assembly quality metrics (whole assembly)                                                                    |
| -    | **SeqKit**        | Annotation FASTA: contigs >= `min_contig_length`, uncompressed, renamed `contig_N`; read by every tool below |
| 2    | **CheckM2**       | Genome completeness and contamination assessment                                                             |
| 3    | **MLST**          | Multi-locus sequence typing                                                                                  |
| 4    | **Bakta**         | Genome annotation                                                                                            |
| 5    | **GTDB-Tk**       | Taxonomic classification _(optional)_                                                                        |
| 6    | **AMRFinderPlus** | Antimicrobial resistance gene detection (uses Bakta + MLST outputs)                                          |
| 7    | **geNomad**       | Mobile genetic element identification                                                                        |
| 8    | **RGI**           | Resistance gene prediction; CARD loaded once per run _(optional)_                                            |
| 9    | **MOB-suite**     | Plasmid detection and reconstruction _(optional)_                                                            |

### Taxa-Specific Analysis

| Species                 | Tool          | Description                           |
| ----------------------- | ------------- | ------------------------------------- |
| _Klebsiella pneumoniae_ | **Kleborate** | Virulence and resistance scoring      |
| _Staphylococcus aureus_ | **sccmec**    | SCCmec cassette typing                |
| _Salmonella enterica_   | **SISTR**     | Serovar and cgMLST subtype prediction |
