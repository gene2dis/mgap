# MGAP Pipeline Workflow Diagram

## Overview

The **gene2dis/mgap** pipeline supports three input modes based on the `--seq_type` parameter: `illumina`, `ont`, or `contig`. Each mode feeds into a shared downstream analysis that includes assembly quality assessment, a per-sample species check, taxonomic classification, gene annotation, antimicrobial resistance detection, plasmid and insertion sequence (IS element) detection, and a run report.

## Full Pipeline Diagram

```mermaid
flowchart TD
    %% ── Input ──
    INPUT["Samplesheet<br>(--input)"] --> SEQ_TYPE{seq_type?}

    %% ── Illumina subworkflow ──
    SEQ_TYPE -- "illumina" --> FASTP["Fastp<br>(quality trimming)"]
    FASTP --> GATE_ILL{"bases ≥ min_bases?"}
    GATE_ILL -- no --> STATUS
    GATE_ILL -- yes --> KR2_ILL{"run_kraken2 + kraken2db?"}
    KR2_ILL -- yes --> KRAKEN2_ILL["Kraken2"] --> BRACKEN["Bracken<br>(Illumina only)"]
    GATE_ILL -- yes --> COV_ILL{"adjust_coverage?"}
    COV_ILL -- yes --> MASH_ILL["Mash Sketch<br>(genome size;<br>coverage = fastp bases / size)"]
    MASH_ILL --> SEQTK_ILL["Seqtk Sample<br>(subsample if > max_coverage)"]
    SEQTK_ILL --> SPADES
    COV_ILL -- no --> SPADES["SPAdes<br>(assembly)"]

    %% ── ONT subworkflow ──
    SEQ_TYPE -- "ont" --> FASTPLONG["fastplong<br>(QC + adapter trimming)"]
    FASTPLONG --> GATE_ONT{"bases ≥ min_bases?"}
    GATE_ONT -- no --> STATUS
    GATE_ONT -- yes --> KR2_ONT{"run_kraken2 + kraken2db?"}
    KR2_ONT -- yes --> KRAKEN2_ONT["Kraken2<br>(ONT: no Bracken)"]
    GATE_ONT -- yes --> ASM_ONT{"ont_assembler?"}

    %% ── Flye branch ──
    ASM_ONT -- "flye" --> ASMCOV_ONT{"flye_asm_coverage?"}
    ASMCOV_ONT -- "yes, no genome_size" --> LRGE["lrge<br>(genome size)"]
    LRGE --> FLYE
    ASMCOV_ONT -- "no, or genome_size given" --> FLYE["Flye<br>(assembly, all reads)"]
    FLYE --> BEROKKA_CHK{"run_berokka?"}
    BEROKKA_CHK -- yes --> BEROKKA["Berokka<br>(end overlaps of<br>circular contigs)"]
    BEROKKA_CHK -- no --> MEDAKA
    BEROKKA --> MEDAKA["Medaka<br>(polishing)"]
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
    ANNOT_IN --> CHECKM2_CHK{"checkm2_db?"}
    CHECKM2_CHK -- yes --> CHECKM2["CheckM2<br>(completeness / contamination)"]
    ANNOT_IN --> MLST["MLST<br>(sequence typing)"]
    ANNOT_IN --> BAKTA_CHK{"bakta_db?"}
    BAKTA_CHK -- yes --> BAKTA["Bakta<br>(annotation)"]

    %% ── Taxonomy ──
    ANNOT_IN --> GTDBTK_CHK{"run_gtdbtk?"}
    GTDBTK_CHK -- yes --> GTDBTK["GTDB-Tk<br>(taxonomic classification)"]

    %% ── AMR detection ──
    ANNOT_IN --> SKANI_CHK{"species_refs?"}
    SKANI_CHK -- yes --> SKANI["skani search<br>(curated NCBI references)"]
    MLST & SKANI --> SPECIES["SPECIES_CHECK<br>(MLST, skani, reads,<br>samplesheet species)"]
    BRACKEN & KRAKEN2_ONT -. "top species" .-> SPECIES
    SPECIES -- "override / reference" --> MLST_SCHEME["MLST --scheme<br>(reference species' scheme)"]
    SPECIES --> SPECIES_TSV["species_check.tsv"]
    MLST_SCHEME -- "re-run ST" --> SPECIES_TSV
    BAKTA --> AMR_CHK{"amrfinder_db?"}
    AMR_CHK -- yes --> AMR_JOIN["Join Bakta outputs<br>+ species code"]
    SPECIES --> AMR_JOIN
    AMR_JOIN --> AMRFINDER["AMRFinderPlus<br>(AMR detection)"]

    %% ── Genomad ──
    ANNOT_IN --> GENOMAD_CHK{"genomad_db?"}
    GENOMAD_CHK -- yes --> GENOMAD["geNomad<br>(mobile genetic elements)"]

    %% ── RGI (optional) ──
    ANNOT_IN --> RGI_CHK{"run_rgi?"}
    RGI_CHK -- yes --> RGI_LOAD["RGI load<br>(CARD, once per run)"]
    RGI_LOAD --> RGI["RGI<br>(AMR gene prediction)"]

    %% ── Plasmids ──
    ANNOT_IN --> PLASMIDFINDER_CHK{"run_plasmidfinder?"}
    PLASMIDFINDER_CHK -- yes --> PLASMIDFINDER["PlasmidFinder<br>(replicons)"]
    ANNOT_IN --> MOBSUITE_CHK{"run_mobsuite?"}
    MOBSUITE_CHK -- yes --> MOBSUITE["MOB-suite recon<br>(plasmid reconstruction)"]
    MOBSUITE_CHK -- yes --> CIRC_CHK{"seq_type ont / contig?"}
    CIRC_CHK -- yes --> PLASMID_CIRC["Circular contigs<br>(shorter than 1 Mb)"]
    PLASMID_CIRC --> MOBTYPER["MOB-typer<br>(each circular contig)"]
    FLYE -. "assembly_info.txt" .-> PLASMID_CIRC
    FLYE -. "assembly_info.txt" .-> PLASMID_CONS
    PLASMIDFINDER & GENOMAD & MOBSUITE & MOBTYPER --> PLASMID_CONS["Plasmid consensus<br>(one call per contig;<br>when any plasmid tool ran)"]

    %% ── ISEScan (optional) ──
    ANNOT_IN --> ISESCAN_CHK{"run_isescan?"}
    ISESCAN_CHK -- yes --> ISESCAN["ISEScan<br>(insertion sequences)"]

    %% ── Taxa-specific tools ──
    SPECIES --> TAXA_BRANCH{"SPECIES_TYPING<br>Species?"}
    ANNOT_IN --> TAXA_BRANCH
    TAXA_BRANCH -- "K. pneumoniae" --> KLEBORATE["Kleborate"]
    TAXA_BRANCH -- "S. aureus" --> SCCMEC["sccmec"]
    TAXA_BRANCH -- "Salmonella" --> SISTR["SISTR<br>(serotype prediction)"]

    %% ── Software versions ──
    QUAST & ANNOT_IN & CHECKM2 & MLST & MLST_SCHEME & SKANI & SPECIES & BAKTA & AMRFINDER & GENOMAD & GTDBTK --> TOOL_VERSIONS["Tool versions<br>(every process, read QC<br>and assembly steps included)"]
    RGI_LOAD & RGI & PLASMIDFINDER & MOBSUITE & PLASMID_CIRC & MOBTYPER & PLASMID_CONS & ISESCAN --> TOOL_VERSIONS
    KLEBORATE & SCCMEC & SISTR --> TOOL_VERSIONS
    FASTP & FASTPLONG & KRAKEN2_ILL & KRAKEN2_ONT & QUAST & BAKTA --> MULTIQC["MULTIQC<br>(aggregate QC report)"]
    TOOL_VERSIONS -- "Software Versions section" --> MULTIQC
    TOOL_VERSIONS & MULTIQC --> VERSIONS["software_versions.yml<br>(tools + MultiQC + report)"]

    %% ── Sample status ──
    ASSEMBLY & ANNOT_IN -- "pass / failed" --> STATUS["sample_status.tsv<br>(pass, skipped, failed)"]

    %% ── Run report ──
    REPORT_CHK{"skip_report?"}
    STATUS & SPECIES_TSV & VERSIONS --> REPORT_CHK
    TOOL_VERSIONS -. "and each step's result tables<br>(PlasmidFinder, MOB-typer via the consensus)" .-> REPORT_CHK
    REPORT_CHK -- no --> REPORT["MGAP_REPORT<br>(report/mgap_report.html,<br>mgap_summary.tsv / .xlsx)"]

    %% ── Styling ──
    classDef decision fill:#ffd966,stroke:#333,color:#000
    classDef optional fill:#d5e8d4,stroke:#82b366,color:#000
    classDef process fill:#dae8fc,stroke:#6c8ebf,color:#000
    classDef input fill:#f8cecc,stroke:#b85450,color:#000
    classDef merge fill:#e1d5e7,stroke:#9673a6,color:#000

    class SEQ_TYPE,GATE_ILL,GATE_ONT,KR2_ILL,KR2_ONT,COV_ILL,ASMCOV_ONT,ASM_ONT,BEROKKA_CHK,DNAAPLER_FLYE,DNAAPLER_AC,CHECKM2_CHK,BAKTA_CHK,AMR_CHK,GENOMAD_CHK,GTDBTK_CHK,SKANI_CHK,RGI_CHK,PLASMIDFINDER_CHK,MOBSUITE_CHK,CIRC_CHK,ISESCAN_CHK,TAXA_BRANCH,REPORT_CHK decision
    class KRAKEN2_ILL,BRACKEN,KRAKEN2_ONT,MASH_ILL,SEQTK_ILL,LRGE,BEROKKA,DNAAPLER_F,DNAAPLER_G,GFA2FASTA,CHECKM2,BAKTA,AMR_JOIN,AMRFINDER,GENOMAD,GTDBTK,SKANI,MLST_SCHEME,RGI_LOAD,RGI,PLASMIDFINDER,MOBSUITE,PLASMID_CIRC,MOBTYPER,PLASMID_CONS,ISESCAN,REPORT optional
    class FASTP,FASTPLONG,SPADES,FLYE,MEDAKA,AC_GSIZE,AC_SUB,AC_ASM,AC_CONS,QUAST,ANNOT_IN,MLST,KLEBORATE,SCCMEC,SISTR,TOOL_VERSIONS,VERSIONS,STATUS,SPECIES,SPECIES_TSV process
    class INPUT input
    class ASSEMBLY merge
```

## Subworkflow Details

### Illumina Subworkflow

| Step | Tool             | Description                                                                                                 |
| ---- | ---------------- | ----------------------------------------------------------------------------------------------------------- |
| 1    | **Fastp**        | Quality trimming of paired-end reads; samples below `--min_bases` stop here (skipped)                       |
| 2    | **Kraken2**      | Contamination detection _(optional, requires `--kraken2db`)_                                                |
| 2b   | **Bracken**      | Abundance re-estimation, read length (`-r`) matched to the reads (`--brackendb`, else the Kraken2 database) |
| 3    | **Mash Sketch**  | Genome size for the coverage estimate (fastp bases ÷ size) _(optional, if `--adjust_coverage`)_             |
| 4    | **Seqtk Sample** | Read subsampling when coverage exceeds `--max_coverage` (default 150x)                                      |
| 5    | **SPAdes**       | _De novo_ genome assembly                                                                                   |

### ONT Subworkflow (Flye mode, default)

| Step | Tool          | Description                                                                                   |
| ---- | ------------- | --------------------------------------------------------------------------------------------- |
| 1    | **fastplong** | Quality filtering and adapter trimming of long reads                                          |
| 2    | **Kraken2**   | Contamination detection _(optional, requires `--kraken2db`; Bracken not run for ONT)_         |
| 3    | **lrge**      | Genome size for Flye _(only with `--flye_asm_coverage` and no `--genome_size`)_               |
| 4    | **Flye**      | _De novo_ long-read assembly from all reads (`--asm-coverage` with `--flye_asm_coverage`)     |
| 5    | **Berokka**   | Trim the end overlap (tandem copies) of Flye's circular contigs _(optional, `--run_berokka`)_ |
| 6    | **Medaka**    | Assembly polishing                                                                            |
| 7    | **Dnaapler**  | Contig reorientation using FASTA input _(optional, `--run_dnaapler`)_                         |

### ONT Subworkflow (Autocycler mode)

| Step | Tool                       | Description                                                                                         |
| ---- | -------------------------- | --------------------------------------------------------------------------------------------------- |
| 1    | **fastplong**              | Quality filtering and adapter trimming of long reads                                                |
| 2    | **Kraken2**                | Contamination detection _(optional, requires `--kraken2db`; Bracken not run for ONT)_               |
| 3    | **Autocycler genome_size** | Genome size estimation via Raven                                                                    |
| 4    | **Autocycler subsample**   | Subsample reads into independent subsets                                                            |
| 5    | **Autocycler assembly**    | Fan-out assembly across multiple assemblers × subsamples (Plassembler only with `--plassembler_db`) |
| 6    | **Autocycler consensus**   | Compress, cluster, trim/resolve, combine (one task, FASTA + GFA)                                    |
| 7    | **Dnaapler**               | Reorient circular contigs using GFA input _(optional, `run_dnaapler`)_                              |
| 8    | **Autocycler gfa2fasta**   | Convert reoriented GFA back to FASTA _(only when Dnaapler is enabled)_                              |

### Shared Downstream Analysis

| Step | Tool                  | Description                                                                                                                                                                      |
| ---- | --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1    | **QUAST**             | Assembly quality metrics (whole assembly)                                                                                                                                        |
| -    | **SeqKit**            | Annotation FASTA: contigs >= `min_contig_length`, uncompressed, renamed `contig_N`; read by every tool below                                                                     |
| 2    | **CheckM2**           | Genome completeness and contamination assessment _(requires `--checkm2_db`)_                                                                                                     |
| 3    | **MLST**              | Multi-locus sequence typing                                                                                                                                                      |
| 3b   | **skani**             | Closest curated NCBI reference _(optional, `--species_refs`)_                                                                                                                    |
| 3c   | **Species check**     | One species decision per sample (MLST, skani, reads, samplesheet): AMRFinderPlus organism and species tools                                                                      |
| 4    | **Bakta**             | Genome annotation _(requires `--bakta_db`)_                                                                                                                                      |
| 5    | **GTDB-Tk**           | Taxonomic classification _(optional)_                                                                                                                                            |
| 6    | **AMRFinderPlus**     | Antimicrobial resistance gene detection (uses Bakta outputs and the species check's organism) _(requires `--bakta_db` and `--amrfinder_db`)_                                     |
| 7    | **geNomad**           | Mobile genetic element identification _(requires `--genomad_db`)_                                                                                                                |
| 8    | **RGI**               | Resistance gene prediction; CARD loaded once per run _(optional)_                                                                                                                |
| 9    | **MOB-suite**         | Plasmid detection and reconstruction _(optional)_                                                                                                                                |
| 9b   | **PlasmidFinder**     | Plasmid replicon detection _(on by default, `--run_plasmidfinder`)_                                                                                                              |
| 9c   | **MOB-typer**         | Types each circular contig shorter than 1 Mb on its own (ONT and contig mode) _(with `--run_mobsuite`)_                                                                          |
| 9d   | **Plasmid consensus** | One plasmid call per contig from PlasmidFinder, geNomad and MOB-suite (and Flye's `assembly_info.txt`) _(when any of them ran)_                                                  |
| 10   | **ISEScan**           | Insertion sequence (IS element) detection _(optional, `--run_isescan`)_                                                                                                          |
| 11   | **MultiQC**           | Aggregate QC report: fastp / fastplong, Kraken2, QUAST and Bakta, plus software versions                                                                                         |
| 12   | **Run report**        | `report/mgap_report.html`, `mgap_summary.tsv` and `.xlsx` from the tool results, `sample_status.tsv`, `species_check.tsv` and `software_versions.yml` _(unless `--skip_report`)_ |

### Taxa-Specific Analysis

| Species                 | Tool          | Description                           |
| ----------------------- | ------------- | ------------------------------------- |
| _Klebsiella pneumoniae_ | **Kleborate** | Virulence and resistance scoring      |
| _Staphylococcus aureus_ | **sccmec**    | SCCmec cassette typing                |
| _Salmonella enterica_   | **SISTR**     | Serovar and cgMLST subtype prediction |
