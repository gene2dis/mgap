# Tool updates benchmark (2026-10-29)

Benchmark for the 2026 review items U1 (Medaka 1.11.3 → 2.x) and U2 (newer builds of the other pinned tools); the
review and roadmap are kept outside the repository, and "Phase 10" below refers to its benchmark and apply steps.
Nothing in the pipeline changed for this benchmark: the baseline is a `dev` run (6586280), and every candidate ran as a
standalone container call with the module's own command line, or as a `-c` override of one process, listed in
[Reproducing](#reproducing). A version is re-pinned only after you approve it, one tool per commit.

Setup: Nextflow 26.04.4, `-profile docker`, executor capped at 32 CPUs / 600 GB, on a shared 80-core server (load
average 35–70 during these runs). Timings are comparable side by side, not with other machines. Each candidate image
was pulled, its `--version` and the `--help` of every subcommand the module calls were diffed against the pinned image,
and its Singularity image (`https://depot.galaxyproject.org/singularity/<tool>:<tag>`) answered HTTP 200. Release notes
are from each project's GitHub releases; tags from quay.io, newest first, on 2026-10-29.

Data (in-house isolates, not distributed): the ONT fixture SCL10571 at ~100x (K. pneumoniae ST11, dorado
`dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, 11 contigs after Flye and Berokka), read headers of the other two ONT fixtures
(SCL25346, MinKNOW `sup-v4.1.0`; SCL4611, a 2021 R10.4 model), and the Illumina routing set (10 samples) from the
Phase 8a results directory (`dev` b11dd97): its annotation FASTAs and Bakta outputs are the inputs of the U2
comparisons. There is no Illumina run of SCL10571, so Medaka's consensus is compared with an independent long-read
consensus instead (the Autocycler 0.8.0 consensus of the same reads from Phase 9b, made without Medaka).

## Summary

| Tool (module)                            | Pinned           | Candidate                          | Output on the fixtures                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | Recommendation                                                                                                                                                                                          |
| ---------------------------------------- | ---------------- | ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Medaka (`MEDAKA`)                        | 1.11.3           | 2.2.2                              | With the same model, 2.2.2's consensus is **byte-identical** to 1.11.3's; 2.2.2 is 2.5× faster (pipeline task 15 min 42 s → 6 min 13 s, peak RSS 8.6 → 7.3 GB). Every difference comes from the model: 1.11.3 cannot read the basecaller model from any of our read headers and always uses its default `r1041_e82_400bps_sup_v4.3.0`; 2.2.2 reads SCL10571's model and picks `sup_v5.0.0`, or with `--bacteria` the bacterial methylation model. On SCL10571 the models disagree at 35–73 sites, 29–67 of them in one 22 kb linear contig, a divergent second copy of a plasmid region; against the Autocycler consensus (99.1 % of the assembly compared) 1.11.3 and `--bacteria` differ at 2 sites, `sup_v5.0.0` at 6. | Update, but only with the models inside the image (see [options](#options-for-medaka)): the 2.2.2 biocontainer bundles 4 models and downloads the rest at run time, which fails inside a pipeline task. |
| RGI (`RGI_LOAD`, `RGI_MAIN`)             | 6.0.3            | 6.0.8 (newest)                     | **Loses Perfect hits**: CTX-M-65 (100.1) and AAC(6')-Ib-cr5 (SCL2333, SCL2461, SCL2535), which AMRFinderPlus finds on the same contigs. The regression starts in 6.0.5.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | **6.0.4**, not 6.0.8: it keeps every 6.0.3 hit in all 10 samples and adds 25 point-mutation hits that 6.0.3 cannot parse (most confirmed by AMRFinderPlus).                                             |
| Dnaapler (`DNAAPLER`)                    | 1.3.0            | 1.4.0                              | Byte-identical reoriented FASTA and summary (SCL10571).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | Update (the Autocycler image already has 1.4.0).                                                                                                                                                        |
| SCCmec (`SCCMEC`)                        | 1.2.0            | —                                  | No newer release or build.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | Nothing to do.                                                                                                                                                                                          |
| Kleborate (`KLEBORATE`)                  | 3.2.4 (`_0`)     | 3.2.4 `_1` rebuild; 3.3.0b1 (beta) | `_1` (Kaptive 3.1.0 → 3.2.2, minimap2 2.30 → 2.31, Python 3.14): every field identical on SCL2333, SCL2535, SCL2461.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | Optional rebuild; skip the beta.                                                                                                                                                                        |
| sistr_cmd (`SISTR`)                      | 1.1.3 `_2`       | —                                  | The pinned build is the newest.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | Nothing to do.                                                                                                                                                                                          |
| AMRFinderPlus (`AMRFINDERPLUS_RUN`)      | 4.2.5            | 4.2.7                              | All 18 tables (10 samples, results and mutations) byte-identical; same database version and organism list.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | Update (bug fix only: a crash when an HMM hit has no BLAST alignment).                                                                                                                                  |
| Mash (`MASH_SKETCH`)                     | 2.3 `he348c14_1` | 2.3 `hf85e966_11` rebuild          | Sketches and stats byte-identical (SCL2535, 100.1).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | Nothing to gain; leave.                                                                                                                                                                                 |
| GTDB-Tk (`GTDBTK`, nf-core)              | 2.5.2            | 2.8.0                              | 2.5.2 cannot read the only database on the server (GTDB r232 needs ≥ 2.7.0). 2.8.0 classified the routing set in Phase 8a (1 h 20 min, peak RSS 176.7 GB of 200 GB).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | Update to 2.8.0 (a fix, not an option): patch the nf-core module's container.                                                                                                                           |
| PlasmidFinder (`PLASMIDFINDER`, nf-core) | 2.1.6 `py310_1`  | 2.1.6 `py314_2` rebuild; 3.0.3     | Rebuild: same 2019 database (identical checksum), identical tables on all 11 Phase 9a fixtures. 3.0.3 is a rewrite (`plasmidfinder` command, BeOne JSON output) and bundles **no** database.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | Stay on 2.1.6. Whether a newer database changes the calls cannot be tested without one: supply a PlasmidFinder database if you want that checked.                                                       |

### Decisions needed

1. Medaka: option A (custom image with the models, pushed by you to Docker Hub), B (stay on 1.11.3) or C (2.2.2
   biocontainer as is, `--medaka_model` required for anything but `sup_v5.2.0` data); and, with A, whether
   `--bacteria` is on by default (it fails samples whose reads do not name a compatible basecaller model, see below).
2. RGI 6.0.4 (recommended), 6.0.8 (not recommended) or stay on 6.0.3.
3. Dnaapler 1.4.0, AMRFinderPlus 4.2.7: update or not.
4. Kleborate 3.2.4 `_1`: update or not.
5. GTDB-Tk 2.8.0 by patching the container (recommended), or `nf-core modules update` (2.7.2 Seqera image, version
   topics) plus a patch to 2.8.0.
6. PlasmidFinder: stay on 2.1.6 (recommended); whether you have a newer database to test.

## U1 · Medaka

### How the model is chosen

`MEDAKA` runs `medaka_consensus` with `-m <model>` only when `--medaka_model` is set (it is `null` by default). Without
`-m`, both versions try to read the basecaller model from the reads (`medaka tools resolve_model --auto_model`) and fall
back to their default model with a warning.

| Read header (fixture)                                                                     | Medaka 1.11.3                                                                                                | Medaka 2.2.2                                                                                   | Medaka 2.2.2 `--bacteria`                                             |
| ----------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| dorado `RG:Z:<run>_dna_r10.4.1_e8.2_400bps_sup@v5.0.0_<kit>_<barcode>` (SCL10571)         | not parsed ("did not contain precisely 1 basecaller model reference"): default `r1041_e82_400bps_sup_v4.3.0` | `r1041_e82_400bps_sup_v5.0.0`                                                                  | `r1041_e82_400bps_bacterial_methylation`                              |
| MinKNOW `model_version_id=dna_r10.4.1_e8.2_400bps_sup-v4.1.0` (SCL25346)                  | not parsed: default `sup_v4.3.0`                                                                             | "Unknown basecaller model" (medaka knows `sup@v4.1.0`, not `sup-v4.1.0`): default `sup_v5.2.0` | exits 1: "--bacteria was specified but input data was not compatible" |
| MinKNOW `model_version_id=2021-11-17_dna_r10.4_minion_promethion_1024_67af0493` (SCL4611) | not parsed: default `sup_v4.3.0`                                                                             | unknown: default `sup_v5.2.0`                                                                  | exits 1                                                               |

So today every MGAP sample is polished with `r1041_e82_400bps_sup_v4.3.0` unless `--medaka_model` is set, whatever
basecaller made the reads; the dev run's log on SCL10571 says so ("WARNING: Failed to detect a model version, will use
default: 'r1041_e82_400bps_sup_v4.3.0'"). 1.11.3 has no v5 models at all.

In 2.2.2 (2.1.1 onward), `--bacteria` replaces the detected model by `r1041_e82_400bps_bacterial_methylation` when the
detected one is a 400 bps hac/sup v4.2.0, v4.3.0, v5.0.0, v5.2.0 or hac v6.0.0 model, and otherwise stops with an
error instead of falling back. `-m` switches off both the detection and `--bacteria`.

**Models are not all in the image.** The 2.2.2 biocontainer bundles four model files (`sup_v5.2.0`, `hac_v6.0.0` and two
variant models); any other model, `sup_v5.0.0` and the bacterial model included, is downloaded from GitHub on first use
into `site-packages/medaka/data` or `~/.medaka/data`. A pipeline task runs as `-u <uid>:<gid>` with `HOME=/`, so neither
is writable: run as Nextflow runs it, `resolve_model` fails (and Medaka's own error message then crashes on a
`PosixPath`). Inside `medaka_consensus` that failure is silenced: without `--bacteria` SCL10571 would quietly get the
default `sup_v5.2.0` instead of its `sup_v5.0.0`, and with `--bacteria` the task exits 1 and the sample is recorded as
failed at polishing (D-7). The same applies to an offline node and to a read-only Singularity image. The model files are
small (about 1.5 MB each; all 25 consensus models that automatic selection can choose plus the bacterial one add 44 MB
to a 2.34 GB image).

### What `--medaka_model` means after an update

- 1.11.3: a medaka model name from 1.11.3's list (newest `r1041_e82_400bps_sup_v4.3.0`); unset means
  `sup_v4.3.0` for every sample.
- 2.2.2: a medaka model name from 2.2.2's list (adds the v5.0.0, v5.2.0 and v6.0.0 models and the bacterial one, drops
  the old r9.4.1 `high_g*` models), a basecaller name with a variety (`dna_r10.4.1_e8.2_400bps_sup@v5.0.0:consensus`),
  or a path to a model file. Unset means the model named in the reads (dorado `RG:Z:` or MinKNOW
  `basecall_model_version_id=`), else `sup_v5.2.0`. When set it also disables `--bacteria`.
- In both, a value from the other version's list may not exist (for example `r941_min_high_g360` is gone in 2.x,
  `r1041_e82_400bps_sup_v5.0.0` is missing in 1.11.3) and stops the task.

### Runs on SCL10571

All four polished the same input: the Flye + Berokka assembly and fastplong reads from the `dev` MEDAKA task, with
`medaka_consensus -t 12` in a 12-CPU / 72 GB container as `process_high` gives the task. Run A reproduces the `dev`
pipeline's polished FASTA exactly.

| Run | Image  | Model                                   |   Wall time |   CPU | Consensus                                  |
| --- | ------ | --------------------------------------- | ----------: | ----: | ------------------------------------------ |
| A   | 1.11.3 | default `sup_v4.3.0` (detection failed) | 15 min 12 s | 184 % | 11 contigs, 5,986,182 bp; = `dev` pipeline |
| B   | 2.2.2  | `-m r1041_e82_400bps_sup_v4.3.0`        |  6 min 49 s | 300 % | **byte-identical to A**                    |
| C   | 2.2.2  | automatic: `sup_v5.0.0`                 |  6 min 19 s | 284 % | 11 contigs, 5,986,201 bp                   |
| D   | 2.2.2  | `--bacteria`: `bacterial_methylation`   |   7 min 8 s | 266 % | 11 contigs, 5,986,260 bp                   |

Peak memory: the `dev` trace gives 8.6 GB for A; GNU `time` in these images over-reports by about 4× (36.5 GB for the
same run), so it is used only as a ratio: B–D used 11–16 % less than A. In the pipeline (run D's settings through a `-c`
override and a scratch image with the models, resumed on the `dev` work directory so Flye and Berokka were reused), the
MEDAKA task took 6 min 13 s and 7.3 GB peak RSS (305 % CPU), against 15 min 42 s and 8.6 GB in `dev`, and its output
was identical to run D. The pipeline then completed normally (Dnaapler, annotation FASTA, report).

### Consensus differences

MUMmer `dnadiff` (one-to-one alignments; consecutive indel bases merged into one event). Homopolymer = length of the run
of the inserted/deleted base at the site; motif = the event touches a Dam (`GATC`) or Dcm (`CCWGG`) site, the motifs the
bacterial model is trained for.

| Comparison                            | Aligned | Events | SNPs | Insertions | Deletions | In homopolymer ≥ 4 | At GATC / CCWGG | Of which in contig_3 |
| ------------------------------------- | ------: | -----: | ---: | ---------: | --------: | -----------------: | --------------: | -------------------: |
| Flye draft → A (what 1.11.3 changed)  |   100 % |     13 |    6 |          4 |         3 |                  1 |               0 |                    0 |
| A vs C (`sup_v4.3.0` vs `sup_v5.0.0`) |   100 % |     35 |    2 |         22 |        11 |                  4 |         1 (Dcm) |                   29 |
| A vs D (`sup_v4.3.0` vs bacterial)    |   100 % |     66 |    6 |         57 |         3 |                  6 |               0 |                   64 |
| C vs D                                |   100 % |     73 |   13 |         52 |         8 |                  2 |         1 (Dcm) |                   67 |
| Autocycler vs A                       |  99.1 % |      2 |    1 |          0 |         1 |                  0 |               0 |                    — |
| Autocycler vs C                       |  99.1 % |      6 |    0 |          1 |         5 |                  1 |               0 |                    — |
| Autocycler vs D                       |  99.1 % |      2 |    1 |          0 |         1 |                  0 |               0 |                    — |

contig_3 is a 21,988 bp linear contig that Flye puts in one alternative group with contig_2 (`assembly_info`: repeat,
graph path `*,3,*`); it matches 21 kb of the 189.6 kb plasmid at 99.2 % identity, a divergent second copy whose reads are
hard to place, and that is where the models disagree. Outside it the three models differ at 2–6 sites. "Aligned" is the
part of the second assembly in one-to-one alignments: against Autocycler, contig_2 (a 99.99 % copy of part of the 110 kb
plasmid), contig_3 and the 4.2 kb plasmid contig_10 (absent from Autocycler) drop out, 51 kb or 0.9 % of the assembly. The reads are
dorado sup v5.0.0 at ~100x and the Flye draft is already within 2 SNPs of the Autocycler consensus, so this fixture can
show that a model change is small, not which model is more accurate; none of the differences is at a Dam or Dcm site.
The bacterial model is aimed at errors at methylated motifs in native DNA, which this basecaller already handles.

Downstream, the three consensuses (A, C, D) went through MGAP in contig mode as three samples (`dev`, Bakta and
AMRFinderPlus databases of section 7; 29 min):

| Consensus              | MLST                 |                   AMRFinderPlus rows | Bakta CDS | Bakta pseudogenes |
| ---------------------- | -------------------- | -----------------------------------: | --------: | ----------------: |
| A, 1.11.3 `sup_v4.3.0` | ST11, same 7 alleles |                                   72 |     5,681 |                14 |
| C, 2.2.2 `sup_v5.0.0`  | ST11, same 7 alleles | 72 (same genes, contigs, identities) |     5,683 |                17 |
| D, 2.2.2 bacterial     | ST11, same 7 alleles |                            72 (same) |     5,679 |                15 |

The Bakta differences cannot be pinned on Medaka: Medaka writes the contigs in a different order on every run (also two
1.11.3 runs on the same input), `SEQKIT_FILTER` numbers them `contig_1..N` in that order, and Bakta's gene caller trains
on the sequences in input order. Compared by original contig name, A and C differ in 12 chromosome CDS lengths where the
two consensuses differ at only 2 sites. Plasmid gene calls in contig_3 (where the models disagree) do change.

Found on the way: because of that order, an ONT sample's annotation contig numbers (`contig_N`) change between reruns of
the same data; within one run every tool still reports the same IDs (D-4).

### Options for Medaka

- **A. Medaka 2.2.2 in a small custom image with the models (recommended).** `docker/medaka/Dockerfile`: `FROM` the
  biocontainer by digest, then `medaka tools download_models` for the 25 consensus models automatic selection can
  choose plus the bacterial model, published like the Autocycler image (you push it to Docker Hub; config pins the
  digest; Singularity pulls it with `docker://`). Tested here as a local scratch image (`HOME=/`, non-root, no network:
  the bacterial model loads from the image). `--medaka_model` keeps its role (now with 2.x names); the default
  becomes "the model the reads name". A new `--medaka_bacteria` (default off) adds `--bacteria`; on by default would
  fail SCL25346 and SCL4611 at polishing (their headers name no compatible model) unless `--medaka_model` is set, and on
  SCL10571 it gave no measurable gain.
- **B. Stay on 1.11.3.** Correct but slow, and `--medaka_model` must be set by hand to match anything but sup v4.3.0.
- **C. 2.2.2 biocontainer as is.** Works offline only for `sup_v5.2.0` (the bundled default): v5.0.0 data silently get
  `sup_v5.2.0`, `--bacteria` fails. Not recommended.

Either way, the docs should say that `--medaka_model` should be set for reads whose header Medaka cannot read (MinKNOW
`sup-v4.1.0`-style names, old R10.4).

## U2 · Other tools

Output comparisons used the module's own command line and inputs: for RGI the CARD database at `/mds_data/dbs/rgi`
(CARD 4.0.1) loaded with each version's `rgi load` + `build_databases()` as `RGI_LOAD` does, then `rgi main -t contig
-a DIAMOND -d wgs` on each annotation FASTA, 6 threads as `process_medium`; for AMRFinderPlus the database
`/mds_data/dbs/amrfinder/data/latest` (2025-12-03.1) with the Phase 8a Bakta outputs and the organism `getTaxaNames()`
gives each sample. The pinned versions reproduced the Phase 8a published tables exactly (RGI 6.0.3 all 10 samples,
AMRFinderPlus 4.2.5 checked on 3).

### RGI

| Version                         | Release notes (affecting outputs)                                                                                   | vs 6.0.3 on 10 samples                                                                                                                                                                                           |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 6.0.4 (2025-02-26)              | newer BLAST, Pyrodigal, KMA, samtools; support for `Var` and `Ter` mutations                                        | every 6.0.3 hit kept (same model, cut-off, identity); +25 point-mutation hits; no errors                                                                                                                         |
| 6.0.5 (2025-05-01)              | `AST_Source` column for mutations; "fixed `contig` column"                                                          | the same +25, but CTX-M-65 (100.1) and AAC(6')-Ib-cr5 (SCL2333, SCL2461, SCL2535) lost                                                                                                                           |
| 6.0.6–6.0.8 (2025-11 → 2026-05) | `split_prodigal_jobs` removed; longer Prodigal ORFs; newer DIAMOND; `--num_threads` renamed `--threads` (`-n` kept) | as 6.0.5, plus: H-NS (Salmonella, Kpn), acrB and soxR (Salmonella), emrA (511.1), QnrB10 (SCL2333), adeF (41.7 %), HelR (CCG_01) lost; tet(38) Perfect found in 2 S. aureus; identities move by up to 2.9 points |

- **6.0.3 drops variant models.** Every 6.0.3 run logs `invalid literal for int() with base 10: '1Va' -> model_type:
variant` once: CARD 4.0.1 has mutation types 6.0.3 cannot parse, so its point-mutation search stops. 6.0.4 adds 25
  hits, all Strict, and AMRFinderPlus (organism mode) reports the same mutation for most: gyrA D87Y (100.1); gyrA S84L,
  parC S80F/S80Y/E84G (the three S. aureus); murA G257D (SCL2273); gyrA D87G and parC S80I (SCL2333); gyrA D87N,
  parC S80I and uhpT E350Q (SCL2461). RGI-only: Haemophilus influenzae PBP3 D350N/S357N at 52 % identity (both
  Salmonella and three Klebsiella), E. coli EF-Tu R234F (pulvomycin; 511.1, SCL2333, SCL2461 ×2) and UhpT E350Q in
  the Salmonella and two Klebsiella that AMRFinderPlus does not flag.
- **6.0.5–6.0.8 lose real genes.** In the four samples that lose CTX-M-65 or AAC(6')-Ib-cr5 (each on a contig of 1.2–2.4 kb) RGI logs `Requested rname contig does not exist! Please check your
FASTA file` from `nudge_strict_to_perfect`, whose exception handler only logs the error, so the ORF's hits are
  dropped. The new `Contig` column is the contig name (`contig_28`) instead of the ORF ID (`contig_28_37`).
- Module: 6.0.4 needs only the container and version strings; 6.0.8 would also need `--num_threads` → `-n`
  (`rgi main` exits 2 on `--num_threads`). `RGI_LOAD`'s `app.Database.build_databases()` works in both.
- Images: `quay.io/biocontainers/rgi:6.0.4--pyh05cac1d_1` (sha256:a1e6ca20…, `rgi main --version` 6.0.4, depot 200).

### Dnaapler

1.4.0 (2026-06-30): reads gzip/bzip2/xz input; fixes reorientation to the wrong coordinate when the chosen CDS is on
the negative strand in the `largest`, `mystery` and `nearest` fallbacks (#102); GFA input also writes a reoriented
FASTA; GFA without circular contigs warns instead of failing. `dnaapler all -h` unchanged. On SCL10571's polished
assembly (`dnaapler all -t 6`, as the module) the reoriented FASTA and the summary are byte-identical to 1.3.0, whose
rerun equals the `dev` pipeline output (repA/dnaA hits; the fallbacks were not used). The Autocycler image already
runs 1.4.0 (Plassembler dependency). Image `dnaapler:1.4.0--pyhdfd78af_0` (sha256:75969182…).

### Kleborate

The newest release is the beta 3.3.0b1 (2026-09-27: new rmp, cgMLST/LIN codes, capsule typing for E. coli, Kleborate
AMR database 3.3.0, Kaptive 3.3.2, different aligner) — a beta that changes outputs and columns, not a candidate. The
3.2.4 `_1` rebuild (2026-08-08) moves Kaptive 3.1.0 → 3.2.2, minimap2 2.30 → 2.31, AMRFinderPlus 4.0.23 → 4.2.7
(not used by the `kpsc` preset) and Python 3.13 → 3.14; with `-p kpsc --trim_headers` both outputs
(`klebsiella_pneumo_complex_output.txt`, 117 columns; hAMRonization, 32 columns) were identical in every field on
SCL2333, SCL2535 and SCL2461. Image sha256:b6542e6e….

### AMRFinderPlus

4.2.7 (2026-01-26): "a bug fix release, it should make no difference in AMRFinderPlus results"; fixes a crash when an
HMM hit has no BLAST alignment (ncbi/amr#180). `amrfinder -h` unchanged; same database version (2025-12-03.1) and
organism list. All 10 result tables and 8 mutation tables byte-identical to 4.2.5. Image sha256:ea12f3f7….

### Mash

No release since 2.3 (2021); `hf85e966_11` (2026-07-22) is a rebuild. `mash sketch -k 32 -m 3 -r` on SCL2535 and 100.1
gave byte-identical sketches and stats; only the help text is re-wrapped. The pinned `he348c14_1` is still on quay.io.

### SCCmec, sistr_cmd

SCCmec 1.2.0 (2024-08-25) and sistr_cmd 1.1.3 `pyhdc42f0e_2` (2025-04-24) are the newest release and build.

### GTDB-Tk

The pinned nf-core module (git sha 5e433c46, patched for a batch file) uses 2.5.2. GTDB r232 (the section 7 database,
`VERSION_DATA=r232`) needs GTDB-Tk 2.7.0 or later (2.7.0 notes: "requires a new reference package"); 2.5.2 accepts
r220/r226 only, so `--run_gtdbtk` cannot run with the server's database today. Changes since 2.5.2 that affect outputs:
2.6.0 fixes skani 0.3.1 and pplacer 1.1.alpha19 and removes the 100-genome limit in dense genera; 2.7.0 uses a
pre-sketched skani database and replaces `--skip_ani_screen` with `--place_species` (the module passes neither) and warns
that divide-and-conquer placement needs ≥ 140 GB; 2.8.0 makes `classify_wf` take the closest reference with AF ≥
`--min_af` before checking its ANI radius (affects fuzzy species in condensed genera) and adds
`ani_screen_unassigned_hits.tsv.gz`. Phase 8a ran 2.8.0 (`gtdbtk:2.8.0--pyhdfd78af_0`, sha256:d886dde1…, `check_install`
OK) through a container override on the routing set with the module's current arguments: 1 h 20 min, peak RSS 176.7 GB
(the label requests 200 GB). The module's version `sed` still yields `2.8.0`. nf-core master now ships 2.7.2 (Seqera
image, versions as topics, a license header); `nf-core modules update` would bring that and need the batch-file patch
re-applied and a 2.8.0 container patch on top.

### PlasmidFinder

The nf-core module still pins 2.1.6. `2.1.6--py314hdfd78af_2` (2025-11-06) is a rebuild with BLAST 2.17.0 and KMA 1.6.7
(the pinned `py310_1` has BLAST 2.12.0, KMA 1.4.3) and the same database: 460 replicons, 2019-08 FASTA files, identical
MD5 of all `.fsa` files. Both builds gave identical `results_tab.tsv` on the 10 routing-set annotation FASTAs and
SCL10571 (0–7 replicons per routing sample, 21 in all; SCL10571's 8 equal the `dev` pipeline's). 3.0.3 (2026-10) is a rewrite:
the command is `plasmidfinder`, output is BeOne JSON (with `--legacy` for the old JSON), the database path comes from
`-p` or `CGE_PLASMIDFINDER_DB`, and no database is bundled; it would need a module rewrite (nf-core has none) and a
database you supply. Whether a newer PlasmidFinder database changes the 9a fixtures' calls cannot be answered from
what is on the server: `--plasmidfinder_db` exists for that, and a database you provide can be compared the same way.

## Reproducing

Paths are the benchmark's scratch directories; `$P` is the Phase 8a results directory.

```bash
# dev baseline (SCL10571, Flye + Berokka + Medaka 1.11.3)
nextflow run . -profile docker --input ss_ont.csv --seq_type ont --run_kraken2 false \
  --outdir dev_ont -work-dir dev_ont.work -c cap.config      # cap.config: resourceLimits 32 CPUs / 600 GB

# Medaka, same input as the dev MEDAKA task (draft copied into the writable dir: medaka writes a .fai beside it)
docker run --rm --cpus 12 --memory 72g -u $(id -u):$(id -g) -e HOME=/w -v $PWD/inputs:/in:ro -v $PWD/run:/w -w /w \
  quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0 \
  medaka_consensus -t 12 [--bacteria | -m r1041_e82_400bps_sup_v4.3.0] \
  -i /in/SCL10571.fastplong.fastq.gz -d /w/draft.fasta -o ./
dnadiff -p a_vs_d A/consensus.fasta D/consensus.fasta

# model selection as a pipeline task sees it (no HOME, non-root)
docker run --rm -u $(id -u):$(id -g) -v $PWD/heads:/h:ro quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0 \
  medaka tools resolve_model --auto_model consensus_bacteria /h/SCL10571.fastq

# RGI (per version): RGI_LOAD, then RGI_MAIN on a copy of the annotation FASTA
rgi load --card_json /mds_data/dbs/rgi/card.json --local
python -c "from app.Database import Database; Database(local_database=True).build_databases()"
rgi main --input_sequence <sample>.fasta --output_file <sample> --local -n 6 -t contig -a DIAMOND -d wgs

# AMRFinderPlus (per version), as AMRFINDERPLUS_RUN
amrfinder --plus -n $P/<s>/annotation/bakta/<s>.fna -p ...faa -g ...gff3 -a bakta \
  [--organism <org> --mutation_all <s>-mutations.tsv] --database /mds_data/dbs/amrfinder/data/latest --threads 8

# Dnaapler, Kleborate, Mash, PlasmidFinder: the module command lines
dnaapler all -i SCL10571_polished.fasta -o output -t 6
kleborate -a <s>.fasta -o . -p kpsc --trim_headers
mash sketch -k 32 -m 3 -o <s> -r <R1> <R2>
plasmidfinder.py -i <s>.fasta -o ./ -x
```
