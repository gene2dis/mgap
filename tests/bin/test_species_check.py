"""Tests for bin/species_check.py.

data/species/: SCL2461's mlst row (ecoli_achtman_4 ST14464, Phase 8a run) and its skani search against the
Phase 12 reference set (K. pneumoniae at 99.11 % ANI), and CCG_01's empty search (Micromonospora, not in
the set). The SCL2333 Bracken table and Kraken2 report are shared with the report tests. The reference
list is the pipeline's own assets/species_refs.tsv.
Run with: pytest tests/bin
"""

import csv
import json
import re
from pathlib import Path

import pytest

import species_check as sc

REPO = Path(__file__).resolve().parents[2]
REFS = REPO / "assets" / "species_refs.tsv"
# amrfinder --list_organisms, AMRFinderPlus 4.2.7 with database 2025-12-03.1 (Phase 3 and 10 results)
AMRFINDER_ORGANISMS = {
    "Acinetobacter_baumannii", "Bordetella_pertussis", "Burkholderia_cepacia", "Burkholderia_mallei",
    "Burkholderia_pseudomallei", "Campylobacter", "Citrobacter_freundii", "Clostridioides_difficile",
    "Corynebacterium_diphtheriae", "Enterobacter_asburiae", "Enterobacter_cloacae", "Enterococcus_faecalis",
    "Enterococcus_faecium", "Escherichia", "Haemophilus_influenzae", "Helicobacter_pylori", "Klebsiella_oxytoca",
    "Klebsiella_pneumoniae", "Neisseria_gonorrhoeae", "Neisseria_meningitidis", "Pseudomonas_aeruginosa",
    "Salmonella", "Serratia_marcescens", "Staphylococcus_aureus", "Staphylococcus_pseudintermedius",
    "Streptococcus_agalactiae", "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae",
    "Vibrio_parahaemolyticus", "Vibrio_vulnificus",
}  # fmt: skip


def taxa_names():
    """getTaxaNames() from workflows/mgap.nf, as the pipeline writes it to mlst_organisms.json."""
    text = (REPO / "workflows" / "mgap.nf").read_text()
    body = re.search(r"def getTaxaNames\(\)\s*\{(.*?)\n\}", text, re.S).group(1)
    return dict(re.findall(r'"([^"]+)":\s*"([^"]+)"', body))


@pytest.fixture
def run(tmp_path, data):
    """Runs species_check.py; files are given as {name: content} or paths; returns the output row."""
    organisms = tmp_path / "mlst_organisms.json"
    organisms.write_text(json.dumps(taxa_names()))

    def _run(mlst, skani=None, taxa=None, taxa_format="bracken", expected=None, extra=()):
        def as_file(name, content):
            if content is None or isinstance(content, Path):
                return content
            path = tmp_path / name
            path.write_text(content)
            return path

        args = ["--sample", "S1", "--seq-type", "illumina", "--refs", str(REFS)]
        args += ["--mlst-organisms", str(organisms), "--out", str(tmp_path / "out.tsv")]
        args += ["--mlst", str(as_file("mlst.tsv", mlst))]
        if skani is not None:
            args += ["--skani", str(as_file("skani.tsv", skani))]
        if taxa is not None:
            args += ["--taxa", str(as_file("taxa.tsv", taxa)), "--taxa-format", taxa_format]
        if expected:
            args += ["--expected", expected]
        sc.main(args + list(extra))
        with open(tmp_path / "out.tsv") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 1
        return rows[0]

    return _run


SKANI_HEADER = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"


def skani(*hits):
    """skani rows from (accession, ANI, AF query)."""
    return SKANI_HEADER + "".join(
        f"genomes/{acc}.fna\tq.fasta\t{ani}\t80.0\t{af}\tref\tquery\n" for acc, ani, af in hits
    )


def mlst(scheme, st="1"):
    return f"S1.fasta\t{scheme}\t{st}\tadk(1)\n"


KPN, KVA, ECO, SAU, SAR = "GCF_006364295.1", "GCF_020525545.1", "GCF_000005845.2", "GCF_026167765.1", "GCF_000236925.1"
BPE, BBR, BMA, BPS = "GCF_004008975.1", "GCF_900636925.1", "GCF_033956065.1", "GCF_030297255.1"


def test_scl2461_override(run, data):
    row = run(data / "species" / "SCL2461.mlst.tsv", data / "species" / "SCL2461.skani.tsv")
    assert row["mlst_scheme"] == "ecoli_achtman_4" and row["mlst_st"] == "14464"
    assert row["mlst_organism"] == "Escherichia"
    assert row["skani_species"] == "Klebsiella pneumoniae" and row["skani_ani"] == "99.11"
    assert row["skani_af_query"] == "86.55"
    assert row["decision"] == "override"
    assert row["organism"] == "Klebsiella_pneumoniae"
    assert row["rerun_scheme"] == "klebsiella"
    assert row["flagged"] == "yes"
    assert "ecoli_achtman_4 does not type" in row["reason"]


def test_agree(run):
    row = run(mlst("klebsiella", "45"), skani((KPN, 99.08, 87.55), (KVA, 94.62, 81.71)))
    assert (row["decision"], row["organism"], row["flagged"], row["rerun_scheme"]) == (
        "agree",
        "Klebsiella_pneumoniae",
        "no",
        "",
    )


def test_relative_of_the_same_organism_agrees(run):
    """K. variicola is typed by the klebsiella scheme and AMRFinderPlus uses Klebsiella_pneumoniae for it."""
    row = run(mlst("klebsiella"), skani((KVA, 99.2, 88.0), (KPN, 94.6, 80.0)))
    assert (row["decision"], row["species"], row["organism"]) == (
        "agree",
        "Klebsiella variicola",
        "Klebsiella_pneumoniae",
    )


def test_refine_to_no_organism_is_flagged(run):
    """S. argenteus is typed by the saureus scheme but is not an AMRFinderPlus S. aureus."""
    row = run(mlst("saureus", "2250"), skani((SAR, 99.5, 90.0), (SAU, 89.4, 60.0)))
    assert (row["decision"], row["organism"], row["flagged"]) == ("refine", "", "yes")


def test_ambiguous_keeps_the_mlst_organism_when_it_is_a_candidate(run):
    row = run(mlst("bpseudomallei"), skani((BMA, 99.9, 95.0), (BPS, 99.6, 90.0)))
    assert row["decision"] == "ambiguous"
    assert row["organism"] == "Burkholderia_pseudomallei"
    assert row["species"] == "Burkholderia mallei / Burkholderia pseudomallei"
    assert row["skani_next"] == "Burkholderia pseudomallei (99.60)"
    assert row["flagged"] == "yes"
    # outside the margin the top hit decides
    clear = run(mlst("bpseudomallei"), skani((BMA, 99.9, 95.0), (BPS, 99.6, 90.0)), extra=["--ani-margin", "0.2"])
    assert (clear["decision"], clear["organism"]) == ("refine", "Burkholderia_mallei")


def test_reference_when_mlst_finds_no_scheme(run):
    row = run(mlst("-", "-"), skani((KPN, 99.0, 88.0)))
    assert (row["decision"], row["organism"], row["rerun_scheme"]) == (
        "reference",
        "Klebsiella_pneumoniae",
        "klebsiella",
    )
    assert row["mlst_scheme"] == "-"


def test_no_reference_hit_keeps_mlst_and_flags(run, data):
    row = run(mlst("saureus"), data / "species" / "CCG_01.skani.tsv")
    assert (row["decision"], row["organism"], row["flagged"]) == ("mlst", "Staphylococcus_aureus", "yes")
    assert row["reason"] == "no reference hit (no reference within the skani screen)"


def test_hits_below_the_thresholds_are_not_a_call(run):
    low_af = run(mlst("klebsiella"), skani((ECO, 99.0, 30.0)))
    assert (low_af["decision"], low_af["organism"], low_af["skani_species"]) == (
        "mlst",
        "Klebsiella_pneumoniae",
        "Escherichia coli",
    )
    assert "best: Escherichia coli, ANI 99.00, AF 30.00 %" in low_af["reason"]
    low_ani = run(mlst("klebsiella"), skani((KVA, 94.6, 80.0)))
    assert low_ani["decision"] == "mlst"


def test_without_reference_set(run):
    row = run(mlst("ecoli_achtman_4", "14464"))
    assert (row["decision"], row["organism"], row["flagged"]) == ("mlst", "Escherichia", "no")
    assert row["reason"] == "no reference set (--species_refs)"


def test_bracken_genus_check_without_reference_call(run, data):
    agrees = run(mlst("klebsiella"), taxa=data / "SCL2333.bracken.tsv")
    assert (agrees["reads_species"], agrees["reads_percent"], agrees["reads_source"]) == (
        "Klebsiella pneumoniae",
        "96.1",
        "Bracken",
    )
    assert agrees["flagged"] == "no"
    assert "Bracken Klebsiella pneumoniae agrees with the MLST genus" in agrees["reason"]
    # SCL2461 without --species_refs: the reads disagree with MLST, which still routes it (flagged)
    differs = run(mlst("ecoli_achtman_4", "14464"), taxa=data / "SCL2333.bracken.tsv")
    assert (differs["decision"], differs["organism"], differs["flagged"]) == ("mlst", "Escherichia", "yes")
    assert "Bracken top species Klebsiella pneumoniae (96.1 %) is another genus" in differs["reason"]


def test_kraken2_report_on_ont(run, data):
    row = run(mlst("klebsiella"), taxa=data / "SCL2333.kraken2.report.txt", taxa_format="kraken2")
    assert (row["reads_species"], row["reads_percent"], row["reads_source"]) == (
        "Klebsiella pneumoniae",
        "27.9",
        "Kraken2",
    )


def test_shigella_reads_count_as_escherichia(run):
    taxa = "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n"
    taxa += "Shigella flexneri\t623\tS\t10\t10\t20\t0.6\nEscherichia coli\t562\tS\t5\t5\t10\t0.4\n"
    row = run(mlst("ecoli_achtman_4"), taxa=taxa)
    assert row["flagged"] == "no"


def test_samplesheet_species(run):
    # used when there is no reference call and no MLST organism
    routed = run(mlst("-", "-"), expected="Salmonella enterica")
    assert (routed["decision"], routed["organism"], routed["species"]) == (
        "samplesheet",
        "Salmonella",
        "Salmonella enterica",
    )
    # recorded and compared otherwise
    differs = run(mlst("klebsiella"), skani((KPN, 99.0, 88.0)), expected="Escherichia coli")
    assert (differs["decision"], differs["organism"], differs["flagged"]) == ("agree", "Klebsiella_pneumoniae", "yes")
    assert "samplesheet species Escherichia coli differs" in differs["reason"]
    genus_only = run(mlst("klebsiella"), skani((KPN, 99.0, 88.0)), expected="Klebsiella")
    assert genus_only["flagged"] == "no"
    against_mlst = run(mlst("saureus"), expected="Klebsiella pneumoniae")
    assert against_mlst["flagged"] == "yes"
    assert "differs from the MLST genus" in against_mlst["reason"]
    unknown = run(mlst("saureus"), expected="Staphylococcus hominis")
    assert "not in the reference list" in unknown["reason"]


def test_empty_mlst_output(run):
    row = run("")
    assert (row["mlst_scheme"], row["mlst_st"], row["decision"], row["organism"]) == ("-", "-", "mlst", "")
