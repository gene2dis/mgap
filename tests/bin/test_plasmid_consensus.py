"""Tests for bin/plasmid_consensus.py.

The tool outputs in data/plasmid/ are excerpts of real MGAP outputs for SCL2535 (Phase 8a run of the
Illumina routing set: geNomad, mob_recon) and PlasmidFinder 2.1.6 run on the same annotation FASTA
(Phase 9a). The FASTA is written by each test with short sequences under the real contig headers.
Run with: pytest tests/bin
"""

import csv

import pytest

import plasmid_consensus as pc

SCL2535_HEADERS = [
    ("contig_1 NODE_1_length_1682762_cov_30.909424", 120),
    ("contig_11 NODE_11_length_110042_cov_17.045271", 40),
    ("contig_16 NODE_16_length_43835_cov_25.974192", 30),
    ("contig_22 NODE_22_length_7041_cov_26.434336", 20),
    ("contig_24 NODE_24_length_5849_cov_144.248864", 20),
    ("contig_27 NODE_27_length_4497_cov_37.351259", 10),
]


def write_fasta(path, headers):
    with open(path, "w") as handle:
        for header, length in headers:
            # two sequence lines, so lengths are summed across lines
            half = length // 2
            handle.write(f">{header}\n{'A' * half}\n{'C' * (length - half)}\n")
    return path


def read_rows(path):
    with open(path, newline="") as handle:
        return {row["contig"]: row for row in csv.DictReader(handle, delimiter="\t")}


def run_consensus(tmp_path, data, *extra, headers=SCL2535_HEADERS, sources=("genomad", "mob", "plasmidfinder")):
    fasta = write_fasta(tmp_path / "SCL2535.fasta", headers)
    out = tmp_path / "consensus.tsv"
    argv = ["consensus", "--sample", "SCL2535", "--fasta", str(fasta), "--out", str(out), *extra]
    plasmid = data / "plasmid"
    if "genomad" in sources:
        argv += ["--genomad", str(plasmid / "SCL2535.genomad_plasmid.tsv")]
    if "mob" in sources:
        argv += ["--mob-contigs", str(plasmid / "SCL2535.contig_report.txt")]
        argv += ["--mob-plasmids", str(plasmid / "SCL2535.mobtyper_results.txt")]
    if "plasmidfinder" in sources:
        argv += ["--plasmidfinder", str(plasmid / "SCL2535.plasmidfinder.tsv")]
    pc.main(argv)
    return read_rows(out)


# --------------------------------------------------------------------------------------------------
# Consensus calls
# --------------------------------------------------------------------------------------------------


def test_consensus_calls_from_three_sources(tmp_path, data):
    rows = run_consensus(tmp_path, data)
    assert list(rows) == ["contig_1", "contig_11", "contig_16", "contig_22", "contig_24", "contig_27"]
    calls = {contig: (row["consensus"], row["support"], row["evidence"]) for contig, row in rows.items()}
    assert calls == {
        "contig_1": ("chromosome", "0", "MOB-suite: chromosome"),
        "contig_11": ("plasmid", "2", "MOB-suite, PlasmidFinder"),
        "contig_16": ("plasmid", "3", "geNomad, MOB-suite, PlasmidFinder"),
        "contig_22": ("candidate", "1", "MOB-suite only"),
        "contig_24": ("plasmid", "2", "geNomad, MOB-suite"),
        "contig_27": ("candidate", "1", "geNomad only"),
    }
    assert rows["contig_1"]["length"] == "120"
    assert rows["contig_1"]["description"] == "NODE_1_length_1682762_cov_30.909424"


def test_consensus_evidence_columns(tmp_path, data):
    rows = run_consensus(tmp_path, data)
    assert rows["contig_11"]["plasmidfinder_replicons"] == "IncFIB(pKPHS1)"
    assert rows["contig_11"]["plasmidfinder_identity"] == "97.32"
    assert rows["contig_11"]["mobsuite_cluster"] == "AC125"
    assert rows["contig_11"]["mobsuite_replicons"] == "IncFIB"
    # mobility of the mob_recon plasmid the contig belongs to
    assert (rows["contig_11"]["predicted_mobility"], rows["contig_11"]["mobility_source"]) == (
        "non-mobilizable",
        "mob_recon cluster",
    )
    assert rows["contig_16"]["predicted_mobility"] == "conjugative"
    assert rows["contig_16"]["mobsuite_relaxases"] == "MOBF"
    assert rows["contig_16"]["genomad_conjugation_genes"].startswith("F_traL;")
    # a novel mob_recon cluster is matched by its primary cluster id
    assert rows["contig_24"]["predicted_mobility"] == "mobilizable"
    # geNomad's direct terminal repeats make the contig circular
    assert (rows["contig_24"]["circular"], rows["contig_24"]["circular_source"]) == ("yes", "geNomad DTR")
    assert rows["contig_16"]["circular"] == "-"
    # mobility is only given for mob_recon plasmids
    assert rows["contig_27"]["predicted_mobility"] == "-"
    assert rows["contig_1"]["genomad_plasmid"] == "no"
    assert rows["contig_1"]["mob_typer"] == "no"


def test_large_contig_with_evidence_is_check(tmp_path, data):
    rows = run_consensus(tmp_path, data, "--large-contig", "40")
    assert rows["contig_11"]["consensus"] == "check"
    assert rows["contig_11"]["evidence"] == "MOB-suite, PlasmidFinder (contig of 40 bp)"
    # no evidence: the length does not matter
    assert rows["contig_1"]["consensus"] == "chromosome"
    assert rows["contig_16"]["consensus"] == "plasmid"


def test_tools_not_run(tmp_path, data):
    rows = run_consensus(tmp_path, data, sources=("genomad",))
    assert rows["contig_1"]["consensus"] == "no call"
    assert rows["contig_1"]["evidence"] == "no plasmid evidence"
    assert rows["contig_16"]["consensus"] == "candidate"
    for column in ("mobsuite_molecule_type", "predicted_mobility", "mob_typer", "plasmidfinder_replicons"):
        assert rows["contig_16"][column] == pc.NOT_RUN
    rows = run_consensus(tmp_path, data, sources=())
    assert {row["consensus"] for row in rows.values()} == {"no call"}
    assert rows["contig_1"]["evidence"] == "-"


def test_plasmidfinder_without_hits(tmp_path, data):
    fasta = write_fasta(tmp_path / "a.fasta", SCL2535_HEADERS[:2])
    out = tmp_path / "out.tsv"
    pc.main(
        [
            "consensus",
            "--sample",
            "A",
            "--fasta",
            str(fasta),
            "--plasmidfinder",
            str(data / "plasmid" / "empty.plasmidfinder.tsv"),
            "--out",
            str(out),
        ]
    )
    rows = read_rows(out)
    assert rows["contig_11"]["plasmidfinder_replicons"] == "-"
    assert rows["contig_11"]["consensus"] == "no call"


def test_mob_typer_on_circular_contigs(tmp_path, data):
    # MOB-typer --multi names each row by the whole FASTA header
    typer = tmp_path / "typer.txt"
    columns = ["sample_id", "num_contigs", "size", "rep_type(s)", "relaxase_type(s)", "mpf_type", "orit_type(s)"]
    columns += ["predicted_mobility", "mash_nearest_neighbor", "mash_neighbor_distance", "mash_neighbor_identification"]
    typer.write_text(
        "\t".join(columns)
        + "\n"
        + "\t".join(["contig_2 circular=true", "1", "40", "rep_cluster_2401", "-", "-", "MOB_unknown"])
        + "\tmobilizable\tCP022918\t0.09\tKlebsiella pneumoniae\n"
        + "\t".join(["contig_3 circular=true", "1", "30", "-", "-", "-", "-"])
        + "\tnon-mobilizable\t-\t-\t-\n"
    )
    headers = [("contig_1 circular=true", 120), ("contig_2 circular=true", 40), ("contig_3 circular=true", 30)]
    fasta = write_fasta(tmp_path / "ont.fasta", headers)
    out = tmp_path / "out.tsv"
    pc.main(["consensus", "--sample", "S", "--fasta", str(fasta), "--mob-typer", str(typer), "--out", str(out)])
    rows = read_rows(out)
    # a replicon or relaxase from MOB-typer is the MOB-suite vote
    assert (rows["contig_2"]["consensus"], rows["contig_2"]["evidence"]) == ("candidate", "MOB-suite only")
    assert (rows["contig_2"]["predicted_mobility"], rows["contig_2"]["mobility_source"]) == ("mobilizable", "MOB-typer")
    assert rows["contig_2"]["mash_neighbor_identification"] == "Klebsiella pneumoniae"
    assert rows["contig_2"]["mobsuite_molecule_type"] == pc.NOT_RUN
    assert (rows["contig_2"]["circular"], rows["contig_2"]["circular_source"]) == ("yes", "header")
    # typed, but neither replicon nor relaxase: no vote
    assert rows["contig_3"]["consensus"] == "no call"
    assert rows["contig_3"]["mob_typer"] == "yes"
    assert rows["contig_1"]["mob_typer"] == "no"


# --------------------------------------------------------------------------------------------------
# Circularity and the MOB-typer input
# --------------------------------------------------------------------------------------------------

FLYE_INFO = (
    "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
    "contig_1\t5333000\t96\tY\tN\t1\t*\t1\n"
    "contig_7\t239000\t120\tY\tN\t1\t*\t7\n"
    "contig_9\t4100\t80\tN\tN\t1\t*\t9\n"
)


def test_circular_from_flye_info(tmp_path):
    # after renaming, the Flye name is the first word of the description; Flye's order is not kept
    headers = [("contig_1 contig_1", 50), ("contig_2 contig_7 rotated=True", 30), ("contig_3 contig_9", 20)]
    fasta = write_fasta(tmp_path / "a.fasta", headers)
    info = tmp_path / "assembly_info.txt"
    info.write_text(FLYE_INFO)
    out = tmp_path / "circular.fasta"
    pc.main(["circular", "--fasta", str(fasta), "--flye-info", str(info), "--out", str(out)])
    records = pc.read_fasta(out)
    assert [(record[0], record[1], record[2]) for record in records] == [
        ("contig_1", "contig_1", 50),
        ("contig_2", "contig_7 rotated=True", 30),
    ]
    shape = pc.topology(pc.read_fasta(fasta), info)
    assert shape["contig_3"] == ("no", "Flye")


@pytest.mark.parametrize(
    "description, expected",
    [
        ("length=5333000 circular=true", "yes"),
        ("length=5333000 circular=True depth=1.00x", "yes"),
        ("[gcode=11] [topology=circular]", "yes"),
        ("topology=circular", "yes"),
        ("circular=false", "no"),
        ("[topology=linear]", "no"),
        ("NODE_1_length_1682762_cov_30.909424", "-"),
        ("notcircular=true", "-"),
    ],
)
def test_circular_header_tags(tmp_path, description, expected):
    fasta = write_fasta(tmp_path / "a.fasta", [(f"contig_1 {description}", 10)])
    assert pc.topology(pc.read_fasta(fasta))["contig_1"][0] == expected


def test_circular_contig_mode_fallback(tmp_path):
    out = tmp_path / "circular.fasta"
    # no topology tag at all: every contig with --all-if-untagged, none without it
    untagged = write_fasta(tmp_path / "untagged.fasta", [("contig_1 chr", 30), ("contig_2 p1", 10)])
    pc.main(["circular", "--fasta", str(untagged), "--out", str(out)])
    assert not out.exists()
    pc.main(["circular", "--fasta", str(untagged), "--all-if-untagged", "--out", str(out)])
    assert [record[0] for record in pc.read_fasta(out)] == ["contig_1", "contig_2"]
    out.unlink()
    # tags present (here only linear ones): the tags decide, even when none is circular
    tagged = write_fasta(tmp_path / "tagged.fasta", [("contig_1 circular=false", 30), ("contig_2 circular=false", 10)])
    pc.main(["circular", "--fasta", str(tagged), "--all-if-untagged", "--out", str(out)])
    assert not out.exists()
