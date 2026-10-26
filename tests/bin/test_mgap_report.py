"""Parser and manifest tests for bin/mgap_report.py.

The files in data/ are excerpts of real MGAP outputs (Phase 8a run of the Illumina routing set,
2026-10-25), with work-directory paths removed. Run with: pytest tests/bin
"""

import json

import pytest

import mgap_report as report

# --------------------------------------------------------------------------------------------------
# Manifest
# --------------------------------------------------------------------------------------------------


def test_manifest_tsv_relative_paths_and_run_level_rows(tmp_path, data):
    (tmp_path / "files").mkdir()
    (tmp_path / "files" / "a.json").write_text("{}")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample\ttool\tfile\n"
        "S1\tfastp\tfiles/a.json\n"
        f"-\tcitations\t{data / 'CITATIONS.md'}\n"
        "S1\tunknown_tool\tfiles/a.json\n"
        "S2\tquast\tfiles/missing.tsv\n"
    )
    parsed = report.read_manifest(manifest)
    assert parsed["fastp"]["S1"] == [tmp_path / "files" / "a.json"]
    assert parsed["citations"]["-"] == [data / "CITATIONS.md"]
    assert "unknown_tool" not in parsed
    assert "S2" not in parsed.get("quast", {})


def test_manifest_json_list_and_object(tmp_path, data):
    rows = [{"sample": "SCL2535", "tool": "mash", "file": str(data / "SCL2535.genome_size")}]
    for content in (rows, {"files": rows}):
        manifest = tmp_path / "manifest.json"
        manifest.write_text(json.dumps(content))
        assert report.read_manifest(manifest)["mash"]["SCL2535"] == [data / "SCL2535.genome_size"]


def test_manifest_empty_sample_is_run_level(tmp_path, data):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"sample\ttool\tfile\n\tsample_status\t{data / 'sample_status.tsv'}\n")
    assert report.read_manifest(manifest)["sample_status"]["-"] == [data / "sample_status.tsv"]


# --------------------------------------------------------------------------------------------------
# Read QC
# --------------------------------------------------------------------------------------------------


def test_parse_fastp(data):
    result = report.parse_fastp(data / "SCL2535.fastp.json")
    assert result["reads_before"] == 2270906
    assert result["bases_before"] == 551085099
    assert result["bases_after"] == 472142971
    assert result["too_short_reads"] == 119844
    assert result["duplication_rate"] == pytest.approx(0.00348231)
    assert result["mean_length_before"] == pytest.approx(242.7, abs=0.05)


def test_parse_fastplong(data):
    result = report.parse_fastp(data / "ont.fastplong.json")
    assert result["reads_before"] == 1000
    assert result["bases_after"] == 33311690
    assert result["duplication_rate"] is None


def test_parse_mash_genome_size(data):
    assert report.parse_mash_genome_size(data / "SCL2535.genome_size") == 5012860


def test_parse_kraken2_species_rows_and_unclassified(data):
    species, unclassified = report.parse_kraken2(data / "SCL2333.kraken2.report.txt")
    assert unclassified == pytest.approx(0.24)
    assert list(species["name"]) == [
        "Klebsiella pneumoniae",
        "Raoultella ornithinolytica",
        "Klebsiella quasipneumoniae",
    ]
    assert species.loc[0, "percent"] == pytest.approx(27.88)
    assert species.loc[0, "reads"] == 329949


def test_parse_bracken_percent_of_classified_reads(data):
    table = report.parse_bracken(data / "SCL2333.bracken.tsv")
    assert table.loc[0, "name"] == "Klebsiella pneumoniae"
    assert table.loc[0, "percent"] == pytest.approx(96.143)
    assert table.loc[0, "reads"] == 1134253


# --------------------------------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------------------------------


def test_parse_quast(data):
    result = report.parse_quast(data / "511.1.quast.tsv")
    assert result["contigs_all"] == 397
    assert result["contigs"] == 48
    assert result["total_length"] == 4830321
    assert result["n50"] == 295713
    assert result["gc"] == pytest.approx(51.92)


def test_parse_filter_stats(data):
    result = report.parse_filter_stats(data / "511.1.filter_stats.tsv")
    assert result == {
        "assembly_contigs": 397,
        "assembly_length": 4907697,
        "assembly_n50": 295713,
        "annotation_contigs": 38,
        "annotation_length": 4823885,
        "annotation_n50": 295713,
    }


# --------------------------------------------------------------------------------------------------
# Typing
# --------------------------------------------------------------------------------------------------


def test_parse_mlst(data):
    result = report.parse_mlst(data / "100.1.mlst.tsv")
    assert result["scheme"] == "senterica_achtman_2"
    assert result["st"] == "32"
    assert result["alleles"].startswith("aroC(17) dnaN(18)")


def test_parse_mlst_without_scheme(tmp_path):
    path = tmp_path / "x.tsv"
    path.write_text("x.fasta\t-\t-\n")
    assert report.parse_mlst(path) == {"scheme": None, "st": None, "alleles": ""}


def test_parse_sistr(data):
    result = report.parse_sistr(data / "511.1.sistr.tab")
    assert result["serovar"] == "Agona"
    assert result["serogroup"] == "B"
    assert result["qc_status"] == "PASS"


def test_parse_kleborate_3_output(data):
    result = report.parse_kleborate(data / "SCL2333.kleborate.txt")
    assert result["species"] == "Klebsiella pneumoniae"
    assert result["ST"] == "ST25"
    assert result["K_locus"] == "KL2"
    assert result["Bla_ESBL_acquired"] == "CTX-M-15"
    assert result["virulence_score"] == "1"


def test_parse_kleborate_strips_module_prefixes(tmp_path):
    path = tmp_path / "klebsiella_pneumo_complex_output.txt"
    path.write_text(
        "strain\tenterobacterales__species__species\tklebsiella_pneumo_complex__mlst__ST\t"
        "klebsiella__abst__ST\n"
        "S1\tKlebsiella pneumoniae\tST258\t0\n"
    )
    result = report.parse_kleborate(path)
    assert result["species"] == "Klebsiella pneumoniae"
    assert result["ST"] == "ST258"


def test_parse_sccmec(data):
    result = report.parse_sccmec(data / "SCL1559.sccmec.tsv")
    assert (result["type"], result["subtype"], result["mecA"]) == ("I", "Ia", "+")


# --------------------------------------------------------------------------------------------------
# Genome quality, AMR and mobile elements
# --------------------------------------------------------------------------------------------------


def test_parse_checkm2(data):
    result = report.parse_checkm2(data / "511.1.checkm2.tsv")
    assert result["completeness"] == 100.0
    assert result["contamination"] == pytest.approx(0.32)
    assert result["genome_size"] == 4823885
    assert result["model"].startswith("Neural Network")


def test_parse_amrfinderplus_4_columns(data):
    table = report.parse_amrfinderplus(data / "511.1.amrfinder.tsv")
    assert {"gene", "type", "subtype", "class", "identity", "coverage", "contig"} <= set(table.columns)
    amr = table[table["type"] == "AMR"]
    assert list(amr["gene"]) == ["mdsB", "mdsA", "fosA7.2"]
    assert amr["identity"].iloc[0] == pytest.approx(99.72)


def test_parse_amrfinderplus_3_columns(tmp_path):
    path = tmp_path / "old.tsv"
    path.write_text(
        "Gene symbol\tSequence name\tElement type\tElement subtype\tClass\tSubclass\t"
        "% Coverage of reference sequence\t% Identity to reference sequence\n"
        "blaKPC-2\tcarbapenem-hydrolyzing class A beta-lactamase KPC-2\tAMR\tAMR\tBETA-LACTAM\tCARBAPENEM\t100.00\t100.00\n"
    )
    table = report.parse_amrfinderplus(path)
    assert table.loc[0, "gene"] == "blaKPC-2"
    assert table.loc[0, "subclass"] == "CARBAPENEM"
    assert table.loc[0, "identity"] == 100.0


def test_parse_rgi(data):
    table = report.parse_rgi(data / "SCL1559.rgi.txt")
    assert list(table["gene"]) == ["norC", "arlR", "arlS", "ANT(9)-Ia", "ErmA"]
    assert list(table["cut_off"][:2]) == ["Strict", "Perfect"]
    assert table.loc[0, "identity"] == pytest.approx(70.64)
    # the ORF number is removed so contigs match the other tools
    assert list(table["contig"]) == ["contig_1", "contig_1", "contig_1", "contig_2", "contig_2"]


def test_parse_genomad(data):
    plasmids = report.parse_genomad(data / "100.1.genomad_plasmid.tsv")
    viruses = report.parse_genomad(data / "100.1.genomad_virus.tsv")
    assert len(plasmids) == 5
    assert plasmids.loc[0, "seq_name"] == "contig_26"
    assert len(viruses) == 2


def test_parse_mobsuite(data):
    contigs = report.parse_mobsuite_contigs(data / "SCL1559.mob_contigs.txt")
    plasmids = report.parse_mobsuite_mobtyper(data / "SCL1559.mobtyper.txt")
    assert (contigs["molecule_type"] == "plasmid").sum() == 7
    assert list(plasmids["primary_cluster_id"]) == ["AA849", "AA769", "AE018", "AE387"]
    assert plasmids.loc[0, "predicted_mobility"] == "mobilizable"


# --------------------------------------------------------------------------------------------------
# Methods
# --------------------------------------------------------------------------------------------------


def test_parse_citations_headings_and_list_items(data):
    entries = report.parse_citations(data / "CITATIONS.md")
    assert [entry["name"] for entry in entries] == ["nf-core", "Nextflow", "Autocycler", "GTDB-Tk", "LRGE"]
    assert entries[3]["url"] == "https://pubmed.ncbi.nlm.nih.gov/31730192/"
    assert entries[3]["text"].startswith("Chaumeil PA")
    assert entries[3]["text"].endswith("PMC7703759.")


def test_parse_software_versions_workflow_first(data):
    rows = report.parse_software_versions(data / "software_versions.yml")
    assert [(row["tool"], row["version"]) for row in rows[:2]] == [("gene2dis/mgap", "2.0.0"), ("Nextflow", "26.04.4")]
    rgi = next(row for row in rows if row["tool"] == "rgi")
    # one row per tool and version, listing every process that used it
    assert rgi["processes"] == ["GENE2DIS_MGAP:MGAP:RGI_LOAD", "GENE2DIS_MGAP:MGAP:RGI_MAIN"]


# --------------------------------------------------------------------------------------------------
# Run-level tables
# --------------------------------------------------------------------------------------------------


def test_parse_gtdbtk(data):
    result = report.parse_gtdbtk(data / "gtdbtk.bac120.summary.tsv")
    assert result["SCL2461"]["species"] == "Klebsiella pneumoniae"
    assert result["SCL2461"]["ani"] == pytest.approx(99.11)
    assert result["SCL2461"]["method"] == "ani_screen"
    # placed by the tree only: genus known, no species and no ANI
    assert result["CCG_02"]["species"] is None
    assert result["CCG_02"]["genus"] == "Micromonospora"
    assert result["CCG_02"]["ani"] is None


def test_parse_sample_status(data):
    table = report.parse_sample_status(data / "sample_status.tsv")
    assert list(table["status"]) == ["pass", "pass", "skipped"]
    low = table[table["sample"] == "SCL2535_low"].iloc[0]
    assert low["stage"] == "read_qc"
    assert low["reason"].startswith("16648045 bases after fastp")


# --------------------------------------------------------------------------------------------------
# Whole report from the fixtures
# --------------------------------------------------------------------------------------------------

FIXTURE_MANIFEST = [
    ("SCL2535", "fastp", "SCL2535.fastp.json"),
    ("SCL2535", "mash", "SCL2535.genome_size"),
    ("SCL2333", "kraken2", "SCL2333.kraken2.report.txt"),
    ("SCL2333", "bracken", "SCL2333.bracken.tsv"),
    ("SCL2333", "kleborate", "SCL2333.kleborate.txt"),
    ("511.1", "quast", "511.1.quast.tsv"),
    ("511.1", "filter_stats", "511.1.filter_stats.tsv"),
    ("511.1", "checkm2", "511.1.checkm2.tsv"),
    ("511.1", "sistr", "511.1.sistr.tab"),
    ("511.1", "amrfinderplus", "511.1.amrfinder.tsv"),
    ("100.1", "mlst", "100.1.mlst.tsv"),
    ("100.1", "genomad_plasmid", "100.1.genomad_plasmid.tsv"),
    ("100.1", "genomad_virus", "100.1.genomad_virus.tsv"),
    ("SCL1559", "rgi", "SCL1559.rgi.txt"),
    ("SCL1559", "sccmec", "SCL1559.sccmec.tsv"),
    ("SCL1559", "mobsuite_contigs", "SCL1559.mob_contigs.txt"),
    ("SCL1559", "mobsuite_mobtyper", "SCL1559.mobtyper.txt"),
    ("-", "gtdbtk", "gtdbtk.bac120.summary.tsv"),
    ("-", "sample_status", "sample_status.tsv"),
    ("-", "software_versions", "software_versions.yml"),
    ("-", "citations", "CITATIONS.md"),
]


def test_report_from_fixtures(tmp_path, data):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample\ttool\tfile\n" + "".join(f"{s}\t{tool}\t{data / name}\n" for s, tool, name in FIXTURE_MANIFEST)
    )
    html, tsv = tmp_path / "report.html", tmp_path / "summary.tsv"
    report.main(["--manifest", str(manifest), "--html", str(html), "--tsv", str(tsv)])
    page = html.read_text()
    for section in ("overview", "read-qc", "assembly", "quality", "typing", "amr", "mobile", "samples", "methods"):
        assert f'<section id="{section}"' in page
    # offline: Plotly is inlined, nothing is loaded from the network
    assert "<script src=" not in page
    summary = tsv.read_text().splitlines()
    assert summary[0].startswith("sample\t")
    samples = {line.split("\t")[0] for line in summary[1:]}
    assert {"SCL2535", "SCL2333", "511.1", "100.1", "SCL1559", "SCL2535_low", "SCL2461"} <= samples


def test_report_sections_need_their_inputs(tmp_path, data):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"sample\ttool\tfile\n511.1\tquast\t{data / '511.1.quast.tsv'}\n")
    html = tmp_path / "report.html"
    report.main(["--manifest", str(manifest), "--html", str(html), "--tsv", str(tmp_path / "summary.tsv")])
    page = html.read_text()
    assert '<section id="assembly"' in page
    for section in ("read-qc", "quality", "typing", "amr", "mobile"):
        assert f'<section id="{section}"' not in page
