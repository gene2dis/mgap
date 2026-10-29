#!/usr/bin/env python3
"""MGAP run report: one self-contained HTML page, one per-sample TSV and an optional XLSX workbook (the
summary, then one sheet per tool) from a run's tool outputs.

The report never looks for files on its own. Every input is listed in a manifest, a TSV with the
columns ``sample``, ``tool`` and ``file`` (or a JSON list of objects with the same keys); in the
pipeline the files are staged into the task directory and the manifest names them. Run-level files
use ``-`` as the sample. Relative paths are read relative to the manifest's directory.

Tools (the ``tool`` column):

per sample
    fastp              fastp JSON (Illumina reads)
    fastplong          fastplong JSON (ONT reads)
    mash               Mash genome size estimate in bases (``*.genome_size``, Illumina)
    kraken2            Kraken2 report (``*.report.txt``)
    bracken            Bracken species table (``*.tsv``)
    quast              QUAST ``<sample>.tsv`` (two columns: metric, value)
    filter_stats       SeqKit stats of the assembly before and after the length filter (``*.filter_stats.tsv``)
    checkm2            CheckM2 quality report (``<sample>_checkm2_report.tsv``)
    mlst               mlst TSV (no header)
    amrfinderplus      AMRFinderPlus TSV
    rgi                RGI ``*.txt``
    genomad_plasmid    geNomad ``*_plasmid_summary.tsv``
    genomad_virus      geNomad ``*_virus_summary.tsv``
    mobsuite_contigs   MOB-suite ``contig_report.txt``
    mobsuite_mobtyper  MOB-suite ``mobtyper_results.txt``
    plasmid_consensus  ``<sample>.plasmid_consensus.tsv`` (bin/plasmid_consensus.py: one row per contig)
    isescan            ISEScan ``<sample>.isescan.tsv`` (one row per IS element; header only when none was found)
    kleborate          Kleborate 3 ``klebsiella_pneumo_complex_output.txt``
    sistr              SISTR ``<sample>.tab``
    sccmec             sccmec ``<sample>.tsv``

run level (sample ``-``)
    gtdbtk             GTDB-Tk ``*.summary.tsv`` (bac120 and ar53; genome names are sample IDs)
    sample_status      ``pipeline_info/sample_status.tsv``
    software_versions  ``pipeline_info/software_versions.yml``
    citations          the pipeline's ``CITATIONS.md``

An unknown tool, a missing file or a file that cannot be parsed gives a warning on stderr and the
report goes on without it; a section appears only when at least one of its inputs was read.

With more than LARGE_RUN (100) samples, the one-bar-per-sample charts become distributions (a box with
one point per sample) on which only the outliers are named; sample cards start collapsed in every run.

``--run-info`` takes an optional JSON object with what no output file holds: ``pipeline_version``,
``nextflow_version``, ``run_name``, ``date``, ``seq_type``, ``command``, ``parameters`` (name to value),
``databases`` (parameter name to path), ``database_versions`` (parameter name to version, read from
the database's own version file) and ``mlst_organisms`` (MLST scheme to the AMRFinderPlus organism,
the pipeline's ``getTaxaNames()``). Database versions that software_versions.yml records
(``amrfinderplus-database``, ``gtdb_db``) are used when run-info has none.
"""

import argparse
import base64
import csv
import html
import json
import math
import re
import sys
from collections import OrderedDict, defaultdict
from datetime import date
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs
from plotly.subplots import make_subplots

try:
    import yaml
except ImportError:  # software_versions.yml is then skipped with a warning
    yaml = None

SAMPLE_TOOLS = (
    "fastp",
    "fastplong",
    "mash",
    "kraken2",
    "bracken",
    "quast",
    "filter_stats",
    "checkm2",
    "mlst",
    "amrfinderplus",
    "rgi",
    "genomad_plasmid",
    "genomad_virus",
    "mobsuite_contigs",
    "mobsuite_mobtyper",
    "plasmid_consensus",
    "isescan",
    "kleborate",
    "sistr",
    "sccmec",
)
RUN_TOOLS = ("gtdbtk", "sample_status", "software_versions", "citations")
RUN_SAMPLE = "-"

# CheckM2 lines drawn on the scatter: the MIMAG high-quality draft thresholds
CHECKM2_MIN_COMPLETENESS = 90.0
CHECKM2_MAX_CONTAMINATION = 5.0
# A second species above this share of the reads is flagged on the sample card
SECOND_SPECIES_WARN = 5.0
# Taxa shown by name in the Kraken2/Bracken chart; the rest are summed as "Other"
TOP_TAXA = 7
HEATMAP_COLUMN = 34  # px per sample in the resistance heatmaps
HEATMAP_LABEL = 50  # characters of a heatmap row label; the hover text has the full name
# Above this many samples, one-bar-per-sample charts become distributions with only the outliers named
LARGE_RUN = 100
OUTLIER_LABELS = 12  # outliers named on a distribution panel; hover shows every point
TABLE_SCROLL_ROWS = 40  # longer tables scroll inside a box of TABLE_SCROLL_HEIGHT px
TABLE_SCROLL_HEIGHT = 560
# software_versions.yml entries that are database versions -> the database parameter they belong to
DATABASE_VERSION_ENTRIES = {"amrfinderplus-database": "amrfinder_db", "gtdb_db": "gtdbtk_db"}
# CheckM2 scatter: genomes outside the limits get their own colour and marker shape
CHECKM2_GROUPS = (
    ("low_completeness", "Low completeness", "#eda100", "triangle-down"),
    ("high_contamination", "High contamination", "#e34948", "diamond"),
    ("both", "Low completeness and high contamination", "#8f1d1b", "x"),
)

# Approximate total length (Mb) of complete genomes of each AMRFinderPlus organism, as a range;
# used only to flag assemblies far outside it. Keyed like the pipeline's getTaxaNames() values.
EXPECTED_GENOME_MB = {
    "Acinetobacter_baumannii": (3.7, 4.3),
    "Bordetella_pertussis": (3.9, 4.2),
    "Burkholderia_cepacia": (7.5, 9.0),
    "Burkholderia_pseudomallei": (7.0, 7.5),
    "Campylobacter": (1.6, 1.9),
    "Citrobacter_freundii": (4.9, 5.5),
    "Clostridioides_difficile": (3.9, 4.4),
    "Corynebacterium_diphtheriae": (2.3, 2.6),
    "Enterobacter_cloacae": (4.6, 5.5),
    "Enterococcus_faecalis": (2.7, 3.3),
    "Enterococcus_faecium": (2.6, 3.2),
    "Escherichia": (4.6, 5.6),
    "Haemophilus_influenzae": (1.8, 2.0),
    "Helicobacter_pylori": (1.55, 1.75),
    "Klebsiella_oxytoca": (5.6, 6.3),
    "Klebsiella_pneumoniae": (5.1, 5.9),
    "Pseudomonas_aeruginosa": (6.2, 7.2),
    "Salmonella": (4.6, 5.1),
    "Staphylococcus_aureus": (2.7, 3.0),
    "Staphylococcus_pseudintermedius": (2.5, 2.8),
    "Streptococcus_agalactiae": (1.9, 2.3),
    "Streptococcus_pneumoniae": (2.0, 2.2),
    "Streptococcus_pyogenes": (1.7, 1.95),
    "Vibrio_cholerae": (3.9, 4.2),
    "Vibrio_parahaemolyticus": (5.0, 5.3),
    "Vibrio_vulnificus": (4.9, 5.2),
}

# Colour-blind-safe categorical palette, used in this fixed order (validated set: adjacent pairs
# stay apart under protan/deutan/tritan simulation); a ninth category is folded into "Other".
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OTHER_COLOUR = "#b4b2a9"
# One-hue sequential ramp (light to dark blue) for the AMR heatmaps
SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"
FONT = 'system-ui, -apple-system, "Segoe UI", Roboto, sans-serif'


def warn(message):
    print(f"WARNING: {message}", file=sys.stderr)


# --------------------------------------------------------------------------------------------------
# Manifest
# --------------------------------------------------------------------------------------------------


def read_manifest(path):
    """Read the manifest into {tool: {sample: [Path, ...]}}; run-level files are under sample "-"."""
    path = Path(path)
    if path.suffix.lower() == ".json":
        with open(path) as handle:
            data = json.load(handle)
        rows = data.get("files", []) if isinstance(data, dict) else data
    else:
        with open(path, newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    manifest = defaultdict(lambda: defaultdict(list))
    for row in rows:
        sample = str(row.get("sample") or RUN_SAMPLE).strip() or RUN_SAMPLE
        tool = str(row.get("tool") or "").strip()
        file_name = str(row.get("file") or "").strip()
        if tool not in SAMPLE_TOOLS and tool not in RUN_TOOLS:
            warn(f"manifest: unknown tool '{tool}' for {file_name}, skipped")
            continue
        if not file_name:
            continue
        file_path = Path(file_name)
        if not file_path.is_absolute():
            file_path = path.parent / file_path
        if not file_path.exists():
            warn(f"manifest: {file_path} ({tool}, {sample}) does not exist, skipped")
            continue
        manifest[tool][sample].append(file_path)
    return manifest


# --------------------------------------------------------------------------------------------------
# Parsers: one per output format; each takes a path and returns plain values or a DataFrame
# --------------------------------------------------------------------------------------------------


def to_number(value):
    """Number from a table cell ("1,234", "12.5", "-", "") or None."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return None if isinstance(value, float) and math.isnan(value) else value
    text = str(value).strip().replace(",", "")
    if text in ("", "-", "NA", "nan", "None"):
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return int(number) if number.is_integer() and "." not in text and "e" not in text.lower() else number


def parse_fastp(path):
    """fastp or fastplong JSON: reads, bases and quality before and after filtering."""
    with open(path) as handle:
        data = json.load(handle)
    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    filtering = data.get("filtering_result", {})
    result = {
        "reads_before": before.get("total_reads"),
        "bases_before": before.get("total_bases"),
        "reads_after": after.get("total_reads"),
        "bases_after": after.get("total_bases"),
        "q20_rate_before": before.get("q20_rate"),
        "q30_rate_before": before.get("q30_rate"),
        "q20_rate_after": after.get("q20_rate"),
        "q30_rate_after": after.get("q30_rate"),
        "gc_after": after.get("gc_content"),
        "low_quality_reads": filtering.get("low_quality_reads"),
        "too_short_reads": filtering.get("too_short_reads"),
        "too_many_n_reads": filtering.get("too_many_N_reads"),
        "duplication_rate": data.get("duplication", {}).get("rate"),
    }
    if result["reads_before"]:
        result["mean_length_before"] = round(result["bases_before"] / result["reads_before"], 1)
    if result["reads_after"]:
        result["mean_length_after"] = round(result["bases_after"] / result["reads_after"], 1)
    return result


def parse_mash_genome_size(path):
    """Mash genome size file: one number, the estimated genome size in bases."""
    with open(path) as handle:
        return int(float(handle.read().split()[0]))


def parse_kraken2(path):
    """Kraken2 report: species rows (name, percent of all reads, reads) and the unclassified percent."""
    species = []
    unclassified = 0.0
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            percent, clade_reads, rank, name = float(fields[0]), int(fields[1]), fields[-3], fields[-1].strip()
            if rank == "U":
                unclassified = percent
            elif rank == "S":
                species.append({"name": name, "percent": percent, "reads": clade_reads})
    table = pd.DataFrame(species, columns=["name", "percent", "reads"])
    return table.sort_values("percent", ascending=False).reset_index(drop=True), unclassified


def parse_bracken(path):
    """Bracken table: species name, estimated reads and percent of the classified reads."""
    table = pd.read_csv(path, sep="\t")
    table = table.rename(columns={"new_est_reads": "reads"})
    table["percent"] = table["fraction_total_reads"].astype(float) * 100
    return table[["name", "percent", "reads"]].sort_values("percent", ascending=False).reset_index(drop=True)


QUAST_METRICS = OrderedDict(
    [
        ("# contigs (>= 0 bp)", "contigs_all"),
        ("# contigs", "contigs"),
        ("Total length (>= 0 bp)", "length_all"),
        ("Total length", "total_length"),
        ("Largest contig", "largest_contig"),
        ("N50", "n50"),
        ("L50", "l50"),
        ("GC (%)", "gc"),
        ("# N's per 100 kbp", "n_per_100kbp"),
    ]
)


def parse_quast(path):
    """QUAST per-sample TSV (metric, value). "# contigs" and "Total length" count contigs >= 500 bp."""
    metrics = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 1 and fields[0] in QUAST_METRICS:
                metrics[QUAST_METRICS[fields[0]]] = to_number(fields[1])
    return metrics


def parse_filter_stats(path):
    """seqkit stats of the assembly (first row) and of the annotation FASTA (second row)."""
    table = pd.read_csv(path, sep="\t")
    result = {}
    for prefix, index in (("assembly", 0), ("annotation", 1)):
        if len(table) > index:
            result[f"{prefix}_contigs"] = int(table.loc[index, "num_seqs"])
            result[f"{prefix}_length"] = int(table.loc[index, "sum_len"])
            result[f"{prefix}_n50"] = to_number(table.loc[index, "N50"])
    return result


def parse_checkm2(path):
    """CheckM2 quality report: completeness and contamination (percent) of the genome."""
    table = pd.read_csv(path, sep="\t")
    row = table.iloc[0]
    return {
        "completeness": float(row["Completeness"]),
        "contamination": float(row["Contamination"]),
        "model": row.get("Completeness_Model_Used"),
        "genome_size": to_number(row.get("Genome_Size")),
    }


def parse_mlst(path):
    """mlst line: file, scheme, ST, alleles. "-" as scheme or ST means none was found."""
    with open(path) as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if not lines:
        return {}
    fields = lines[-1]
    scheme = fields[1] if len(fields) > 1 and fields[1] != "-" else None
    st = fields[2] if len(fields) > 2 and fields[2] != "-" else None
    return {"scheme": scheme, "st": st, "alleles": " ".join(fields[3:])}


def parse_gtdbtk(path):
    """GTDB-Tk summary: one row per genome; returns {genome: {...}}."""
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    ani_column = next((c for c in ("closest_genome_ani", "fastani_ani") if c in table.columns), None)
    af_column = next((c for c in ("closest_genome_af", "fastani_af") if c in table.columns), None)
    result = {}
    for _, row in table.iterrows():
        classification = row.get("classification", "")
        ranks = dict(part.split("__", 1) for part in classification.split(";") if "__" in part)
        result[row["user_genome"]] = {
            "classification": classification,
            "species": ranks.get("s") or None,
            "genus": ranks.get("g") or None,
            "ani": to_number(row.get(ani_column)) if ani_column else None,
            "af": to_number(row.get(af_column)) if af_column else None,
            "method": row.get("classification_method") or None,
            "warnings": row.get("warnings") if row.get("warnings") not in ("", "N/A") else None,
        }
    return result


def parse_kleborate(path):
    """Kleborate 3 output for one genome; module prefixes (klebsiella_pneumo_complex__ ...) removed."""
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    row = table.iloc[0].to_dict()
    result = {}
    for column, value in row.items():
        short = column.split("__")[-1] if "__" in column else column
        # keep the first occurrence of a short name (e.g. ST from the KpSC MLST module)
        result.setdefault(short, value)
    return result


def parse_sistr(path):
    """SISTR serovar prediction (one row)."""
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return table.iloc[0].to_dict()


def parse_sccmec(path):
    """sccmec typing result (one row)."""
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return table.iloc[0].to_dict()


def sccmec_problem(row):
    """Why an sccmec call is not clear-cut, or None (one type from the target genes, or no SCCmec at all)."""
    sccmec_type = (row.get("type") or "").strip()
    if sccmec_type in ("", "-"):
        if row.get("mecA") == "+" or (row.get("targets") or "").strip() not in ("", "-"):
            return "mecA or cassette genes found, but no SCCmec type"
        return None  # no mecA and no cassette genes: no SCCmec (MSSA)
    if re.search(r"[,;]", sccmec_type):
        return f"several SCCmec types ({sccmec_type})"
    if row.get("target_comment"):
        # sccmec falls back to the closest cassette when the target genes do not match one type
        return "SCCmec type from cassette coverage only; the target genes gave no type"
    return None


AMRFINDER_COLUMNS = {
    # AMRFinderPlus 4.x name: (3.x name, short name)
    "Element symbol": ("Gene symbol", "gene"),
    "Element name": ("Sequence name", "name"),
    "Type": ("Element type", "type"),
    "Subtype": ("Element subtype", "subtype"),
    "Class": ("Class", "class"),
    "Subclass": ("Subclass", "subclass"),
    "Method": ("Method", "method"),
    "Scope": ("Scope", "scope"),
    "Contig id": ("Contig id", "contig"),
    "% Coverage of reference": ("% Coverage of reference sequence", "coverage"),
    "% Identity to reference": ("% Identity to reference sequence", "identity"),
    "Closest reference accession": ("Accession of closest sequence", "accession"),
}


def parse_amrfinderplus(path):
    """AMRFinderPlus hits with short column names (gene, class, subclass, type, subtype, identity ...)."""
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    columns = {}
    for new_name, (old_name, short) in AMRFINDER_COLUMNS.items():
        if new_name in table.columns:
            columns[new_name] = short
        elif old_name in table.columns:
            columns[old_name] = short
    table = table.rename(columns=columns)[list(columns.values())]
    for column in ("identity", "coverage"):
        if column in table:
            table[column] = pd.to_numeric(table[column], errors="coerce")
    return table


def parse_rgi(path):
    """RGI hits: gene (Best_Hit_ARO), drug class, mechanism, family, cut-off, identity, coverage."""
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    columns = {
        "Best_Hit_ARO": "gene",
        "Cut_Off": "cut_off",
        "Model_type": "model",
        "Drug Class": "drug_class",
        "Resistance Mechanism": "mechanism",
        "AMR Gene Family": "family",
        "Best_Identities": "identity",
        "Percentage Length of Reference Sequence": "coverage",
        "Contig": "contig",
        "SNPs_in_Best_Hit_ARO": "snps",
    }
    table = table.rename(columns=columns)[[c for c in columns.values() if c in table.rename(columns=columns)]]
    for column in ("identity", "coverage"):
        if column in table:
            table[column] = pd.to_numeric(table[column], errors="coerce")
    if "contig" in table:
        # RGI appends the ORF number to the contig name (contig_3_12)
        table["contig"] = table["contig"].str.replace(r"_\d+$", "", regex=True)
    return table


def parse_genomad(path):
    """geNomad plasmid or virus summary: one row per sequence."""
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def parse_mobsuite_contigs(path):
    """MOB-suite contig report: one row per contig with molecule_type chromosome or plasmid."""
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def parse_mobsuite_mobtyper(path):
    """MOB-typer results: one row per reconstructed plasmid."""
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def parse_plasmid_consensus(path):
    """Plasmid consensus: one row per contig with each tool's evidence and the call (the manifest gives the sample)."""
    return pd.read_csv(path, sep="\t", dtype=str).fillna("").drop(columns="sample", errors="ignore")


def parse_isescan(path):
    """ISEScan table: one row per IS element (seqID, family, isLen, type c complete or p partial, ...)."""
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for name in ("isBegin", "isEnd", "isLen", "ncopy4is"):
        if name in table:
            table[name] = pd.to_numeric(table[name], errors="coerce")
    return table


def parse_sample_status(path):
    """pipeline_info/sample_status.tsv (Phase 7 columns)."""
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def parse_software_versions(path):
    """software_versions.yml: {process: {tool: version}} -> rows (tool, version, processes)."""
    if yaml is None:
        raise RuntimeError("PyYAML is not installed")
    with open(path) as handle:
        # BaseLoader keeps every value as text: safe_load reads a version such as 2.20 as the number 2.2
        data = yaml.load(handle, Loader=yaml.BaseLoader) or {}
    by_tool = OrderedDict()
    for process, tools in data.items():
        if not isinstance(tools, dict):
            continue
        for tool, version in tools.items():
            key = (str(tool), str(version))
            by_tool.setdefault(key, []).append(str(process))
    rows = [
        {"tool": tool, "version": version, "processes": processes} for (tool, version), processes in by_tool.items()
    ]
    return sorted(rows, key=lambda row: (row["processes"] != ["Workflow"], row["tool"].lower()))


def parse_citations(path):
    """CITATIONS.md: [(name, url, reference text)] for every "[Name](url)" heading or list item."""
    entries = []
    current = None
    with open(path) as handle:
        for line in handle:
            match = re.match(r"^\s*(?:##|-)\s*\[([^\]]+)\]\(([^)]+)\)", line)
            if match:
                current = {"name": match.group(1), "url": match.group(2), "text": ""}
                entries.append(current)
                continue
            quote = re.match(r"^\s*>\s?(.*)", line)
            if quote and current is not None:
                current["text"] = (current["text"] + " " + quote.group(1)).strip()
    return entries


# --------------------------------------------------------------------------------------------------
# Run data: everything parsed, keyed by sample
# --------------------------------------------------------------------------------------------------


def natural_key(text):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", str(text))]


class RunData:
    """Parsed outputs of one run. Attributes stay empty when their inputs were not given."""

    def __init__(self, manifest, run_info):
        self.run_info = run_info or {}
        self.reads = {}  # sample -> fastp/fastplong dict
        self.read_tool = {}  # sample -> "fastp" or "fastplong"
        self.genome_size = {}  # sample -> Mash estimate (bases)
        self.taxa = {}  # sample -> DataFrame(name, percent, reads)
        self.taxa_source = {}  # sample -> "Bracken" or "Kraken2"
        self.unclassified = {}  # sample -> percent (Kraken2)
        self.quast = {}
        self.filter_stats = {}
        self.checkm2 = {}
        self.mlst = {}
        self.gtdbtk = {}
        self.kleborate = {}
        self.sistr = {}
        self.sccmec = {}
        self.amrfinder = {}  # sample -> DataFrame
        self.rgi = {}
        self.genomad_plasmid = {}
        self.genomad_virus = {}
        self.mobsuite_contigs = {}
        self.mobsuite_mobtyper = {}
        self.plasmid_consensus = {}
        self.isescan = {}  # sample -> DataFrame, empty when ISEScan found no IS element
        self.status = None  # DataFrame
        self.versions = []
        self.citations = []
        self.problems = []  # (tool, file, error) for the methods section
        self._load(manifest)

    def _parse(self, tool, parser, path):
        try:
            return parser(path)
        except Exception as error:  # a broken file must not stop the report
            warn(f"{tool}: could not parse {path}: {error}")
            self.problems.append((tool, str(path), str(error)))
            return None

    def _load(self, manifest):
        per_sample = {
            "mash": (parse_mash_genome_size, self.genome_size),
            "quast": (parse_quast, self.quast),
            "filter_stats": (parse_filter_stats, self.filter_stats),
            "checkm2": (parse_checkm2, self.checkm2),
            "mlst": (parse_mlst, self.mlst),
            "kleborate": (parse_kleborate, self.kleborate),
            "sistr": (parse_sistr, self.sistr),
            "sccmec": (parse_sccmec, self.sccmec),
        }
        concatenated = {
            "amrfinderplus": (parse_amrfinderplus, self.amrfinder),
            "rgi": (parse_rgi, self.rgi),
            "genomad_plasmid": (parse_genomad, self.genomad_plasmid),
            "genomad_virus": (parse_genomad, self.genomad_virus),
            "mobsuite_contigs": (parse_mobsuite_contigs, self.mobsuite_contigs),
            "mobsuite_mobtyper": (parse_mobsuite_mobtyper, self.mobsuite_mobtyper),
            "plasmid_consensus": (parse_plasmid_consensus, self.plasmid_consensus),
            "isescan": (parse_isescan, self.isescan),
        }
        for tool in ("fastp", "fastplong"):
            for sample, paths in manifest.get(tool, {}).items():
                parsed = self._parse(tool, parse_fastp, paths[0])
                if parsed is not None:
                    self.reads[sample] = parsed
                    self.read_tool[sample] = tool
        for sample, paths in manifest.get("bracken", {}).items():
            parsed = self._parse("bracken", parse_bracken, paths[0])
            if parsed is not None:
                self.taxa[sample] = parsed
                self.taxa_source[sample] = "Bracken"
        for sample, paths in manifest.get("kraken2", {}).items():
            parsed = self._parse("kraken2", parse_kraken2, paths[0])
            if parsed is None:
                continue
            species, unclassified = parsed
            self.unclassified[sample] = unclassified
            if sample not in self.taxa:
                self.taxa[sample] = species
                self.taxa_source[sample] = "Kraken2"
        for tool, (parser, target) in per_sample.items():
            for sample, paths in manifest.get(tool, {}).items():
                parsed = self._parse(tool, parser, paths[0])
                if parsed is not None:
                    target[sample] = parsed
        for tool, (parser, target) in concatenated.items():
            for sample, paths in manifest.get(tool, {}).items():
                frames = [frame for frame in (self._parse(tool, parser, path) for path in paths) if frame is not None]
                if frames:
                    target[sample] = pd.concat(frames, ignore_index=True)
        for path in manifest.get("gtdbtk", {}).get(RUN_SAMPLE, []):
            parsed = self._parse("gtdbtk", parse_gtdbtk, path)
            if parsed:
                self.gtdbtk.update(parsed)
        for path in manifest.get("sample_status", {}).get(RUN_SAMPLE, [])[:1]:
            self.status = self._parse("sample_status", parse_sample_status, path)
        for path in manifest.get("software_versions", {}).get(RUN_SAMPLE, [])[:1]:
            self.versions = self._parse("software_versions", parse_software_versions, path) or []
        for path in manifest.get("citations", {}).get(RUN_SAMPLE, [])[:1]:
            self.citations = self._parse("citations", parse_citations, path) or []

        samples = set()
        for tool in SAMPLE_TOOLS:
            samples.update(sample for sample in manifest.get(tool, {}) if sample != RUN_SAMPLE)
        samples.update(self.gtdbtk)
        if self.status is not None and "sample" in self.status:
            samples.update(self.status["sample"])
        # natural order of the names (S2 before S10); sample_status.tsv is sorted as plain text
        self.samples = sorted(samples, key=natural_key)

    # ---- per-sample derived values

    def status_row(self, sample):
        if self.status is None or "sample" not in self.status:
            return {}
        rows = self.status[self.status["sample"] == sample]
        return rows.iloc[0].to_dict() if len(rows) else {}

    def organism(self, sample):
        """AMRFinderPlus organism for the sample's MLST scheme (run-info map), or None."""
        scheme = self.mlst.get(sample, {}).get("scheme")
        return self.run_info.get("mlst_organisms", {}).get(scheme) if scheme else None

    def expected_length(self, sample):
        organism = self.organism(sample)
        return EXPECTED_GENOME_MB.get(organism) if organism else None

    def checkm2_group(self, sample):
        """CHECKM2_GROUPS key for a genome outside the CheckM2 limits, or None."""
        checkm2 = self.checkm2.get(sample)
        if not checkm2:
            return None
        low = checkm2["completeness"] < CHECKM2_MIN_COMPLETENESS
        high = checkm2["contamination"] > CHECKM2_MAX_CONTAMINATION
        return "both" if low and high else "low_completeness" if low else "high_contamination" if high else None

    def amr_genes(self, sample):
        """Acquired AMR genes, as in the heatmap; point mutations (subtype POINT) are listed separately."""
        table = self.amrfinder.get(sample)
        if table is None or table.empty or "type" not in table:
            return []
        acquired = (table["type"] == "AMR") & (column(table, "subtype") != "POINT")
        return sorted(set(table.loc[acquired, "gene"]), key=str.lower)

    def warnings(self, sample):
        """Plain-language warnings for the sample card and the TSV."""
        notes = []
        status = self.status_row(sample)
        if status.get("status") in ("skipped", "failed"):
            notes.append(f"{status['status'].capitalize()} at {status.get('stage') or '?'}: {status.get('reason')}")
        elif status.get("reason"):
            notes.append(status["reason"])
        checkm2 = self.checkm2.get(sample)
        if checkm2:
            if checkm2["completeness"] < CHECKM2_MIN_COMPLETENESS:
                notes.append(f"CheckM2 completeness {checkm2['completeness']:.1f} % (< {CHECKM2_MIN_COMPLETENESS:g} %)")
            if checkm2["contamination"] > CHECKM2_MAX_CONTAMINATION:
                notes.append(
                    f"CheckM2 contamination {checkm2['contamination']:.1f} % (> {CHECKM2_MAX_CONTAMINATION:g} %)"
                )
        expected = self.expected_length(sample)
        length = self.quast.get(sample, {}).get("total_length")
        if expected and length:
            low, high = expected
            if not low * 0.9 <= length / 1e6 <= high * 1.1:
                notes.append(
                    f"Assembly length {length / 1e6:.2f} Mb is outside the usual {low:.1f}-{high:.1f} Mb "
                    f"for {self.organism(sample).replace('_', ' ')}"
                )
        taxa = self.taxa.get(sample)
        if taxa is not None and len(taxa) > 1 and taxa.iloc[1]["percent"] > SECOND_SPECIES_WARN:
            notes.append(
                f"{self.taxa_source[sample]}: second species {taxa.iloc[1]['name']} at "
                f"{taxa.iloc[1]['percent']:.1f} % of reads (possible contamination)"
            )
        disagreement = self.species_disagreement(sample)
        if disagreement:
            notes.append(disagreement)
        if sample in self.sccmec and sccmec_problem(self.sccmec[sample]):
            notes.append(sccmec_problem(self.sccmec[sample]))
        return notes

    def species_checked(self, sample):
        """True when a second species call (GTDB-Tk or reads) exists to compare with MLST."""
        taxa = self.taxa.get(sample)
        return bool(self.gtdbtk.get(sample, {}).get("genus")) or (taxa is not None and len(taxa) > 0)

    def mlst_genus(self, sample):
        organism = self.organism(sample)
        return organism.split("_")[0] if organism else None

    def species_disagreement(self, sample):
        """A note when the MLST genus differs from the GTDB-Tk or read-based genus."""
        mlst_genus = self.mlst_genus(sample)
        if not mlst_genus:
            return None
        others = []
        gtdb = self.gtdbtk.get(sample, {}).get("genus")
        if gtdb:
            others.append(("GTDB-Tk", gtdb.split("_")[0]))
        taxa = self.taxa.get(sample)
        if taxa is not None and len(taxa):
            others.append((self.taxa_source[sample], taxa.iloc[0]["name"].split()[0]))
        differing = [f"{source} {genus}" for source, genus in others if genus != mlst_genus]
        if differing:
            return f"MLST scheme {self.mlst[sample]['scheme']} ({mlst_genus}) disagrees with " + ", ".join(differing)
        return None


# --------------------------------------------------------------------------------------------------
# Formatting helpers
# --------------------------------------------------------------------------------------------------


def esc(value):
    return html.escape("" if value is None else str(value))


def fmt_int(value):
    return "" if value is None else f"{int(value):,}"


def fmt_float(value, digits=1):
    return "" if value is None or (isinstance(value, float) and math.isnan(value)) else f"{value:,.{digits}f}"


def fmt_mb(value, digits=2):
    return "" if value is None else f"{value / 1e6:,.{digits}f}"


def fmt_pct(rate, digits=1):
    return "" if rate is None else f"{rate * 100:.{digits}f}"


def sample_link(sample):
    return f'<a href="#sample-{slug(sample)}">{esc(sample)}</a>'


def slug(text):
    """HTML id part for a sample name: other characters, "_" included, become _<hex>_, so two names never share one."""
    return re.sub(r"[^A-Za-z0-9-]", lambda match: f"_{ord(match.group()):x}_", str(text))


def species_name(organism):
    return organism.replace("_", " ") if organism else ""


def render_table(columns, rows, table_id=None, numeric=(), raw=(), max_height=None, caption=None):
    """HTML table. columns: [(key, label)]; values in `raw` keys are already HTML."""
    if not rows:
        return ""
    attributes = f' id="{esc(table_id)}"' if table_id else ""
    if max_height is None and len(rows) > TABLE_SCROLL_ROWS:
        max_height = TABLE_SCROLL_HEIGHT
    style = f' style="max-height:{max_height}px"' if max_height else ""
    parts = [f'<div class="table-wrap"{style}><table class="data sortable"{attributes}>']
    if caption:
        parts.append(f"<caption>{esc(caption)}</caption>")
    parts.append("<thead><tr>")
    for key, label in columns:
        cls = ' class="num"' if key in numeric else ""
        parts.append(f'<th{cls} scope="col">{esc(label)}</th>')
    parts.append("</tr></thead><tbody>")
    for row in rows:
        parts.append("<tr>")
        for key, _ in columns:
            value = row.get(key)
            cell = value if key in raw else esc(value)
            if isinstance(cell, str) and key not in raw:
                cell = cell.replace("/", "/<wbr>")
            cls = ' class="num"' if key in numeric else ""
            parts.append(f"<td{cls}>{'' if cell is None else cell}</td>")
        parts.append("</tr>")
    parts.append("</tbody></table></div>")
    return "".join(parts)


def figure_html(fig, height, legend_y=1.0):
    """Plotly figure as a div; plotly.js itself is inlined once in the page head."""
    fig.update_layout(
        height=height,
        font=dict(family=FONT, size=12, color=INK),
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        hoverlabel=dict(font=dict(family=FONT)),
        legend=dict(orientation="h", yanchor="bottom", y=legend_y, xanchor="left", x=0, title_text=""),
    )
    fig.update_xaxes(gridcolor=GRID, linecolor="#c3c2b7", zerolinecolor="#c3c2b7", title_font=dict(color=INK_SECONDARY))
    fig.update_yaxes(gridcolor=GRID, linecolor="#c3c2b7", zerolinecolor="#c3c2b7", title_font=dict(color=INK_SECONDARY))
    config = {"displaylogo": False, "responsive": True, "modeBarButtonsToRemove": ["select2d", "lasso2d"]}
    return '<div class="figure">' + fig.to_html(include_plotlyjs=False, full_html=False, config=config) + "</div>"


def bar_height(n_rows, per_row=22, minimum=260, extra=120):
    """Height of a horizontal bar chart: one row per sample, so 200 samples stay readable (scroll)."""
    return max(minimum, n_rows * per_row + extra)


def outliers(values, low=None, high=None):
    """Samples outside Tukey's fences (1.5 x IQR beyond the quartiles) or outside the given limits."""
    series = pd.Series(values, dtype=float).dropna()
    if series.empty:
        return set()
    q1, q3 = series.quantile(0.25), series.quantile(0.75)
    spread = 1.5 * (q3 - q1)
    outside = (series < q1 - spread) | (series > q3 + spread) if spread > 0 else pd.Series(False, index=series.index)
    if low is not None:
        outside |= series < low
    if high is not None:
        outside |= series > high
    return set(series[outside].index)


class Panel:
    """One metric of a distribution figure: {sample: value} with optional limits and extra flagged samples."""

    def __init__(self, title, unit, values, low=None, high=None, flagged=(), digits=1):
        self.title, self.unit, self.digits = title, unit, digits
        self.values = {s: v for s, v in values.items() if v is not None and not pd.isna(v)}
        self.low, self.high = low, high
        self.outliers = outliers(self.values, low, high) | {s for s in flagged if s in self.values}


def outlier_labels(values, flagged):
    """Labels for the most extreme outliers: points closer than 2 % of the axis share one label, 'name +N'."""
    spread = (max(values.values()) - min(values.values())) or 1
    groups = []  # [first sample, value, count]
    for s in flagged:
        group = next((g for g in groups if abs(values[s] - g[1]) <= 0.02 * spread), None)
        if group:
            group[2] += 1
        elif len(groups) < OUTLIER_LABELS:
            groups.append([s, values[s], 1])
    return {s: f"{s} +{count - 1}" if count > 1 else s for s, _, count in groups}


def distribution_html(panels):
    """Large runs: one box with every sample as a point per metric; outliers in red, the most extreme named."""
    panels = [panel for panel in panels if panel.values]
    if not panels:
        return ""
    fig = make_subplots(
        rows=len(panels), cols=1, subplot_titles=[p.title for p in panels], vertical_spacing=0.6 / len(panels)
    )
    for row, panel in enumerate(panels, 1):
        samples = list(panel.values)
        values = [panel.values[s] for s in samples]
        hover = f"%{{text}}<br>%{{x:,.{panel.digits}f}} {panel.unit}<extra></extra>"
        fig.add_trace(
            go.Box(
                x=values,
                y=[panel.title] * len(values),
                text=samples,
                orientation="h",
                boxpoints="all",
                jitter=0.6,
                pointpos=0,
                marker=dict(size=5, color=PALETTE[0], opacity=0.55),
                line=dict(color=INK_SECONDARY, width=1),
                fillcolor="rgba(0,0,0,0)",
                name="Samples",
                showlegend=False,
                hovertemplate=hover,
            ),
            row=row,
            col=1,
        )
        flagged = sorted(panel.outliers, key=lambda s: abs(panel.values[s] - pd.Series(values).median()), reverse=True)
        labels = outlier_labels(panel.values, flagged)
        if flagged:
            fig.add_trace(
                go.Scatter(
                    x=[panel.values[s] for s in flagged],
                    y=[panel.title] * len(flagged),
                    text=[labels.get(s, "") for s in flagged],
                    customdata=flagged,
                    mode="markers+text",
                    textposition="top center",
                    textfont=dict(size=10, color=INK_SECONDARY),
                    marker=dict(size=9, color=PALETTE[7], symbol="diamond", line=dict(color=SURFACE, width=1)),
                    name="Outliers (most extreme named)",
                    showlegend=not any(t.legendgroup == "outliers" for t in fig.data),
                    legendgroup="outliers",
                    hovertemplate=f"%{{customdata}}<br>%{{x:,.{panel.digits}f}} {panel.unit}<extra></extra>",
                ),
                row=row,
                col=1,
            )
        for limit in (panel.low, panel.high):
            if limit is not None:
                fig.add_vline(x=limit, line_dash="dash", line_color=INK_SECONDARY, row=row, col=1)
        fig.update_xaxes(title_text=panel.unit, row=row, col=1)
        fig.update_yaxes(showticklabels=False, row=row, col=1)
    fig.update_layout(margin=dict(l=10, r=20, t=70, b=40))
    fig.update_annotations(font=dict(size=13, color=INK))
    return figure_html(fig, 150 * len(panels) + 90, legend_y=1.0 + 40 / (150 * len(panels) + 90))


# --------------------------------------------------------------------------------------------------
# Sections
# --------------------------------------------------------------------------------------------------


class Section:
    def __init__(self, key, title, intro, body):
        self.key, self.title, self.intro, self.body = key, title, intro, body

    def html(self):
        intro = "".join(f"<p>{paragraph}</p>" for paragraph in self.intro)
        return (
            f'<section id="{self.key}"><h2>{esc(self.title)}</h2>'
            f'<div class="intro">{intro}</div>{self.body}</section>'
        )


def section_overview(run):
    info = run.run_info
    facts = []
    for key, label in (
        ("pipeline_version", "Pipeline version"),
        ("date", "Run date"),
        ("seq_type", "Sequencing type"),
        ("run_name", "Run name"),
        ("nextflow_version", "Nextflow"),
    ):
        if info.get(key):
            facts.append((label, info[key]))
    facts.append(("Samples", str(len(run.samples))))
    body = '<dl class="facts">' + "".join(f"<dt>{esc(k)}</dt><dd>{esc(v)}</dd>" for k, v in facts) + "</dl>"

    if run.status is not None and len(run.status):
        counts = run.status["status"].value_counts()
        chips = "".join(
            f'<span class="chip status-{esc(name)}">{esc(name)}: {int(counts.get(name, 0))}</span>'
            for name in ("pass", "skipped", "failed")
        )
        body += f'<p class="chips">{chips}</p>'
        not_passed = run.status[run.status["status"] != "pass"]
        if len(not_passed):
            rows = [
                {
                    "sample": sample_link(row["sample"]),
                    "status": row["status"],
                    "stage": row.get("stage", ""),
                    "reason": row.get("reason", ""),
                }
                for _, row in not_passed.iterrows()
            ]
            body += "<h3>Samples that did not pass</h3>" + render_table(
                [("sample", "Sample"), ("status", "Status"), ("stage", "Stage"), ("reason", "Reason")],
                rows,
                raw=("sample",),
            )

    parameters = info.get("parameters") or {}
    if parameters:
        rows = [{"name": name, "value": value} for name, value in parameters.items()]
        body += "<h3>Key thresholds</h3>" + render_table([("name", "Parameter"), ("value", "Value")], rows)
    databases = info.get("databases") or {}
    if databases:
        versions = {
            DATABASE_VERSION_ENTRIES[row["tool"]]: row["version"]
            for row in run.versions
            if row["tool"] in DATABASE_VERSION_ENTRIES
        }
        versions.update(info.get("database_versions") or {})
        rows = [
            {"name": name, "version": versions.get(name) or "not recorded", "path": path}
            for name, path in databases.items()
        ]
        body += "<h3>Databases</h3>" + render_table(
            [("name", "Parameter"), ("version", "Version"), ("path", "Path")], rows
        )

    body += (
        '<p class="note">Reading the charts: with up to '
        f"{LARGE_RUN} samples, bar charts have one row per sample, in name order (S2 before S10); with more, each "
        "metric is a box with one point per sample, and only the outliers (outside 1.5 times the interquartile "
        "range or a QC limit) are marked, the most extreme by name. Heatmaps have one column per sample and scroll "
        "sideways when they are wider than the page. Hover a bar, point or cell for its value; click a table header"
        " to sort; long tables scroll inside their box; sample names open the sample cards. Colours follow one "
        "colour-blind-safe palette, and every colour is also named in a legend or label.</p>"
    )
    intro = [
        "This page summarises one run of MGAP (Microbial Genome Analysis Pipeline): read quality, assembly, "
        "genome quality, typing, antimicrobial resistance and mobile elements for every sample.",
        "Each sample is <b>pass</b> (assembled and annotated), <b>skipped</b> (too few bases after read "
        "filtering to assemble) or <b>failed</b> (a step failed; the stage and reason are listed). "
        "Sections appear only for the tools that ran; every sample has a card at the end with its key numbers "
        "and warnings.",
    ]
    return Section("overview", "Run overview", intro, body)


def section_read_qc(run):
    samples = [s for s in run.samples if s in run.reads]
    taxa_samples = [s for s in run.samples if s in run.taxa]
    if not samples and not taxa_samples:
        return None
    body = ""
    min_bases = to_number((run.run_info.get("parameters") or {}).get("min_bases"))
    if len(samples) > LARGE_RUN:
        reads = run.reads
        body += distribution_html(
            [
                Panel(
                    "Reads kept",
                    "million reads",
                    {s: (reads[s].get("reads_after") or 0) / 1e6 for s in samples},
                    digits=2,
                ),
                Panel(
                    "Bases kept",
                    "Mb",
                    {s: (reads[s].get("bases_after") or 0) / 1e6 for s in samples},
                    low=min_bases / 1e6 if min_bases else None,
                ),
                Panel(
                    "Bases kept by read filtering",
                    "%",
                    {
                        s: 100 * reads[s]["bases_after"] / reads[s]["bases_before"]
                        for s in samples
                        if reads[s].get("bases_before") and reads[s].get("bases_after") is not None
                    },
                ),
                Panel(
                    "Q30 after filtering",
                    "%",
                    {
                        s: 100 * reads[s]["q30_rate_after"]
                        for s in samples
                        if reads[s].get("q30_rate_after") is not None
                    },
                ),
            ]
        )
    elif samples:
        fig = make_subplots(rows=1, cols=2, shared_yaxes=True, horizontal_spacing=0.04)
        ordered = list(reversed(samples))  # first sample at the top
        for col, (key, scale, unit) in enumerate((("reads", 1e6, "Reads (millions)"), ("bases", 1e6, "Bases (Mb)")), 1):
            for slot, (stage, label) in enumerate((("before", "Before filtering"), ("after", "After filtering"))):
                values = [(run.reads[s].get(f"{key}_{stage}") or 0) / scale for s in ordered]
                fig.add_trace(
                    go.Bar(
                        y=ordered,
                        x=values,
                        orientation="h",
                        name=label,
                        marker_color=PALETTE[slot],
                        legendgroup=label,
                        showlegend=col == 1,
                        hovertemplate=f"%{{y}}<br>{label}: %{{x:,.2f}} {unit.split('(')[1][:-1]}<extra></extra>",
                    ),
                    row=1,
                    col=col,
                )
            fig.update_xaxes(title_text=unit, row=1, col=col)
        if min_bases:
            fig.add_vline(x=min_bases / 1e6, line_dash="dash", line_color=INK_SECONDARY, row=1, col=2)
            fig.add_annotation(
                x=min_bases / 1e6,
                xref="x2",
                y=1,
                yref="paper",
                yanchor="bottom",
                xanchor="left",
                showarrow=False,
                text=f"--min_bases {min_bases / 1e6:g} Mb",
                font=dict(color=INK_SECONDARY, size=11),
            )
        # first trace (before) drawn below the second within a group; list it first in the legend too
        fig.update_layout(
            barmode="group",
            bargap=0.25,
            bargroupgap=0.08,
            margin=dict(l=10, r=20, t=50, b=50),
            legend_traceorder="normal",
        )
        body += figure_html(fig, bar_height(len(samples), per_row=26))

    if samples:
        rows = []
        for s in samples:
            r = run.reads[s]
            rows.append(
                {
                    "sample": sample_link(s),
                    "tool": run.read_tool[s],
                    "reads_before": fmt_int(r.get("reads_before")),
                    "reads_after": fmt_int(r.get("reads_after")),
                    "bases_before": fmt_mb(r.get("bases_before")),
                    "bases_after": fmt_mb(r.get("bases_after")),
                    "kept": fmt_float(100 * r["bases_after"] / r["bases_before"]) if r.get("bases_before") else "",
                    "mean_length": fmt_float(r.get("mean_length_before")),
                    "coverage": (
                        fmt_float(r["bases_after"] / run.genome_size[s])
                        if run.genome_size.get(s) and r.get("bases_after")
                        else ""
                    ),
                    "q30_before": fmt_pct(r.get("q30_rate_before")),
                    "q30_after": fmt_pct(r.get("q30_rate_after")),
                    "gc": fmt_pct(r.get("gc_after")),
                    "duplication": fmt_pct(r.get("duplication_rate"), 2),
                }
            )
        columns = [
            ("sample", "Sample"),
            ("tool", "Tool"),
            ("reads_before", "Reads in"),
            ("reads_after", "Reads kept"),
            ("bases_before", "Mb in"),
            ("bases_after", "Mb kept"),
            ("kept", "Bases kept (%)"),
            ("mean_length", "Mean read length in (bp)"),
            ("coverage", "Coverage kept (x)"),
            ("q30_before", "Q30 in (%)"),
            ("q30_after", "Q30 kept (%)"),
            ("gc", "GC (%)"),
            ("duplication", "Duplication (%)"),
        ]
        if not any(run.reads[s].get("duplication_rate") is not None for s in samples):
            columns = [c for c in columns if c[0] != "duplication"]
        if not run.genome_size:
            columns = [c for c in columns if c[0] != "coverage"]
        body += render_table(columns, rows, numeric={c for c, _ in columns[2:]}, raw=("sample",))

    if taxa_samples:
        body += "<h3>Species in the reads</h3>"
        body += taxa_chart(run, taxa_samples)
        rows = []
        for s in taxa_samples:
            taxa = run.taxa[s]
            top = taxa.iloc[0] if len(taxa) else None
            second = taxa.iloc[1] if len(taxa) > 1 else None
            rows.append(
                {
                    "sample": sample_link(s),
                    "source": run.taxa_source[s],
                    "top": "" if top is None else f"<i>{esc(top['name'])}</i>",
                    "top_pct": "" if top is None else fmt_float(top["percent"]),
                    "second": "" if second is None else f"<i>{esc(second['name'])}</i>",
                    "second_pct": "" if second is None else fmt_float(second["percent"]),
                    "unclassified": fmt_float(run.unclassified.get(s)),
                }
            )
        columns = [
            ("sample", "Sample"),
            ("source", "Source"),
            ("top", "Top species"),
            ("top_pct", "Top (%)"),
            ("second", "Second species"),
            ("second_pct", "Second (%)"),
            ("unclassified", "Unclassified by Kraken2 (% of reads)"),
        ]
        body += render_table(
            columns, rows, numeric=("top_pct", "second_pct", "unclassified"), raw=("sample", "top", "second")
        )

    intro = [
        "Reads are cleaned before assembly: fastp (Illumina) or fastplong (Oxford Nanopore) trims adapters and "
        "removes short and low-quality reads. The chart compares reads and bases before and after this step; "
        "a large loss, or a sample far below the others, points to a sequencing or library problem. "
        "Q30 is the share of bases with a quality of 30 or more (one error in 1,000 bases). For Illumina, "
        "coverage is the bases kept divided by Mash's genome size estimate; deeper samples are subsampled "
        "before assembly.",
    ]
    if taxa_samples:
        intro.append(
            "Kraken2 assigns each read to a taxon and Bracken re-estimates the species shares; an isolate should "
            "show one species with most of the reads. A second species above "
            f"{SECOND_SPECIES_WARN:g} % may mean a mixed or contaminated culture. "
            + (
                f"The first chart counts the samples by their top species ({TOP_TAXA} species by name, the rest as "
                "Other); the boxes below give each sample's share of its top and second species."
                if len(taxa_samples) > LARGE_RUN
                else f"The chart shows the {TOP_TAXA} most abundant species of the run by name and sums the rest "
                "as Other."
            )
        )
    return Section("read-qc", "Read quality and species in the reads", intro, body)


def top_taxa(run, samples, count=TOP_TAXA):
    """Taxa named in the chart: each sample's top species first (most samples first), then the largest
    summed share, so one sample's many minor species do not push out another sample's main species."""
    leaders = defaultdict(int)
    totals = defaultdict(float)
    for s in samples:
        taxa = run.taxa[s]
        if len(taxa):
            leaders[taxa.iloc[0]["name"]] += 1
        for _, row in taxa.iterrows():
            totals[row["name"]] += row["percent"]
    ranked = sorted(totals, key=lambda name: (-leaders.get(name, 0), -totals[name], name))
    return ranked[:count]


def taxa_chart(run, samples):
    if len(samples) > LARGE_RUN:
        return taxa_distribution(run, samples)
    top = top_taxa(run, samples)
    ordered = list(reversed(samples))
    fig = go.Figure()
    for index, name in enumerate(top):
        values = []
        for s in ordered:
            taxa = run.taxa[s]
            match = taxa.loc[taxa["name"] == name, "percent"]
            values.append(float(match.iloc[0]) if len(match) else 0.0)
        fig.add_trace(
            go.Bar(
                y=ordered,
                x=values,
                orientation="h",
                name=name,
                marker=dict(color=PALETTE[index], line=dict(color=SURFACE, width=1)),
                hovertemplate=f"%{{y}}<br><i>{esc(name)}</i>: %{{x:.2f}} %<extra></extra>",
            )
        )
    other = []
    for s in ordered:
        taxa = run.taxa[s]
        other.append(float(taxa.loc[~taxa["name"].isin(top), "percent"].sum()))
    fig.add_trace(
        go.Bar(
            y=ordered,
            x=other,
            orientation="h",
            name="Other species",
            marker=dict(color=OTHER_COLOUR, line=dict(color=SURFACE, width=1)),
            hovertemplate="%{y}<br>Other species: %{x:.2f} %<extra></extra>",
        )
    )
    # Kraken2 percentages are of all reads: the rest are unclassified or assigned above species level
    rest = [
        max(0.0, 100 - float(run.taxa[s]["percent"].sum())) if run.taxa_source[s] == "Kraken2" else 0.0 for s in ordered
    ]
    if any(value > 0.05 for value in rest):
        fig.add_trace(
            go.Bar(
                y=ordered,
                x=rest,
                orientation="h",
                name="Unclassified or not to species",
                marker=dict(color="#e1e0d9", line=dict(color=SURFACE, width=1)),
                hovertemplate="%{y}<br>Unclassified or not to species: %{x:.2f} %<extra></extra>",
            )
        )
    sources = sorted(set(run.taxa_source[s] for s in samples))
    basis = (
        "of classified reads (Bracken)"
        if sources == ["Bracken"]
        else "of reads (Kraken2) or classified reads (Bracken)"
    )
    if sources == ["Kraken2"]:
        basis = "of reads (Kraken2)"
    fig.update_layout(barmode="stack", bargap=0.3, margin=dict(l=10, r=20, t=60, b=50), legend_traceorder="normal")
    fig.update_xaxes(title_text=f"Share {basis} (%)", range=[0, 100])
    return figure_html(fig, bar_height(len(samples), extra=150))


def taxa_distribution(run, samples):
    """Large runs: samples per top species (one bar per species), then the top and second species shares."""
    leaders = pd.Series([run.taxa[s].iloc[0]["name"] for s in samples if len(run.taxa[s])], dtype=str)
    counts = leaders.value_counts()
    named = counts.iloc[:TOP_TAXA]
    if len(counts) > TOP_TAXA:
        named = pd.concat([named, pd.Series({"Other species": int(counts.iloc[TOP_TAXA:].sum())})])
    ordered = list(reversed(named.index))
    fig = go.Figure(
        go.Bar(
            y=ordered,
            x=[int(named[name]) for name in ordered],
            orientation="h",
            marker_color=[OTHER_COLOUR if name == "Other species" else PALETTE[0] for name in ordered],
            hovertemplate="%{y}: %{x} samples<extra></extra>",
            showlegend=False,
        )
    )
    fig.update_layout(bargap=0.3, margin=dict(l=10, r=20, t=30, b=50))
    fig.update_xaxes(title_text="Samples whose top species in the reads is this one")
    second = {s: float(run.taxa[s].iloc[1]["percent"]) for s in samples if len(run.taxa[s]) > 1}
    return figure_html(fig, bar_height(len(ordered), per_row=28)) + distribution_html(
        [
            Panel(
                "Share of the top species",
                "%",
                {s: float(run.taxa[s].iloc[0]["percent"]) for s in samples if len(run.taxa[s])},
            ),
            Panel("Share of the second species", "%", second, high=SECOND_SPECIES_WARN),
        ]
    )


def section_assembly(run):
    samples = [s for s in run.samples if s in run.quast or s in run.filter_stats]
    if not samples:
        return None
    body = ""
    quast_samples = [s for s in samples if s in run.quast]
    if len(quast_samples) > LARGE_RUN:
        unusual = set()
        for s in quast_samples:
            expected, length = run.expected_length(s), run.quast[s].get("total_length")
            if expected and length and not expected[0] * 0.9 <= length / 1e6 <= expected[1] * 1.1:
                unusual.add(s)
        body += distribution_html(
            [
                Panel(
                    "Total length (outliers include lengths outside the usual range for the MLST species)",
                    "Mb",
                    {s: (run.quast[s].get("total_length") or 0) / 1e6 for s in quast_samples},
                    flagged=unusual,
                    digits=2,
                ),
                Panel(
                    "Contigs (>= 500 bp)", "contigs", {s: run.quast[s].get("contigs") for s in quast_samples}, digits=0
                ),
                Panel("N50", "kb", {s: (run.quast[s].get("n50") or 0) / 1e3 for s in quast_samples}),
            ]
        )
    elif quast_samples:
        ordered = list(reversed(quast_samples))
        fig = make_subplots(
            rows=1,
            cols=3,
            shared_yaxes=True,
            horizontal_spacing=0.03,
            subplot_titles=("Total length", "Contigs (>= 500 bp)", "N50"),
        )
        lengths = [(run.quast[s].get("total_length") or 0) / 1e6 for s in ordered]
        fig.add_trace(
            go.Bar(
                y=ordered,
                x=lengths,
                orientation="h",
                marker_color=PALETTE[0],
                name="Total length",
                showlegend=False,
                hovertemplate="%{y}<br>%{x:.3f} Mb<extra></extra>",
            ),
            row=1,
            col=1,
        )
        ranges = [(s, run.expected_length(s)) for s in ordered]
        if any(expected for _, expected in ranges):
            xs, ys = [], []
            for s, expected in ranges:
                if expected:
                    xs += [expected[0], expected[1], None]
                    ys += [s, s, None]
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="lines+markers",
                    line=dict(color=INK, width=2),
                    marker=dict(symbol="line-ns", size=10, line=dict(width=2, color=INK)),
                    name="Usual length for the MLST species",
                    hoverinfo="skip",
                ),
                row=1,
                col=1,
            )
        fig.add_trace(
            go.Bar(
                y=ordered,
                x=[run.quast[s].get("contigs") or 0 for s in ordered],
                orientation="h",
                marker_color=PALETTE[0],
                showlegend=False,
                hovertemplate="%{y}<br>%{x:,} contigs<extra></extra>",
            ),
            row=1,
            col=2,
        )
        fig.add_trace(
            go.Bar(
                y=ordered,
                x=[(run.quast[s].get("n50") or 0) / 1e3 for s in ordered],
                orientation="h",
                marker_color=PALETTE[0],
                showlegend=False,
                hovertemplate="%{y}<br>N50 %{x:,.1f} kb<extra></extra>",
            ),
            row=1,
            col=3,
        )
        fig.update_xaxes(title_text="Mb", row=1, col=1)
        fig.update_xaxes(title_text="Contigs", row=1, col=2)
        fig.update_xaxes(title_text="kb", row=1, col=3)
        fig.update_layout(bargap=0.3, margin=dict(l=10, r=20, t=90, b=50))
        fig.update_annotations(font=dict(size=13, color=INK))
        # the legend goes above the subplot titles
        body += figure_html(fig, bar_height(len(quast_samples), extra=180), legend_y=1.06)

    rows = []
    for s in samples:
        q = run.quast.get(s, {})
        f = run.filter_stats.get(s, {})
        expected = run.expected_length(s)
        rows.append(
            {
                "sample": sample_link(s),
                "contigs": fmt_int(q.get("contigs")),
                "length": fmt_mb(q.get("total_length"), 3),
                "expected": f"{expected[0]:.1f}-{expected[1]:.1f}" if expected else "",
                "largest": fmt_int(q.get("largest_contig")),
                "n50": fmt_int(q.get("n50")),
                "gc": fmt_float(q.get("gc"), 2),
                "kept_contigs": fmt_int(f.get("annotation_contigs")),
                "kept_length": fmt_mb(f.get("annotation_length"), 3),
            }
        )
    columns = [
        ("sample", "Sample"),
        ("contigs", "Contigs (>= 500 bp)"),
        ("length", "Total length (Mb)"),
        ("expected", "Usual length (Mb)"),
        ("largest", "Largest contig (bp)"),
        ("n50", "N50 (bp)"),
        ("gc", "GC (%)"),
        ("kept_contigs", "Contigs annotated"),
        ("kept_length", "Length annotated (Mb)"),
    ]
    if not run.filter_stats:
        columns = columns[:-2]
    if not any(run.expected_length(s) for s in samples):
        columns = [c for c in columns if c[0] != "expected"]
    body += render_table(columns, rows, numeric={c for c, _ in columns[1:]}, raw=("sample",))

    min_length = (run.run_info.get("parameters") or {}).get("min_contig_length")
    threshold = f"{to_number(min_length):,} bp" if to_number(min_length) else "--min_contig_length"
    intro = [
        "The cleaned reads are assembled into contigs (SPAdes for Illumina; Flye with Medaka, or Autocycler, "
        "for Nanopore). QUAST measures the assembly: total length should be close to the usual genome size of "
        "the species, and fewer contigs with a larger N50 (half of the assembly is in contigs at least this "
        "long) mean a less fragmented assembly. QUAST counts contigs of at least 500 bp.",
        f"Annotation and typing use only contigs of at least {esc(threshold)}; the last two columns give what "
        "was kept. The usual length range is approximate and is shown for the species given by the MLST "
        "scheme.",
    ]
    return Section("assembly", "Assembly", intro, body)


def section_quality(run):
    samples = [s for s in run.samples if s in run.checkm2]
    if not samples:
        return None
    groups = [
        (
            f">= {CHECKM2_MIN_COMPLETENESS:g} % complete and <= {CHECKM2_MAX_CONTAMINATION:g} % contaminated",
            PALETTE[0],
            "circle",
            [s for s in samples if run.checkm2_group(s) is None],
        )
    ]
    groups += [
        (name, colour, symbol, [s for s in samples if run.checkm2_group(s) == key])
        for key, name, colour, symbol in CHECKM2_GROUPS
    ]
    outside = sum(len(group) for _, _, _, group in groups[1:])
    fig = go.Figure()
    for slot, (name, colour, symbol, group) in enumerate(groups):
        if not group:
            continue
        fig.add_trace(
            go.Scatter(
                x=[run.checkm2[s]["completeness"] for s in group],
                y=[run.checkm2[s]["contamination"] for s in group],
                text=group,
                # passing genomes sit together near 100 / 0, so only the ones to look at are labelled
                mode="markers+text" if slot and outside <= 25 else "markers",
                textposition="top center",
                textfont=dict(size=10, color=INK_SECONDARY),
                marker=dict(
                    size=12 if slot else 10,
                    color=colour,
                    symbol=symbol,
                    line=dict(color=colour if symbol == "x" else SURFACE, width=1 if symbol == "x" else 2),
                ),
                name=name,
                hovertemplate="%{text}<br>Completeness %{x:.2f} %<br>Contamination %{y:.2f} %<extra></extra>",
            )
        )
    low_x = min(min(run.checkm2[s]["completeness"] for s in samples) - 2, 85)
    high_y = max(max(run.checkm2[s]["contamination"] for s in samples) + 1, 8)
    fig.add_vline(x=CHECKM2_MIN_COMPLETENESS, line_dash="dash", line_color=INK_SECONDARY)
    fig.add_hline(y=CHECKM2_MAX_CONTAMINATION, line_dash="dash", line_color=INK_SECONDARY)
    fig.update_xaxes(title_text="Completeness (%)", range=[low_x, 100.5])
    fig.update_yaxes(title_text="Contamination (%)", range=[-0.3, high_y])
    fig.update_layout(margin=dict(l=10, r=20, t=60, b=50), showlegend=True)
    body = figure_html(fig, 460)
    rows = [
        {
            "sample": sample_link(s),
            "completeness": fmt_float(run.checkm2[s]["completeness"], 2),
            "contamination": fmt_float(run.checkm2[s]["contamination"], 2),
            "model": run.checkm2[s].get("model") or "",
            "check": (
                f'<span class="flag">{esc(dict((k, n) for k, n, _, _ in CHECKM2_GROUPS)[run.checkm2_group(s)].lower())}</span>'
                if run.checkm2_group(s)
                else "within limits"
            ),
        }
        # genomes outside the limits first, then name order
        for s in sorted(samples, key=lambda s: run.checkm2_group(s) is None)
    ]
    body += render_table(
        [
            ("sample", "Sample"),
            ("completeness", "Completeness (%)"),
            ("contamination", "Contamination (%)"),
            ("model", "CheckM2 model"),
            ("check", "Check"),
        ],
        rows,
        numeric=("completeness", "contamination"),
        raw=("sample", "check"),
    )
    intro = [
        "CheckM2 estimates how complete the genome is and how much of it looks like it comes from another "
        "organism, from the genes found in the assembly. Points in the lower right corner, at least "
        f"{CHECKM2_MIN_COMPLETENESS:g} % complete and at most {CHECKM2_MAX_CONTAMINATION:g} % contaminated "
        "(the dashed lines; the usual high-quality draft limits), are good genomes. Low completeness points to "
        "a poor assembly; high contamination to a mixed culture. Genomes outside the limits have their own colour "
        "and marker shape (low completeness, high contamination, or both) and are labelled; hover a point for "
        "its sample. The table lists every value, genomes outside the limits first.",
    ]
    return Section("quality", "Genome quality (CheckM2)", intro, body)


def section_typing(run):
    samples = [s for s in run.samples if s in run.mlst or s in run.gtdbtk]
    if not samples and not (run.kleborate or run.sistr or run.sccmec):
        return None
    body = ""
    rows = []
    for s in samples:
        m = run.mlst.get(s, {})
        g = run.gtdbtk.get(s, {})
        taxa = run.taxa.get(s)
        disagreement = run.species_disagreement(s)
        if disagreement:
            check = '<span class="flag">differs</span>'
        elif not run.organism(s):
            check = ""
        else:
            check = "agrees" if run.species_checked(s) else "no second call"
        rows.append(
            {
                "_flagged": bool(disagreement),
                "sample": sample_link(s),
                "scheme": m.get("scheme") or "none",
                "st": m.get("st") or "",
                "mlst_species": f"<i>{esc(species_name(run.organism(s)))}</i>" if run.organism(s) else "",
                "gtdb": gtdb_name(g),
                "ani": fmt_float(g.get("ani"), 2),
                "reads": f"<i>{esc(taxa.iloc[0]['name'])}</i>" if taxa is not None and len(taxa) else "",
                "check": check,
            }
        )
    # disagreements first, then name order (sorted() keeps the order within each group)
    rows = sorted(rows, key=lambda row: not row.pop("_flagged"))
    if rows:
        columns = [
            ("sample", "Sample"),
            ("scheme", "MLST scheme"),
            ("st", "ST"),
            ("mlst_species", "Species from MLST"),
            ("gtdb", "GTDB-Tk species"),
            ("ani", "GTDB-Tk ANI (%)"),
            ("reads", "Top species in reads"),
            ("check", "Genus check"),
        ]
        if not run.gtdbtk:
            columns = [c for c in columns if c[0] not in ("gtdb", "ani")]
        if not run.taxa:
            columns = [c for c in columns if c[0] != "reads"]
        if not run.gtdbtk and not run.taxa:
            columns = [c for c in columns if c[0] != "check"]
            body += (
                '<p class="note">Neither GTDB-Tk nor a read classification ran, so the MLST species could not be '
                "checked against a second call.</p>"
            )
        body += render_table(columns, rows, numeric=("ani",), raw=("sample", "mlst_species", "gtdb", "reads", "check"))

    if run.kleborate:
        rows = []
        for s in [s for s in run.samples if s in run.kleborate]:
            k = run.kleborate[s]
            rows.append(
                {
                    "sample": sample_link(s),
                    "species": f"<i>{esc(k.get('species', ''))}</i>",
                    "st": k.get("ST", ""),
                    "virulence": k.get("virulence_score", ""),
                    "resistance": k.get("resistance_score", ""),
                    "k": k.get("K_locus", ""),
                    "o": k.get("O_type", k.get("O_locus", "")),
                    "ybt": k.get("Yersiniabactin", ""),
                    "iuc": k.get("Aerobactin", ""),
                    "carbapenemase": k.get("Bla_Carb_acquired", ""),
                    "esbl": k.get("Bla_ESBL_acquired", ""),
                }
            )
        body += "<h3>Kleborate (<i>Klebsiella pneumoniae</i> species complex)</h3>" + render_table(
            [
                ("sample", "Sample"),
                ("species", "Species"),
                ("st", "ST"),
                ("virulence", "Virulence score (0-5)"),
                ("resistance", "Resistance score (0-3)"),
                ("k", "K locus"),
                ("o", "O type"),
                ("ybt", "Yersiniabactin"),
                ("iuc", "Aerobactin"),
                ("carbapenemase", "Carbapenemases"),
                ("esbl", "ESBLs"),
            ],
            rows,
            numeric=("virulence", "resistance"),
            raw=("sample", "species"),
        )
    if run.sistr:
        rows = []
        for s in [s for s in run.samples if s in run.sistr]:
            r = run.sistr[s]
            rows.append(
                {
                    "sample": sample_link(s),
                    "serovar": r.get("serovar", ""),
                    "serogroup": r.get("serogroup", ""),
                    "antigens": f"{r.get('o_antigen', '')}:{r.get('h1', '')}:{r.get('h2', '')}",
                    "cgmlst": r.get("serovar_cgmlst", ""),
                    "subspecies": r.get("cgmlst_subspecies", ""),
                    "qc": r.get("qc_status", ""),
                }
            )
        body += "<h3>SISTR (<i>Salmonella</i>)</h3>" + render_table(
            [
                ("sample", "Sample"),
                ("serovar", "Serovar"),
                ("serogroup", "Serogroup"),
                ("antigens", "O:H1:H2"),
                ("cgmlst", "Serovar (cgMLST)"),
                ("subspecies", "Subspecies"),
                ("qc", "SISTR QC"),
            ],
            rows,
            raw=("sample",),
        )
    if run.sccmec:
        rows = []
        # calls that are not clear-cut first, then name order
        for s in sorted([s for s in run.samples if s in run.sccmec], key=lambda s: not sccmec_problem(run.sccmec[s])):
            r = run.sccmec[s]
            problem = sccmec_problem(r)
            rows.append(
                {
                    "sample": sample_link(s),
                    "type": r.get("type", ""),
                    "subtype": r.get("subtype", ""),
                    "mecA": r.get("mecA", ""),
                    "call": f'<span class="flag">{esc(problem)}</span>' if problem else "clear",
                    "comment": r.get("comment", ""),
                }
            )
        body += "<h3>SCCmec (<i>Staphylococcus aureus</i>)</h3>" + render_table(
            [
                ("sample", "Sample"),
                ("type", "SCCmec type"),
                ("subtype", "Subtype"),
                ("mecA", "mecA"),
                ("call", "Call"),
                ("comment", "Comment"),
            ],
            rows,
            raw=("sample", "call"),
        )
    intro = [
        "MLST compares seven housekeeping genes with a scheme for the species and gives a sequence type (ST); "
        "the scheme it picks also names the species and decides which species tools run (Kleborate for "
        "<i>Klebsiella pneumoniae</i>, SISTR for <i>Salmonella</i>, sccmec for <i>Staphylococcus aureus</i>). "
        "An empty ST means a new or incomplete allele combination.",
        "GTDB-Tk (when it ran) and the read classification (Bracken, or Kraken2 without it) are independent "
        "species calls: when either genus differs from the MLST species the row is flagged and listed first, "
        "and the species tools may have been chosen wrongly.",
    ]
    if run.sccmec:
        intro.append(
            "sccmec types the cassette from its target genes (mec and ccr complexes); when they do not give one "
            "type it falls back to the closest whole cassette. Calls without a clear type from the target genes "
            "are flagged and listed first."
        )
    return Section("typing", "Species and typing", intro, body)


def column(table, name):
    """A table column, or an empty-string column when the tool version does not write it."""
    return table[name] if name in table else pd.Series("", index=table.index)


def shorten(text, limit=HEATMAP_LABEL):
    return text if len(text) <= limit else text[: limit - 1] + "…"


def amr_heatmap(tables, samples, group_column):
    """Gene-by-sample heatmap; cell colour = % identity of the best hit, rows grouped by class."""
    records = []
    for s in samples:
        table = tables[s]
        for _, row in table.iterrows():
            records.append(
                {
                    "sample": s,
                    "gene": row["gene"],
                    "group": row.get(group_column) or "other",
                    "identity": row.get("identity"),
                    "coverage": row.get("coverage"),
                }
            )
    if not records:
        return ""
    frame = pd.DataFrame(records)
    best = frame.sort_values("identity", ascending=False).drop_duplicates(["sample", "gene"])
    groups = best.groupby("gene")["group"].agg(lambda values: sorted(set(values))[0])
    genes = sorted(groups.index, key=lambda gene: (groups[gene].lower(), gene.lower()))
    z, hover = [], []
    lookup = {(r.sample, r.gene): r for r in best.itertuples()}
    for gene in genes:
        z_row, hover_row = [], []
        for s in samples:
            hit = lookup.get((s, gene))
            if hit is None:
                # no hover text on empty cells (hoverongaps=False): it doubles the page size on large runs
                z_row.append(None)
                hover_row.append("")
            else:
                z_row.append(hit.identity)
                hover_row.append(
                    f"{s}<br>{gene} ({esc(groups[gene])})<br>identity {fmt_float(hit.identity, 2)} %"
                    f"<br>coverage {fmt_float(hit.coverage, 2)} %"
                )
        z.append(z_row)
        hover.append(hover_row)
    # rows are keyed by the full name, so two genes whose shortened labels match stay two rows
    keys = [f"{groups[gene].lower()} · {gene}" for gene in genes]
    labels = [shorten(key) for key in keys]
    scale = [[index / (len(SEQUENTIAL) - 1), colour] for index, colour in enumerate(SEQUENTIAL)]
    fig = go.Figure(
        go.Heatmap(
            z=z,
            x=samples,
            y=keys,
            text=hover,
            hovertemplate="%{text}<extra></extra>",
            hoverongaps=False,
            colorscale=scale,
            zmin=min(80.0, best["identity"].min()) if best["identity"].notna().any() else 80.0,
            zmax=100,
            xgap=2,
            ygap=2,
            colorbar=dict(
                title=dict(text="Identity (%)", side="right"), thickness=12, len=min(1.0, 300 / (len(genes) * 18 + 1))
            ),
        )
    )
    # fixed column width; the left margin follows the longest gene label, and sample names turn when
    # they are wider than a column
    left = int(6.5 * max(len(label) for label in labels)) + 20
    turned = 7 * max(len(s) for s in samples) > HEATMAP_COLUMN
    top = 60 + (int(5 * max(len(s) for s in samples)) if turned else 20)
    fig.update_yaxes(
        autorange="reversed", tickfont=dict(size=11), showgrid=False, tickmode="array", tickvals=keys, ticktext=labels
    )
    fig.update_xaxes(side="top", tickangle=-60 if turned else 0, showgrid=False, title_text="")
    width = left + HEATMAP_COLUMN * len(samples) + 110
    fig.update_layout(width=width, margin=dict(l=left, r=10, t=top, b=10))
    return '<div class="scroll-x">' + figure_html(fig, max(300, 18 * len(genes) + top + 30)) + "</div>"


def section_amr(run):
    amr_samples = [s for s in run.samples if s in run.amrfinder]
    rgi_samples = [s for s in run.samples if s in run.rgi]
    if not amr_samples and not rgi_samples:
        return None
    body = ""
    if amr_samples:
        genes = {}
        for s in amr_samples:
            table = run.amrfinder[s]
            genes[s] = table[(column(table, "type") == "AMR") & (column(table, "subtype") != "POINT")]
        body += "<h3>AMRFinderPlus: acquired resistance genes</h3>"
        heat = amr_heatmap(genes, amr_samples, "class")
        body += heat or "<p>No acquired resistance genes were found.</p>"

        points = []
        for s in amr_samples:
            table = run.amrfinder[s]
            for _, row in table[column(table, "subtype") == "POINT"].iterrows():
                points.append(
                    {
                        "sample": sample_link(s),
                        "mutation": row["gene"],
                        "class": row["class"],
                        "subclass": row["subclass"],
                        "name": row["name"],
                    }
                )
        if points:
            body += "<h3>Point mutations</h3>" + render_table(
                [
                    ("sample", "Sample"),
                    ("mutation", "Mutation"),
                    ("class", "Class"),
                    ("subclass", "Subclass"),
                    ("name", "Description"),
                ],
                points,
                raw=("sample",),
            )
        hits = []
        for s in amr_samples:
            for _, row in run.amrfinder[s].iterrows():
                hits.append(
                    {
                        "sample": sample_link(s),
                        "gene": row["gene"],
                        "type": row.get("type", ""),
                        "class": row.get("class", ""),
                        "subclass": row.get("subclass", ""),
                        "method": row.get("method", ""),
                        "identity": fmt_float(row.get("identity"), 2),
                        "coverage": fmt_float(row.get("coverage"), 2),
                        "contig": row.get("contig", ""),
                    }
                )
        body += "<h3>All AMRFinderPlus hits</h3>" + render_table(
            [
                ("sample", "Sample"),
                ("gene", "Gene"),
                ("type", "Type"),
                ("class", "Class"),
                ("subclass", "Subclass"),
                ("method", "Method"),
                ("identity", "Identity (%)"),
                ("coverage", "Coverage (%)"),
                ("contig", "Contig"),
            ],
            hits,
            numeric=("identity", "coverage"),
            raw=("sample",),
            max_height=480,
        )
    if rgi_samples:
        strict = {}
        for s in rgi_samples:
            table = run.rgi[s]
            table = table[table["cut_off"].isin(["Perfect", "Strict"])].copy()
            table["drug_group"] = table["drug_class"].str.split(";").str[0].str.strip()
            strict[s] = table
        body += "<h3>RGI (CARD): Perfect and Strict hits</h3>"
        heat = amr_heatmap(strict, rgi_samples, "drug_group")
        body += heat or "<p>No Perfect or Strict RGI hits were found.</p>"
    intro = [
        "AMRFinderPlus (NCBI) finds acquired resistance genes and, for supported species, resistance point "
        "mutations; the species comes from the MLST scheme. In the heatmap each row is a gene, grouped by drug "
        "class, and each column a sample; a coloured cell means the gene was found, darker blue meaning closer "
        "to the reference sequence. A gene found does not by itself prove phenotypic resistance.",
    ]
    if rgi_samples:
        intro.append(
            "RGI searches the CARD database and adds efflux pumps and regulatory genes that AMRFinderPlus "
            "reports less often; only Perfect and Strict hits are shown, grouped by the first drug class CARD "
            "lists for the gene."
        )
    intro.append(
        "Layout: one column per sample (34 px each; with many samples the chart scrolls sideways) and one row per "
        "gene; hover a cell for identity and coverage."
    )
    return Section("amr", "Antimicrobial resistance", intro, body)


def section_mobile(run):
    genomad_samples = [s for s in run.samples if s in run.genomad_plasmid or s in run.genomad_virus]
    mob_samples = [s for s in run.samples if s in run.mobsuite_contigs or s in run.mobsuite_mobtyper]
    if not genomad_samples and not mob_samples:
        return None
    body = ""
    if len(genomad_samples) > LARGE_RUN:
        body += "<h3>geNomad</h3>" + distribution_html(
            [
                Panel(
                    "Plasmid contigs",
                    "contigs",
                    {s: len(run.genomad_plasmid.get(s, [])) for s in genomad_samples},
                    digits=0,
                ),
                Panel(
                    "Virus sequences",
                    "sequences",
                    {s: len(run.genomad_virus.get(s, [])) for s in genomad_samples},
                    digits=0,
                ),
            ]
        )
    if genomad_samples and len(genomad_samples) <= LARGE_RUN:
        ordered = list(reversed(genomad_samples))
        fig = go.Figure()
        for slot, (label, source) in enumerate(
            (("Plasmids", run.genomad_plasmid), ("Viruses and proviruses", run.genomad_virus))
        ):
            fig.add_trace(
                go.Bar(
                    y=ordered,
                    x=[len(source.get(s, [])) for s in ordered],
                    orientation="h",
                    name=label,
                    marker_color=PALETTE[slot],
                    hovertemplate=f"%{{y}}<br>{label}: %{{x}}<extra></extra>",
                )
            )
        fig.update_layout(barmode="group", bargap=0.25, margin=dict(l=10, r=20, t=40, b=50))
        most = max(len(source.get(s, [])) for source in (run.genomad_plasmid, run.genomad_virus) for s in ordered)
        # whole-number ticks for small counts; Plotly picks the step for larger ones
        fig.update_xaxes(title_text="Sequences called by geNomad", dtick=1 if most <= 10 else None)
        body += "<h3>geNomad</h3>" + figure_html(fig, bar_height(len(genomad_samples), per_row=26))
    if genomad_samples:
        rows = []
        for s in genomad_samples:
            plasmids = run.genomad_plasmid.get(s, pd.DataFrame())
            viruses = run.genomad_virus.get(s, pd.DataFrame())
            plasmid_bp = pd.to_numeric(plasmids.get("length", pd.Series(dtype=float)), errors="coerce").sum()
            circular = (
                int((plasmids.get("topology", pd.Series(dtype=str)) != "No terminal repeats").sum())
                if len(plasmids)
                else 0
            )
            amr = (
                int((plasmids.get("amr_genes", pd.Series(dtype=str)).replace("NA", "") != "").sum())
                if len(plasmids)
                else 0
            )
            rows.append(
                {
                    "sample": sample_link(s),
                    "plasmids": len(plasmids),
                    "plasmid_kb": fmt_float(plasmid_bp / 1e3) if len(plasmids) else "",
                    "circular": circular if len(plasmids) else "",
                    "amr": amr if len(plasmids) else "",
                    "viruses": len(viruses),
                    "proviruses": (
                        int((viruses.get("topology", pd.Series(dtype=str)) == "Provirus").sum()) if len(viruses) else ""
                    ),
                }
            )
        body += render_table(
            [
                ("sample", "Sample"),
                ("plasmids", "Plasmid contigs"),
                ("plasmid_kb", "Plasmid length (kb)"),
                ("circular", "With terminal repeats"),
                ("amr", "Plasmid contigs with AMR genes"),
                ("viruses", "Virus sequences"),
                ("proviruses", "Proviruses"),
            ],
            rows,
            numeric=("plasmids", "plasmid_kb", "circular", "amr", "viruses", "proviruses"),
            raw=("sample",),
        )
    if mob_samples:
        rows = []
        for s in mob_samples:
            typer = run.mobsuite_mobtyper.get(s, pd.DataFrame())
            contigs = run.mobsuite_contigs.get(s, pd.DataFrame())
            if len(typer):
                for _, row in typer.iterrows():
                    rows.append(
                        {
                            "sample": sample_link(s),
                            "plasmid": row.get("sample_id", "").split(":")[-1],
                            "cluster": row.get("primary_cluster_id", ""),
                            "size": fmt_float(to_number(row.get("size")) / 1e3) if to_number(row.get("size")) else "",
                            "contigs": row.get("num_contigs", ""),
                            "replicons": row.get("rep_type(s)", "").replace(",", ", "),
                            "relaxase": row.get("relaxase_type(s)", "").replace(",", ", "),
                            "mpf": row.get("mpf_type", ""),
                            "mobility": row.get("predicted_mobility", ""),
                            "host": row.get(
                                "observed_host_range_ncbi_name", row.get("predicted_host_range_overall_name", "")
                            ),
                        }
                    )
            elif len(contigs):
                rows.append({"sample": sample_link(s), "plasmid": "no plasmid reconstructed"})
        body += "<h3>MOB-suite</h3>" + render_table(
            [
                ("sample", "Sample"),
                ("plasmid", "Plasmid"),
                ("cluster", "Cluster"),
                ("size", "Size (kb)"),
                ("contigs", "Contigs"),
                ("replicons", "Replicons"),
                ("relaxase", "Relaxase"),
                ("mpf", "Mating pair formation"),
                ("mobility", "Predicted mobility"),
                ("host", "Host range"),
            ],
            rows,
            numeric=("size", "contigs"),
            raw=("sample",),
            max_height=480,
        )
    intro = [
        "Mobile genetic elements move resistance and virulence genes between bacteria. geNomad classifies "
        "each contig as chromosome, plasmid or virus (including prophages inside the chromosome); the chart "
        "counts the plasmid and virus sequences per sample. In a short-read assembly a plasmid is often split "
        "into several contigs, so counts are contigs, not plasmids.",
    ]
    if mob_samples:
        intro.append(
            "MOB-suite groups contigs into reconstructed plasmids, names their replicon and relaxase types and "
            "predicts mobility: conjugative plasmids can transfer themselves, mobilizable ones need the help of "
            "another plasmid, non-mobilizable ones cannot move."
        )
    return Section("mobile", "Plasmids and other mobile elements", intro, body)


# consensus calls that point to a plasmid, in display order (bin/plasmid_consensus.py)
PLASMID_CALL_LABELS = OrderedDict((("plasmid", "plasmid"), ("candidate", "candidate"), ("check", "to check")))


def plasmid_calls(table):
    """{call: contigs} for the plasmid calls of one sample's consensus table."""
    if table is None or "consensus" not in table:
        return {}
    counts = table["consensus"].value_counts()
    return {call: int(counts[call]) for call in PLASMID_CALL_LABELS if call in counts}


def unique_items(values):
    items = []
    for value in values:
        for item in str(value).split(","):
            item = item.strip()
            if item and item not in ("-", "not run") and item not in items:
                items.append(item)
    return items


def plasmid_replicons(table):
    """Replicons (PlasmidFinder, then MOB-suite) on the contigs called plasmid."""
    if table is None or "consensus" not in table:
        return ""
    plasmids = table[table["consensus"] == "plasmid"]
    return ", ".join(
        unique_items(list(column(plasmids, "plasmidfinder_replicons")) + list(column(plasmids, "mobsuite_replicons")))
    )


def section_plasmids(run):
    samples = [s for s in run.samples if s in run.plasmid_consensus]
    if not samples:
        return None
    body = ""
    counts = {s: plasmid_calls(run.plasmid_consensus[s]) for s in samples}
    if len(samples) > LARGE_RUN:
        body += distribution_html(
            [
                Panel(f"Contigs called {label}", "contigs", {s: counts[s].get(call, 0) for s in samples}, digits=0)
                for call, label in PLASMID_CALL_LABELS.items()
            ]
        )
    else:
        ordered = list(reversed(samples))
        fig = go.Figure()
        for slot, (call, label) in enumerate(PLASMID_CALL_LABELS.items()):
            fig.add_trace(
                go.Bar(
                    y=ordered,
                    x=[counts[s].get(call, 0) for s in ordered],
                    orientation="h",
                    name=label,
                    marker_color=PALETTE[slot],
                    hovertemplate=f"%{{y}}<br>{label}: %{{x}} contigs<extra></extra>",
                )
            )
        # legend in call order (Plotly lists stacked horizontal traces last-first)
        fig.update_layout(
            barmode="stack", bargap=0.25, margin=dict(l=10, r=20, t=40, b=50), legend=dict(traceorder="normal")
        )
        most = max((sum(c.values()) for c in counts.values()), default=0)
        fig.update_xaxes(title_text="Contigs by consensus call", dtick=1 if most <= 10 else None)
        body += figure_html(fig, bar_height(len(samples), per_row=26))

    summary, contigs = [], []
    for s in samples:
        table = run.plasmid_consensus[s]
        plasmids = table[table["consensus"] == "plasmid"] if "consensus" in table else table.iloc[0:0]
        length = pd.to_numeric(column(plasmids, "length"), errors="coerce").sum()
        mobility = pd.Series(column(plasmids, "predicted_mobility")).replace({"-": "", "not run": ""})
        mobility = mobility[mobility != ""].value_counts()
        summary.append(
            {
                "sample": sample_link(s),
                "plasmid": counts[s].get("plasmid", 0),
                "plasmid_kb": fmt_float(length / 1e3) if len(plasmids) else "",
                "circular": int((pd.Series(column(plasmids, "circular")) == "yes").sum()) if len(plasmids) else "",
                "candidate": counts[s].get("candidate", 0),
                "check": counts[s].get("check", 0),
                "replicons": plasmid_replicons(table),
                "mobility": ", ".join(f"{n} {name}" for name, n in mobility.items()),
            }
        )
        if "consensus" not in table:
            continue
        for _, row in table[table["consensus"].isin(list(PLASMID_CALL_LABELS))].iterrows():
            contigs.append(
                {
                    "sample": sample_link(s),
                    "contig": row.get("contig", ""),
                    "length": fmt_float(to_number(row.get("length")) / 1e3) if to_number(row.get("length")) else "",
                    "circular": row.get("circular", ""),
                    "call": PLASMID_CALL_LABELS.get(row.get("consensus", ""), row.get("consensus", "")),
                    "evidence": row.get("evidence", ""),
                    # one name per replicon or relaxase type; the consensus TSV keeps every hit
                    "plasmidfinder": ", ".join(unique_items([row.get("plasmidfinder_replicons", "")])) or "-",
                    "mobsuite": ", ".join(unique_items([row.get("mobsuite_replicons", "")])) or "-",
                    "relaxase": ", ".join(unique_items([row.get("mobsuite_relaxases", "")])) or "-",
                    "mobility": row.get("predicted_mobility", ""),
                    "cluster": row.get("mobsuite_cluster", ""),
                }
            )
    body += render_table(
        [
            ("sample", "Sample"),
            ("plasmid", "Plasmid contigs"),
            ("plasmid_kb", "Plasmid length (kb)"),
            ("circular", "Circular plasmid contigs"),
            ("candidate", "Candidate contigs"),
            ("check", "Contigs to check"),
            ("replicons", "Replicons on plasmid contigs"),
            ("mobility", "Mobility of plasmid contigs"),
        ],
        summary,
        numeric=("plasmid", "plasmid_kb", "circular", "candidate", "check"),
        raw=("sample",),
    )
    if contigs:
        body += "<h3>Contigs with plasmid evidence</h3>" + render_table(
            [
                ("sample", "Sample"),
                ("contig", "Contig"),
                ("length", "Length (kb)"),
                ("circular", "Circular"),
                ("call", "Call"),
                ("evidence", "Evidence"),
                ("plasmidfinder", "PlasmidFinder replicons"),
                ("mobsuite", "MOB-suite replicons"),
                ("relaxase", "Relaxase"),
                ("mobility", "Predicted mobility"),
                ("cluster", "MOB-suite cluster"),
            ],
            contigs,
            numeric=("length",),
            raw=("sample",),
            max_height=480,
        )
    intro = [
        "Each contig is checked by up to three tools that find plasmids in different ways: geNomad (genes "
        "typical of plasmids), MOB-suite (similarity to known plasmids, replicon and relaxase genes; on long-read "
        "and finished assemblies MOB-typer types each circular contig) and PlasmidFinder (replicon sequences). "
        "A contig is called plasmid when at least two of them agree, candidate when only one does, and "
        "to check when the evidence sits on a chromosome-sized contig (1 Mb or more).",
        "A short-read assembly splits a plasmid into several contigs, so the counts are contigs, not plasmids; "
        "the MOB-suite cluster groups the contigs of one plasmid. Circular contigs from a long-read assembly are "
        "usually whole plasmids. Mobility: conjugative plasmids can transfer themselves, mobilizable ones need "
        "another plasmid's help, non-mobilizable ones cannot move. The rules are in docs/output.md.",
    ]
    return Section("plasmids", "Plasmid consensus", intro, body)


# IS families named in the chart; the rest are summed as "other families" (PALETTE has 8 colours)
TOP_IS_FAMILIES = 7


def is_families(table):
    """{family: elements} for one sample's ISEScan table, most elements first."""
    if table is None or "family" not in table or not len(table):
        return {}
    counts = table["family"].replace("", "unknown").value_counts()
    return {family: int(n) for family, n in sorted(counts.items(), key=lambda item: (-item[1], natural_key(item[0])))}


def is_complete(table):
    """Complete IS elements (ISEScan type c); the others are partial (type p)."""
    return int((column(table, "type") == "c").sum()) if table is not None and len(table) else 0


def is_family_text(families, limit=None):
    """ "IS3 12, IS5 6, ..." with the most elements first; limit cuts the list to that many families."""
    named = list(families.items())[:limit] if limit else list(families.items())
    text = ", ".join(f"{family} {n}" for family, n in named)
    return text + (", ..." if limit and len(families) > limit else "")


def is_on_plasmids(run, sample, call="plasmid"):
    """IS elements on contigs the plasmid consensus gives this call, or None without a consensus table."""
    table = run.isescan.get(sample)
    consensus = run.plasmid_consensus.get(sample)
    if table is None or consensus is None or not {"consensus", "contig"} <= set(consensus) or "seqID" not in table:
        return None
    plasmids = set(consensus.loc[consensus["consensus"] == call, "contig"])
    return int(table["seqID"].isin(plasmids).sum())


def section_isescan(run):
    samples = [s for s in run.samples if s in run.isescan]
    if not samples:
        return None
    families = {s: is_families(run.isescan[s]) for s in samples}
    body = ""
    if len(samples) > LARGE_RUN:
        body += distribution_html(
            [
                Panel("IS elements", "elements", {s: len(run.isescan[s]) for s in samples}, digits=0),
                Panel("Complete IS elements", "elements", {s: is_complete(run.isescan[s]) for s in samples}, digits=0),
                Panel("IS families", "families", {s: len(families[s]) for s in samples}, digits=0),
            ]
        )
    elif any(families.values()):
        totals = defaultdict(int)
        for counts in families.values():
            for family, n in counts.items():
                totals[family] += n
        ranked = sorted(totals, key=lambda family: (-totals[family], natural_key(family)))
        top = ranked[:TOP_IS_FAMILIES] if len(ranked) > TOP_IS_FAMILIES + 1 else ranked
        ordered = list(reversed(samples))
        fig = go.Figure()
        for index, family in enumerate(top):
            fig.add_trace(
                go.Bar(
                    y=ordered,
                    x=[families[s].get(family, 0) for s in ordered],
                    orientation="h",
                    name=family,
                    marker=dict(color=PALETTE[index], line=dict(color=SURFACE, width=1)),
                    hovertemplate=f"%{{y}}<br>{esc(family)}: %{{x}} elements<extra></extra>",
                )
            )
        if len(top) < len(ranked):
            other = [sum(n for family, n in families[s].items() if family not in top) for s in ordered]
            fig.add_trace(
                go.Bar(
                    y=ordered,
                    x=other,
                    orientation="h",
                    name="Other families",
                    marker=dict(color=OTHER_COLOUR, line=dict(color=SURFACE, width=1)),
                    hovertemplate="%{y}<br>Other families: %{x} elements<extra></extra>",
                )
            )
        # legend in family order (Plotly lists stacked horizontal traces last-first)
        fig.update_layout(
            barmode="stack", bargap=0.25, margin=dict(l=10, r=20, t=60, b=50), legend=dict(traceorder="normal")
        )
        most = max(len(run.isescan[s]) for s in samples)
        fig.update_xaxes(title_text="IS elements by family", dtick=1 if most <= 10 else None)
        body += figure_html(fig, bar_height(len(samples), per_row=26, extra=150))

    on_plasmids = {s: is_on_plasmids(run, s) for s in samples}
    on_candidates = {s: is_on_plasmids(run, s, "candidate") for s in samples}
    rows = []
    for s in samples:
        table = run.isescan[s]
        length = run.filter_stats.get(s, {}).get("annotation_length")
        is_bp = float(table["isLen"].sum()) if "isLen" in table else 0.0
        complete = is_complete(table)
        rows.append(
            {
                "sample": sample_link(s),
                "elements": len(table),
                "complete": complete,
                "partial": len(table) - complete,
                "families": len(families[s]),
                "kb": fmt_float(is_bp / 1e3),
                "percent": fmt_float(100 * is_bp / length, 2) if length else "",
                "plasmids": "" if on_plasmids[s] is None else on_plasmids[s],
                "candidates": "" if on_candidates[s] is None else on_candidates[s],
                "by_family": is_family_text(families[s]) or "none found",
            }
        )
    columns = [
        ("sample", "Sample"),
        ("elements", "IS elements"),
        ("complete", "Complete"),
        ("partial", "Partial"),
        ("families", "Families"),
        ("kb", "IS length (kb)"),
        ("percent", "Share of annotated contigs (%)"),
        ("plasmids", "On plasmid contigs"),
        ("candidates", "On candidate plasmid contigs"),
        ("by_family", "Elements by family"),
    ]
    if all(value is None for value in on_plasmids.values()):
        columns.remove(("plasmids", "On plasmid contigs"))
    # candidates (one plasmid tool only, e.g. PlasmidFinder without geNomad and MOB-suite) only when they carry any
    if not any(on_candidates.values()):
        columns.remove(("candidates", "On candidate plasmid contigs"))
    if not any(row["percent"] for row in rows):
        columns.remove(("percent", "Share of annotated contigs (%)"))
    body += render_table(
        columns,
        rows,
        numeric=("elements", "complete", "partial", "families", "kb", "percent", "plasmids", "candidates"),
        raw=("sample",),
    )
    intro = [
        "Insertion sequences (IS elements) are short mobile elements, usually 0.7 to 2.5 kb, that carry little "
        "more than the transposase gene that moves them. They copy or move themselves within and between the "
        "chromosome and plasmids, can switch genes off or on where they land, and pairs of them can move "
        "resistance genes together (IS26 of the IS6 family, for example, around many plasmid resistance genes). "
        "ISEScan finds them from their transposase genes and terminal inverted repeats and names the IS family "
        "as in ISfinder; ISNCY is ISfinder's group of elements not yet assigned to a family, and new marks an "
        "element whose transposase matches one of ISEScan's own groups that has no ISfinder family.",
        "Partial elements are shorter than the smallest known member of their family, or single copies with a "
        "weak transposase match or without clear inverted repeats at their ends; in a short-read assembly they "
        "are often cut at the end of a contig. Repeated IS copies break short-read assemblies, so an Illumina assembly can "
        "show fewer complete elements than a long-read one of the same strain. The XLSX workbook lists every "
        "element with its contig and position.",
    ]
    return Section("isescan", "Insertion sequences", intro, body)


def section_cards(run):
    if not run.samples:
        return None
    cards, chips = [], []
    for s in run.samples:
        status = run.status_row(s).get("status") or ("pass" if s in run.quast or s in run.mlst else "")
        facts = []
        reads = run.reads.get(s)
        if reads:
            facts.append(
                ("Bases kept", f"{fmt_mb(reads.get('bases_after'), 1)} Mb of {fmt_mb(reads.get('bases_before'), 1)} Mb")
            )
        taxa = run.taxa.get(s)
        if taxa is not None and len(taxa):
            facts.append(
                ("Top species in reads", f"<i>{esc(taxa.iloc[0]['name'])}</i> {taxa.iloc[0]['percent']:.1f} %")
            )
        q = run.quast.get(s)
        if q:
            facts.append(
                (
                    "Assembly",
                    f"{fmt_mb(q.get('total_length'))} Mb, {fmt_int(q.get('contigs'))} contigs, N50 {fmt_int(q.get('n50'))} bp",
                )
            )
        c = run.checkm2.get(s)
        if c:
            facts.append(("CheckM2", f"{c['completeness']:.1f} % complete, {c['contamination']:.2f} % contamination"))
        m = run.mlst.get(s)
        if m:
            species = f" (<i>{esc(species_name(run.organism(s)))}</i>)" if run.organism(s) else ""
            if m.get("scheme"):
                facts.append(("MLST", f"{esc(m['scheme'])} ST {esc(m.get('st') or 'not assigned')}{species}"))
            else:
                facts.append(("MLST", "no scheme matched"))
        g = run.gtdbtk.get(s)
        if g:
            facts.append(("GTDB-Tk", gtdb_name(g)))
        k = run.kleborate.get(s)
        if k:
            facts.append(
                (
                    "Kleborate",
                    f"virulence {esc(k.get('virulence_score', ''))}, resistance {esc(k.get('resistance_score', ''))}, {esc(k.get('K_locus', ''))}",
                )
            )
        r = run.sistr.get(s)
        if r:
            facts.append(("SISTR serovar", esc(r.get("serovar", ""))))
        r = run.sccmec.get(s)
        if r:
            facts.append(("SCCmec", esc(r.get("type") or "none") + (" (not clear-cut)" if sccmec_problem(r) else "")))
        if s in run.amrfinder:
            genes = run.amr_genes(s)
            facts.append(("AMR genes", f"{len(genes)}: " + esc(", ".join(genes)) if genes else "none"))
        if s in run.genomad_plasmid or s in run.genomad_virus:
            facts.append(
                (
                    "geNomad",
                    f"{len(run.genomad_plasmid.get(s, []))} plasmid, {len(run.genomad_virus.get(s, []))} virus sequences",
                )
            )
        if s in run.mobsuite_mobtyper:
            facts.append(("MOB-suite plasmids", str(len(run.mobsuite_mobtyper[s]))))
        if s in run.plasmid_consensus:
            calls = plasmid_calls(run.plasmid_consensus[s])
            facts.append(
                (
                    "Plasmid consensus",
                    ", ".join(f"{calls[call]} {label}" for call, label in PLASMID_CALL_LABELS.items() if call in calls)
                    or "no plasmid contig",
                )
            )
        if s in run.isescan:
            table = run.isescan[s]
            facts.append(
                (
                    "IS elements",
                    (
                        f"{len(table)} ({is_complete(table)} complete): {esc(is_family_text(is_families(table), 3))}"
                        if len(table)
                        else "none found"
                    ),
                )
            )
        warnings = run.warnings(s)
        fact_html = "".join(f"<dt>{esc(label)}</dt><dd>{value}</dd>" for label, value in facts)
        warning_html = "".join(f"<li>{esc(w)}</li>" for w in warnings)
        badge = f'<span class="chip status-{esc(status)}">{esc(status)}</span>' if status else ""
        count = (
            f'<span class="flag">{len(warnings)} warning{"s" if len(warnings) > 1 else ""}</span>' if warnings else ""
        )
        name = esc(s.lower())
        cards.append(
            f'<details class="card" id="sample-{slug(s)}" data-name="{name}"><summary>'
            f'<span class="card-name">{esc(s)}</span> {badge} {count}</summary>'
            f'<dl class="facts">{fact_html}</dl>'
            + (f'<ul class="warnings">{warning_html}</ul>' if warnings else "")
            + "</details>"
        )
        title = esc("; ".join(warnings)) if warnings else "no warnings"
        chips.append(
            (
                not warnings,
                f'<li data-name="{name}"><a class="sample-chip{" flagged" if warnings else ""}" '
                f'href="#sample-{slug(s)}" title="{title}">{esc(s)}</a></li>',
            )
        )
    flagged = sum(1 for clean, _ in chips if not clean)
    # samples with warnings first (sorted() keeps name order within each group)
    listing = "".join(chip for _, chip in sorted(chips, key=lambda item: item[0]))
    body = (
        '<div class="sample-finder"><label for="sample-search">Find a sample</label> '
        '<input type="search" id="sample-search" placeholder="Sample name" autocomplete="off"> '
        f'<span class="note">{len(run.samples)} samples, {flagged} with warnings (marked and listed first)</span></div>'
        f'<ul class="sample-list">{listing}</ul>'
        f'<div class="cards">{"".join(cards)}</div>'
    )
    intro = [
        "One card per sample with its key numbers; click a sample in the list, or a sample name anywhere in the "
        "report, to open its card. Samples with warnings are marked with a warning sign and a red outline and "
        "come first in the list. Warnings list what deserves a second look: a skipped or failed step, low "
        "completeness or high contamination, an assembly length far from the usual genome size, a second "
        "species in the reads, species calls that disagree, or an SCCmec call that is not clear-cut.",
    ]
    return Section("samples", "Samples", intro, body)


def section_methods(run):
    if not run.versions and not run.citations:
        return None
    body = ""
    if run.versions:
        rows = [
            {"tool": row["tool"], "version": row["version"], "processes": ", ".join(row["processes"])}
            for row in run.versions
        ]
        body += "<h3>Software versions</h3>" + render_table(
            [("tool", "Tool"), ("version", "Version"), ("processes", "Pipeline steps")], rows
        )
    if run.citations:
        used = {
            normalise_tool(row["tool"])
            for row in run.versions
            if row["processes"] != ["Workflow"] and not re.search(r"(_db|-database)$", row["tool"])
        }
        items = []
        cited = set()
        for entry in run.citations:
            key = normalise_tool(entry["name"])
            if run.versions and key not in ALWAYS_CITED and key not in used:
                continue
            cited.add(key)
            items.append(f'<li><a href="{esc(entry["url"])}">{esc(entry["name"])}</a>. {esc(entry["text"])}</li>')
        body += '<h3>Citations</h3><ul class="citations">' + "".join(items) + "</ul>"
        missing = sorted(tool for tool in used if tool not in cited)
        if missing:
            body += (
                '<p class="note">No reference in CITATIONS.md for: '
                + esc(", ".join(row["tool"] for row in run.versions if normalise_tool(row["tool"]) in missing))
                + ".</p>"
            )
    if run.problems:
        body += "<h3>Inputs that could not be read</h3>" + render_table(
            [("tool", "Tool"), ("file", "File"), ("error", "Error")],
            [{"tool": t, "file": f, "error": e} for t, f, e in run.problems],
        )
    intro = [
        "The exact versions of every tool the run used, read from the pipeline's software_versions.yml, and the "
        "publications to cite when these results are reported.",
    ]
    return Section("methods", "Methods, versions and citations", intro, body)


# Cited for every run (framework and packaging), whatever the versions file lists
ALWAYS_CITED = {"nfcore", "nextflow", "bioconda", "biocontainers", "docker", "singularity", "anaconda"}
# software_versions.yml tool name -> CITATIONS.md name, where they differ beyond punctuation and case
TOOL_ALIASES = {"amrfinder": "amrfinderplus", "mobrecon": "mobsuite", "mobsuite": "mobsuite"}


def gtdb_name(call):
    """GTDB-Tk species in italics; a genome placed only to genus reads "<genus> sp. (no species match)"."""
    if call.get("species"):
        return f"<i>{esc(call['species'])}</i>"
    if call.get("genus"):
        return f"<i>{esc(call['genus'])}</i> sp. (no species match)"
    return esc(call.get("classification", ""))


def normalise_tool(name):
    key = re.sub(r"[^a-z0-9]", "", str(name).lower())
    return TOOL_ALIASES.get(key, key)


# --------------------------------------------------------------------------------------------------
# Per-sample summary table
# --------------------------------------------------------------------------------------------------


def summary_rows(run):
    rows = []
    for s in run.samples:
        status = run.status_row(s)
        reads = run.reads.get(s, {})
        q = run.quast.get(s, {})
        f = run.filter_stats.get(s, {})
        c = run.checkm2.get(s, {})
        m = run.mlst.get(s, {})
        g = run.gtdbtk.get(s, {})
        k = run.kleborate.get(s, {})
        taxa = run.taxa.get(s)
        row = OrderedDict(sample=s)
        row["status"] = status.get("status", "")
        row["stage"] = status.get("stage", "")
        row["reason"] = status.get("reason", "")
        row["reads_before"] = reads.get("reads_before")
        row["bases_before"] = reads.get("bases_before")
        row["reads_after"] = reads.get("reads_after")
        row["bases_after"] = reads.get("bases_after")
        row["q30_rate_after"] = reads.get("q30_rate_after")
        row["mash_genome_size"] = run.genome_size.get(s)
        row["coverage_after"] = (
            round(reads["bases_after"] / run.genome_size[s], 1)
            if run.genome_size.get(s) and reads.get("bases_after")
            else None
        )
        row["reads_top_species"] = taxa.iloc[0]["name"] if taxa is not None and len(taxa) else ""
        row["reads_top_species_percent"] = (
            round(float(taxa.iloc[0]["percent"]), 2) if taxa is not None and len(taxa) else None
        )
        row["reads_second_species"] = taxa.iloc[1]["name"] if taxa is not None and len(taxa) > 1 else ""
        row["reads_second_species_percent"] = (
            round(float(taxa.iloc[1]["percent"]), 2) if taxa is not None and len(taxa) > 1 else None
        )
        row["assembly_contigs"] = q.get("contigs")
        row["assembly_length"] = q.get("total_length")
        row["assembly_n50"] = q.get("n50")
        row["assembly_gc"] = q.get("gc")
        row["annotation_contigs"] = f.get("annotation_contigs")
        row["annotation_length"] = f.get("annotation_length")
        row["checkm2_completeness"] = c.get("completeness")
        row["checkm2_contamination"] = c.get("contamination")
        row["mlst_scheme"] = m.get("scheme") or ""
        row["mlst_st"] = m.get("st") or ""
        row["mlst_organism"] = run.organism(s) or ""
        row["gtdbtk_classification"] = g.get("classification", "")
        row["gtdbtk_ani"] = g.get("ani")
        row["kleborate_species"] = k.get("species", "")
        row["kleborate_virulence_score"] = k.get("virulence_score", "")
        row["kleborate_resistance_score"] = k.get("resistance_score", "")
        row["kleborate_k_locus"] = k.get("K_locus", "")
        row["sistr_serovar"] = run.sistr.get(s, {}).get("serovar", "")
        row["sccmec_type"] = run.sccmec.get(s, {}).get("type", "")
        genes = run.amr_genes(s)
        row["amrfinder_amr_genes"] = ",".join(genes) if s in run.amrfinder else ""
        row["amrfinder_amr_gene_count"] = len(genes) if s in run.amrfinder else None
        amr = run.amrfinder.get(s)
        row["amrfinder_point_mutations"] = (
            ",".join(amr.loc[amr["subtype"] == "POINT", "gene"]) if amr is not None and "subtype" in amr else ""
        )
        rgi = run.rgi.get(s)
        row["rgi_strict_perfect_hits"] = (
            int(rgi["cut_off"].isin(["Perfect", "Strict"]).sum()) if rgi is not None else None
        )
        row["genomad_plasmids"] = len(run.genomad_plasmid[s]) if s in run.genomad_plasmid else None
        row["genomad_viruses"] = len(run.genomad_virus[s]) if s in run.genomad_virus else None
        row["mobsuite_plasmids"] = len(run.mobsuite_mobtyper[s]) if s in run.mobsuite_mobtyper else None
        consensus = run.plasmid_consensus.get(s)
        calls = plasmid_calls(consensus) if consensus is not None else {}
        for call in PLASMID_CALL_LABELS:
            row[f"plasmid_consensus_{call}"] = calls.get(call, 0) if consensus is not None else None
        row["plasmid_replicons"] = plasmid_replicons(consensus) if consensus is not None else ""
        is_table = run.isescan.get(s)
        row["isescan_elements"] = len(is_table) if is_table is not None else None
        row["isescan_complete"] = is_complete(is_table) if is_table is not None else None
        row["isescan_families"] = (
            ",".join(f"{family}:{n}" for family, n in is_families(is_table).items()) if is_table is not None else ""
        )
        row["warnings"] = "; ".join(run.warnings(s))
        rows.append(row)
    return pd.DataFrame(rows)


def tidy_summary(frame):
    """Drop columns no sample has a value for; whole-number columns stay integers when a value is missing."""
    keep = [c for c in frame.columns if c == "sample" or frame[c].replace("", pd.NA).notna().any()]
    frame = frame[keep].copy()
    for column in frame.columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        present = frame[column].notna() & (frame[column] != "")
        if present.any() and values[present].notna().all() and (values[present] % 1 == 0).all():
            if frame[column].map(lambda v: isinstance(v, (int, float)) or v is None).all():
                frame[column] = values.astype("Int64")
    return frame


# --------------------------------------------------------------------------------------------------
# Workbook
# --------------------------------------------------------------------------------------------------

# Excel's cell limit is 32,767 characters, and openpyxl refuses control characters
XLSX_CELL = 32000
XLSX_ILLEGAL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
# RGI sequence columns: long and of little use in a spreadsheet
RGI_SEQUENCE_COLUMNS = ("Predicted_DNA", "Predicted_Protein", "CARD_Protein_Sequence")


def per_sample(tables):
    """{sample: DataFrame} -> one DataFrame with a leading sample column."""
    frames = [frame.assign(sample=s)[["sample"] + list(frame.columns)] for s, frame in tables.items() if len(frame)]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def per_sample_dicts(tables):
    """{sample: dict} -> one row per sample."""
    return pd.DataFrame([{"sample": s, **values} for s, values in tables.items()])


def tool_tables(run):
    """The workbook sheets after the summary: one per tool that ran, in sample name order."""
    order = {s: index for index, s in enumerate(run.samples)}
    reads = per_sample_dicts(run.reads)
    if len(reads):
        reads.insert(1, "tool", reads["sample"].map(run.read_tool))
        if run.genome_size:
            reads["mash_genome_size"] = reads["sample"].map(run.genome_size)
    taxa = per_sample({s: frame.assign(source=run.taxa_source[s]) for s, frame in run.taxa.items()})
    gtdbtk = per_sample_dicts({s: {k: v for k, v in g.items()} for s, g in run.gtdbtk.items()})
    rgi = per_sample(run.rgi)
    rgi = rgi.drop(columns=[c for c in RGI_SEQUENCE_COLUMNS if c in rgi], errors="ignore")
    versions = pd.DataFrame(
        [{"tool": r["tool"], "version": r["version"], "processes": ", ".join(r["processes"])} for r in run.versions]
    )
    sheets = OrderedDict(
        [
            ("Sample status", run.status if run.status is not None else pd.DataFrame()),
            ("Reads", reads),
            ("Species in reads", taxa),
            ("QUAST", per_sample_dicts(run.quast)),
            ("Length filter", per_sample_dicts(run.filter_stats)),
            ("CheckM2", per_sample_dicts(run.checkm2)),
            ("MLST", per_sample_dicts(run.mlst)),
            ("GTDB-Tk", gtdbtk),
            ("Kleborate", per_sample_dicts(run.kleborate)),
            ("SISTR", per_sample_dicts(run.sistr)),
            ("SCCmec", per_sample_dicts(run.sccmec)),
            ("AMRFinderPlus", per_sample(run.amrfinder)),
            ("RGI", rgi),
            ("geNomad plasmids", per_sample(run.genomad_plasmid)),
            ("geNomad viruses", per_sample(run.genomad_virus)),
            ("MOB-suite contigs", per_sample(run.mobsuite_contigs)),
            ("MOB-suite plasmids", per_sample(run.mobsuite_mobtyper)),
            ("Plasmid consensus", per_sample(run.plasmid_consensus)),
            ("ISEScan", per_sample(run.isescan)),
            ("Software versions", versions),
        ]
    )
    for name, frame in sheets.items():
        if len(frame) and "sample" in frame and name != "Sample status":
            # a stable sort keeps each sample's rows in the tool's order
            sheets[name] = frame.sort_values("sample", key=lambda c: c.map(order).fillna(len(order)), kind="stable")
    return OrderedDict((name, frame) for name, frame in sheets.items() if len(frame))


def xlsx_value(value):
    if isinstance(value, str):
        return XLSX_ILLEGAL.sub("", value)[:XLSX_CELL]
    if isinstance(value, (list, tuple, dict)):
        return xlsx_value(json.dumps(value) if isinstance(value, dict) else ", ".join(map(str, value)))
    return value


def write_xlsx(path, summary, sheets):
    """Summary sheet first, then one sheet per tool; header row frozen and filterable."""
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, frame in [("Summary", summary)] + list(sheets.items()):
            # every cell: pandas 3 gives text columns the "str" dtype, not object
            frame = frame.map(xlsx_value)
            frame.to_excel(writer, sheet_name=name[:31], index=False)
            sheet = writer.sheets[name[:31]]
            sheet.freeze_panes = "B2" if "sample" in frame else "A2"
            sheet.auto_filter.ref = sheet.dimensions
            for cells in sheet.iter_cols(min_row=1, max_row=min(sheet.max_row, 200)):
                width = max(len(str(cell.value)) if cell.value is not None else 0 for cell in cells)
                sheet.column_dimensions[cells[0].column_letter].width = min(60, max(8, width + 2))


# --------------------------------------------------------------------------------------------------
# Page
# --------------------------------------------------------------------------------------------------

CSS = """
:root { color-scheme: light; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781; --line: #e1e0d9;
    --surface: #fcfcfb; --page: #f9f9f7; --accent: #2a78d6; --good: #0ca30c; --warn: #fab219; --bad: #d03b3b; }
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink); font: 15px/1.5 system-ui, -apple-system, "Segoe UI",
    Roboto, sans-serif; }
header { background: var(--surface); border-bottom: 1px solid var(--line); padding: 20px 24px 8px; }
header h1 { margin: 0 0 4px; font-size: 24px; }
header p { margin: 0; color: var(--ink2); }
nav { position: sticky; top: 0; z-index: 10; background: var(--surface); border-bottom: 1px solid var(--line);
    padding: 8px 24px; display: flex; flex-wrap: wrap; gap: 4px 16px; font-size: 14px; }
nav a { color: var(--ink2); text-decoration: none; } nav a:hover { color: var(--accent); }
main { max-width: 1280px; margin: 0 auto; padding: 8px 24px 48px; }
section { background: var(--surface); border: 1px solid var(--line); border-radius: 8px; padding: 16px 20px;
    margin: 20px 0; scroll-margin-top: 80px; }
h2 { margin: 0 0 8px; font-size: 20px; } h3 { margin: 20px 0 8px; font-size: 16px; }
.intro p { margin: 0 0 8px; max-width: 80ch; color: var(--ink2); }
a { color: var(--accent); }
.table-wrap { overflow: auto; margin: 8px 0 12px; border: 1px solid var(--line); border-radius: 6px; }
table.data { border-collapse: collapse; width: 100%; font-size: 13px; }
table.data th, table.data td { padding: 5px 10px; border-bottom: 1px solid var(--line); text-align: left;
    vertical-align: top; }
table.data th { position: sticky; top: 0; background: #f0efec; font-weight: 600; cursor: pointer;
    white-space: nowrap; }
table.data th.num, table.data td.num { text-align: right; font-variant-numeric: tabular-nums; }
table.data th[aria-sort="ascending"]::after { content: " \\25B2"; }
table.data th[aria-sort="descending"]::after { content: " \\25BC"; }
table.data tbody tr:hover { background: #f4f3ef; }
caption { text-align: left; padding: 6px 10px; color: var(--ink2); }
.figure { margin: 8px 0; } .scroll-x { overflow-x: auto; }
dl.facts { display: grid; grid-template-columns: max-content 1fr; gap: 2px 16px; margin: 8px 0; font-size: 14px; }
dl.facts dt { color: var(--ink2); } dl.facts dd { margin: 0; }
.chips { display: flex; gap: 8px; flex-wrap: wrap; }
.chip { display: inline-block; border-radius: 999px; padding: 1px 10px; font-size: 13px; font-weight: 600;
    border: 1px solid var(--line); background: #f0efec; vertical-align: middle; }
.chip.status-pass::before { content: "\\2713  "; color: var(--good); }
.chip.status-skipped::before { content: "\\2013  "; color: #8a6100; }
.chip.status-failed::before { content: "\\2715  "; color: var(--bad); }
.flag { color: var(--bad); font-weight: 600; } .flag::before { content: "\\26A0  "; }
.cards { display: grid; grid-template-columns: repeat(auto-fill, minmax(380px, 1fr)); gap: 12px;
    align-items: start; }
.card { border: 1px solid var(--line); border-radius: 8px; padding: 8px 14px; background: #fff;
    scroll-margin-top: 80px; }
.card:target { outline: 2px solid var(--accent); }
.card summary { cursor: pointer; font-weight: 600; }
.card summary .chip, .card summary .flag { font-weight: 600; font-size: 13px; margin-left: 6px; }
.card[open] summary { margin-bottom: 4px; }
.sample-finder { display: flex; flex-wrap: wrap; align-items: center; gap: 8px; margin: 8px 0; }
.sample-finder input { font: inherit; padding: 4px 8px; border: 1px solid #c3c2b7; border-radius: 6px;
    min-width: 220px; }
ul.sample-list { list-style: none; padding: 0; margin: 8px 0 16px; display: flex; flex-wrap: wrap; gap: 6px;
    max-height: 240px; overflow: auto; }
.sample-chip { display: inline-block; padding: 1px 10px; border: 1px solid var(--line); border-radius: 999px;
    font-size: 13px; text-decoration: none; color: var(--ink); background: #f0efec; }
.sample-chip.flagged { border: 2px solid var(--bad); background: #fbeaea; color: #8a2b2b; font-weight: 600; }
.sample-chip.flagged::before { content: "\\26A0  "; }
.sample-chip:hover { border-color: var(--accent); }
[hidden] { display: none !important; }
.card dl.facts { font-size: 13px; }
ul.warnings { margin: 6px 0 0; padding-left: 20px; font-size: 13px; color: #8a2b2b; }
ul.citations { font-size: 13px; padding-left: 20px; } .note { font-size: 13px; color: var(--ink2); } ul.citations li { margin-bottom: 6px; }
footer { color: var(--muted); font-size: 12px; text-align: center; padding: 16px; }
@media (max-width: 600px) { main, header, nav { padding-left: 12px; padding-right: 12px; }
    nav { position: static; } section, .card { scroll-margin-top: 8px; }
    .cards { grid-template-columns: 1fr; } }
@media print { nav, .sample-finder, ul.sample-list { display: none; } section { break-inside: avoid-page; } }
"""

CARDS_JS = """
function openCard(id) {
    var card = document.getElementById(id);
    if (card && card.tagName === "DETAILS") { card.open = true; }
}
document.querySelectorAll('a[href^="#sample-"]').forEach(function (link) {
    link.addEventListener("click", function () { openCard(link.getAttribute("href").slice(1)); });
});
window.addEventListener("hashchange", function () { openCard(location.hash.slice(1)); });
if (location.hash) { openCard(location.hash.slice(1)); }
var search = document.getElementById("sample-search");
if (search) {
    search.addEventListener("input", function () {
        var query = search.value.trim().toLowerCase();
        document.querySelectorAll(".sample-list li, .cards .card").forEach(function (item) {
            item.hidden = query !== "" && item.getAttribute("data-name").indexOf(query) === -1;
        });
    });
}
"""

SORT_JS = """
document.querySelectorAll("table.sortable").forEach(function (table) {
    table.querySelectorAll("th").forEach(function (th, index) {
        th.addEventListener("click", function () {
            var body = table.tBodies[0];
            var ascending = th.getAttribute("aria-sort") !== "ascending";
            table.querySelectorAll("th").forEach(function (other) { other.removeAttribute("aria-sort"); });
            th.setAttribute("aria-sort", ascending ? "ascending" : "descending");
            var numeric = th.classList.contains("num");
            var rows = Array.prototype.slice.call(body.rows);
            rows.sort(function (a, b) {
                var x = a.cells[index].textContent.trim(), y = b.cells[index].textContent.trim();
                if (numeric) {
                    x = x === "" ? -Infinity : parseFloat(x.replace(/,/g, ""));
                    y = y === "" ? -Infinity : parseFloat(y.replace(/,/g, ""));
                    return ascending ? x - y : y - x;
                }
                return ascending ? x.localeCompare(y, undefined, {numeric: true}) : y.localeCompare(x, undefined, {numeric: true});
            });
            rows.forEach(function (row) { body.appendChild(row); });
        });
    });
});
"""

# A plain ring as the page icon, so the browser does not ask the network for /favicon.ico
FAVICON = (
    "data:image/svg+xml;base64,"
    + base64.b64encode(
        b'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16"><circle cx="8" cy="8" r="6" fill="none" '
        b'stroke="#2a78d6" stroke-width="3"/></svg>'
    ).decode()
)


def build_html(run, sections):
    info = run.run_info
    subtitle = " · ".join(
        str(part)
        for part in (
            f"MGAP {info['pipeline_version']}" if info.get("pipeline_version") else "MGAP",
            info.get("date"),
            info.get("seq_type"),
            f"{len(run.samples)} samples",
        )
        if part
    )
    nav = "".join(f'<a href="#{s.key}">{esc(s.title)}</a>' for s in sections)
    content = "".join(section.html() for section in sections)
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f'<title>MGAP run report</title><link rel="icon" href="{FAVICON}">'
        f"<style>{CSS}</style><script>{get_plotlyjs()}</script></head><body>"
        f"<header><h1>MGAP run report</h1><p>{esc(subtitle)}</p></header>"
        f"<nav>{nav}</nav><main>{content}</main>"
        f"<footer>Generated {date.today().isoformat()} by mgap_report.py. The page works offline.</footer>"
        f"<script>{SORT_JS}{CARDS_JS}</script></body></html>"
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--manifest", required=True, help="TSV or JSON: sample, tool, file")
    parser.add_argument("--run-info", help="JSON with run metadata (version, date, seq_type, parameters ...)")
    parser.add_argument("--html", default="mgap_report.html", help="output HTML (default: %(default)s)")
    parser.add_argument("--tsv", default="mgap_summary.tsv", help="output per-sample TSV (default: %(default)s)")
    parser.add_argument("--xlsx", help="output XLSX workbook: the summary, then one sheet per tool (needs openpyxl)")
    args = parser.parse_args(argv)

    run_info = {}
    if args.run_info:
        with open(args.run_info) as handle:
            run_info = json.load(handle)
    run = RunData(read_manifest(args.manifest), run_info)
    builders = (
        section_overview,
        section_read_qc,
        section_assembly,
        section_quality,
        section_typing,
        section_amr,
        section_mobile,
        section_plasmids,
        section_isescan,
        section_cards,
        section_methods,
    )
    sections = [section for section in (builder(run) for builder in builders) if section is not None]
    Path(args.html).write_text(build_html(run, sections), encoding="utf-8")
    summary = tidy_summary(summary_rows(run))
    summary.to_csv(args.tsv, sep="\t", index=False)
    outputs = [args.tsv]
    if args.xlsx:
        sheets = tool_tables(run)
        write_xlsx(args.xlsx, summary, sheets)
        outputs.append(f"{args.xlsx} ({len(sheets) + 1} sheets)")
    print(f"{args.html}: {len(sections)} sections, {len(run.samples)} samples; {'; '.join(outputs)}", file=sys.stderr)


if __name__ == "__main__":
    main()
