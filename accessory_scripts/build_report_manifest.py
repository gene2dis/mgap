#!/usr/bin/env python3
"""Build the manifest and run-info JSON that bin/mgap_report.py reads, from a finished MGAP results folder.

The pipeline gives the report its files directly; this helper is for re-running the report on an
existing results folder (and for developing it). It looks for each tool's files in the published
layout ({outdir}/<sample>/..., gtdbtk/, pipeline_info/) and writes:

- a manifest TSV (sample, tool, file) with absolute paths
- a run-info JSON: pipeline and Nextflow versions (software_versions.yml), run date (execution trace
    timestamp), seq_type (sample_status.tsv), the parameters and databases given with --param, the
    database versions read from each database folder's own version file (Bakta, AMRFinderPlus, geNomad,
    GTDB-Tk, CARD for RGI; Kraken2, CheckM2 and MOB-suite databases record none), and the MLST scheme to
    AMRFinderPlus organism map read from the pipeline's workflows/mgap.nf

Example:
    build_report_manifest.py --results results --pipeline . --param min_bases=50000000 \\
        --param bakta_db=/dbs/bakta --manifest manifest.tsv --run-info run_info.json
    ../bin/mgap_report.py --manifest manifest.tsv --run-info run_info.json --html report.html
"""

import argparse
import csv
import glob
import json
import re
from pathlib import Path

import yaml

# tool -> glob below {outdir}/<sample>/; {sample} is the sample name (where a folder holds other files too)
SAMPLE_PATTERNS = {
    "fastp": ["read_processing/fastp/*.fastp.json"],
    "fastplong": ["read_processing/fastplong/*.json"],
    "mash": ["read_processing/mash/*.genome_size"],
    "kraken2": ["read_processing/kraken2/*.report.txt"],
    "bracken": ["read_processing/bracken/*.tsv"],
    "quast": ["qc/quast/*.tsv"],
    "filter_stats": ["assemblies/annotation_input/*.filter_stats.tsv"],
    "checkm2": ["annotation/checkm2/*_checkm2_report.tsv"],
    "mlst": ["annotation/mlst/*.tsv"],
    "amrfinderplus": ["annotation/amrfinder/{sample}.tsv"],
    "rgi": ["annotation/rgi/{sample}.txt"],
    "genomad_plasmid": ["annotation/genomad/**/*_plasmid_summary.tsv"],
    "genomad_virus": ["annotation/genomad/**/*_virus_summary.tsv"],
    "mobsuite_contigs": ["annotation/mobsuite/contig_report.txt"],
    "mobsuite_mobtyper": ["annotation/mobsuite/mobtyper_results.txt"],
    "plasmid_consensus": ["annotation/plasmid_consensus/*.plasmid_consensus.tsv"],
    "isescan": ["annotation/isescan/*.isescan.tsv"],
    "kleborate": ["annotation/kleborate/*_output.txt"],
    "sistr": ["annotation/sistr/*.tab"],
    "sccmec": ["annotation/sccmec/{sample}.tsv"],
}
# run-level tool -> glob below {outdir}
RUN_PATTERNS = {
    "gtdbtk": ["gtdbtk/*/*.summary.tsv"],
    "sample_status": ["pipeline_info/sample_status.tsv"],
    "species_check": ["pipeline_info/species_check.tsv"],
    "software_versions": ["pipeline_info/software_versions.yml"],
}
EXCLUDE = re.compile(r"hAMRonization")
NON_SAMPLE_DIRS = {"pipeline_info", "gtdbtk", "multiqc", "report"}
# database paths, as the pipeline lists them in run_info.json (mlst_datadir included)
DATABASE_PARAMS = re.compile(r"(db|_datadir|^species_refs)$")


def bakta_version(base):
    data = json.loads((base / "version.json").read_text())
    return f"{data['major']}.{data['minor']} {data.get('type', '')} ({data['date']})".replace("  ", " ")


def gtdbtk_version(base):
    match = re.search(r"VERSION_DATA=(\S+)", (base / "metadata" / "metadata.txt").read_text())
    return match.group(1) if match else None


def card_version(base):
    names = [
        re.search(r"card_database_v([\d.]+?)(?:_all)?\.fasta$", path.name)
        for path in base.glob("card_database_v*.fasta")
    ]
    # numeric order (4.0.10 after 4.0.9), as the pipeline's databaseVersion() sorts them
    versions = sorted({match.group(1) for match in names if match}, key=lambda v: [int(p) for p in v.split(".") if p])
    return f"CARD {versions[-1]}" if versions else None


def species_refs_version(base):
    """version.txt of accessory_scripts/build_species_refs.py: built date, genome count, skani version."""
    fields = dict(line.split(": ", 1) for line in (base / "version.txt").read_text().splitlines() if ": " in line)
    if not fields.get("built"):
        return None
    return f"{fields['built']} ({fields.get('genomes')} genomes, skani {fields.get('skani_version')})"


# database parameter -> reader of the version file in its folder
DATABASE_VERSIONS = {
    "bakta_db": bakta_version,
    "amrfinder_db": lambda base: (base / "version.txt").read_text().strip(),
    "genomad_db": lambda base: (base / "version.txt").read_text().strip(),
    "species_refs": species_refs_version,
    "gtdbtk_db": gtdbtk_version,
    "rgi_db": card_version,
}


def database_versions(databases):
    versions = {}
    for param, path in databases.items():
        base = Path(path)
        if param not in DATABASE_VERSIONS or not base.exists():
            continue
        try:
            version = DATABASE_VERSIONS[param](base if base.is_dir() else base.parent)
        except (OSError, ValueError, KeyError):
            version = None
        if version:
            versions[param] = version
    return versions


def find(base, patterns, sample=""):
    files = []
    for pattern in patterns:
        files += [
            path
            for path in sorted(base.glob(pattern.format(sample=glob.escape(sample))))
            if path.is_file() and not EXCLUDE.search(path.name)
        ]
    return files


def manifest_rows(results, citations):
    rows = []
    for sample_dir in sorted(path for path in results.iterdir() if path.is_dir()):
        if sample_dir.name in NON_SAMPLE_DIRS:
            continue
        for tool, patterns in SAMPLE_PATTERNS.items():
            files = find(sample_dir, patterns, sample_dir.name)
            if tool == "mlst":
                # after a species-check override, <sample>.<scheme>.tsv (the re-run) is what the report shows
                rerun = [path for path in files if path.name != f"{sample_dir.name}.tsv"]
                files = rerun or files
            for path in files:
                rows.append({"sample": sample_dir.name, "tool": tool, "file": str(path.resolve())})
    for tool, patterns in RUN_PATTERNS.items():
        for path in find(results, patterns):
            rows.append({"sample": "-", "tool": tool, "file": str(path.resolve())})
    if citations and citations.exists():
        rows.append({"sample": "-", "tool": "citations", "file": str(citations.resolve())})
    return rows


def mlst_organisms(workflow_file):
    """The getTaxaNames() map ("scheme": "Organism") from workflows/mgap.nf."""
    text = workflow_file.read_text()
    body = re.search(r"def getTaxaNames\(\)\s*\{(.*?)\n\}", text, re.S)
    return dict(re.findall(r'"([^"]+)"\s*:\s*"([^"]+)"', body.group(1))) if body else {}


def run_info(results, pipeline, params):
    info = {}
    versions = results / "pipeline_info" / "software_versions.yml"
    if versions.exists():
        workflow = (yaml.load(versions.read_text(), Loader=yaml.BaseLoader) or {}).get("Workflow", {})
        for name, version in workflow.items():
            if name.lower() == "nextflow":
                info["nextflow_version"] = str(version)
            else:
                info["pipeline_version"] = str(version)
    traces = sorted((results / "pipeline_info").glob("execution_trace_*.txt"))
    if traces:
        stamp = re.search(r"(\d{4}-\d{2}-\d{2})", traces[-1].name)
        if stamp:
            info["date"] = stamp.group(1)
    status = results / "pipeline_info" / "sample_status.tsv"
    if status.exists():
        with open(status, newline="") as handle:
            seq_types = sorted({row["seq_type"] for row in csv.DictReader(handle, delimiter="\t")})
        info["seq_type"] = ", ".join(seq_types)
    info["parameters"] = {k: v for k, v in params.items() if not DATABASE_PARAMS.search(k)}
    info["databases"] = {k: v for k, v in params.items() if DATABASE_PARAMS.search(k)}
    info["database_versions"] = database_versions(info["databases"])
    workflow_file = pipeline / "workflows" / "mgap.nf"
    if workflow_file.exists():
        info["mlst_organisms"] = mlst_organisms(workflow_file)
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--results", required=True, type=Path, help="MGAP --outdir of a finished run")
    parser.add_argument("--pipeline", type=Path, default=Path(__file__).resolve().parent.parent, help="MGAP checkout")
    parser.add_argument("--param", action="append", default=[], help="name=value shown in the run overview")
    parser.add_argument("--manifest", default="manifest.tsv", help="output manifest TSV")
    parser.add_argument("--run-info", default="run_info.json", help="output run-info JSON")
    args = parser.parse_args()

    params = dict(item.split("=", 1) for item in args.param)
    rows = manifest_rows(args.results, args.pipeline / "CITATIONS.md")
    with open(args.manifest, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "tool", "file"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with open(args.run_info, "w") as handle:
        json.dump(run_info(args.results, args.pipeline, params), handle, indent=2)
        handle.write("\n")
    print(f"{args.manifest}: {len(rows)} files; {args.run_info}")


if __name__ == "__main__":
    main()
