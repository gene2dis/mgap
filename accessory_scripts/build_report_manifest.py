#!/usr/bin/env python3
"""Build the manifest and run-info JSON that bin/mgap_report.py reads, from a finished MGAP results folder.

The pipeline gives the report its files directly; this helper is for re-running the report on an
existing results folder (and for developing it). It looks for each tool's files in the published
layout ({outdir}/<sample>/..., gtdbtk/, pipeline_info/) and writes:

- a manifest TSV (sample, tool, file) with absolute paths
- a run-info JSON: pipeline and Nextflow versions (software_versions.yml), run date (execution trace
    timestamp), seq_type (sample_status.tsv), the parameters and databases given with --param, and the
    MLST scheme to AMRFinderPlus organism map read from the pipeline's workflows/mgap.nf

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
    "kleborate": ["annotation/kleborate/*_output.txt"],
    "sistr": ["annotation/sistr/*.tab"],
    "sccmec": ["annotation/sccmec/{sample}.tsv"],
}
# run-level tool -> glob below {outdir}
RUN_PATTERNS = {
    "gtdbtk": ["gtdbtk/*/*.summary.tsv"],
    "sample_status": ["pipeline_info/sample_status.tsv"],
    "software_versions": ["pipeline_info/software_versions.yml"],
}
EXCLUDE = re.compile(r"hAMRonization")
NON_SAMPLE_DIRS = {"pipeline_info", "gtdbtk", "multiqc", "report"}
DATABASE_PARAMS = re.compile(r"(_db|db)$")


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
            for path in find(sample_dir, patterns, sample_dir.name):
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
        workflow = (yaml.safe_load(versions.read_text()) or {}).get("Workflow", {})
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
