#!/usr/bin/env python3
"""Build the skani species reference database that MGAP's --species_refs points to.

Reads the accession list (assets/species_refs.tsv: accession, species, strain, amrfinder_organism,
mlst_schemes), downloads each genome from NCBI with the `datasets` command-line tool, sketches them with
`skani sketch` and writes, into the output folder:

- genomes/<accession>.fna   the downloaded genomes (one FASTA per accession)
- sketch/                   the skani database (`skani search -d`)
- species_refs.tsv          a copy of the accession list: the pipeline reads species, organism and schemes
                            from this copy, so a database always matches the list it was built from
- version.txt               build date, skani and datasets versions, accession count, MD5 of the list

The database must be sketched by the skani version the pipeline runs (0.4.0, the SKANI_SEARCH container):
sketches written by skani 0.4 cannot be read by 0.3 or older. The script stops when `skani --version` is
not that version, unless --skani-version says otherwise. skani does not have to be installed: --skani takes
any command, for example the pipeline's container (the output folder must be mounted, and the script runs
skani from inside it):

    build_species_refs.py --refs assets/species_refs.tsv --outdir /dbs/skani/species_refs/2026-10-29 \\
        --skani "docker run --rm -u $(id -u):$(id -g) -v /dbs/skani/species_refs/2026-10-29:/db -w /db \\
            quay.io/biocontainers/skani:0.4.0--hec9b1f2_0 skani"

Nothing is downloaded by the pipeline itself; this script is run once, by hand, from an approved list.
"""

import argparse
import csv
import datetime
import hashlib
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

COLUMNS = ["accession", "species", "strain", "amrfinder_organism", "mlst_schemes"]
SKANI_VERSION = "0.4.0"


def read_refs(path):
    """The accession list as a list of dicts; stops on a missing column or a repeated accession or species."""
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"{path}: missing column(s) {', '.join(missing)}")
        rows = [row for row in reader if row["accession"].strip()]
    for column in ("accession", "species"):
        values = [row[column] for row in rows]
        repeated = sorted({value for value in values if values.count(value) > 1})
        if repeated:
            sys.exit(f"{path}: repeated {column}: {', '.join(repeated)}")
    return rows


def run(command, **kwargs):
    # the NCBI API key is not echoed
    shown = ["***" if index and command[index - 1] == "--api-key" else part for index, part in enumerate(command)]
    print("+ " + " ".join(shlex.quote(part) for part in shown), file=sys.stderr, flush=True)
    return subprocess.run(command, check=True, **kwargs)


def tool_version(command):
    output = run(command + ["--version"], capture_output=True, text=True)
    return (output.stdout or output.stderr).strip().split()[-1]


def md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(accessions, genome_dir, datasets, api_key=None):
    """Download the genomes with NCBI datasets into genome_dir/<accession>.fna."""
    with tempfile.TemporaryDirectory(dir=genome_dir.parent) as tmp:
        tmp = Path(tmp)
        (tmp / "accessions.txt").write_text("\n".join(accessions) + "\n")
        command = datasets + ["download", "genome", "accession", "--inputfile", str(tmp / "accessions.txt")]
        command += ["--include", "genome", "--filename", str(tmp / "genomes.zip")]
        if api_key:
            command += ["--api-key", api_key]
        run(command)
        with zipfile.ZipFile(tmp / "genomes.zip") as archive:
            archive.extractall(tmp / "unzipped")
        for accession in accessions:
            found = sorted((tmp / "unzipped" / "ncbi_dataset" / "data" / accession).glob("*_genomic.fna"))
            if len(found) != 1:
                sys.exit(f"{accession}: expected one *_genomic.fna in the NCBI download, found {len(found)}")
            shutil.move(str(found[0]), genome_dir / f"{accession}.fna")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--refs", required=True, type=Path, help="accession list (assets/species_refs.tsv)")
    parser.add_argument("--outdir", required=True, type=Path, help="new database folder (must not exist)")
    parser.add_argument("--skani", default="skani", help="skani command, run inside --outdir (default: skani)")
    parser.add_argument("--skani-version", default=SKANI_VERSION, help=f"required skani version ({SKANI_VERSION})")
    parser.add_argument("--datasets", default="datasets", help="NCBI datasets command (default: datasets)")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"), help="NCBI API key (or NCBI_API_KEY)")
    parser.add_argument("--threads", type=int, default=4, help="skani sketch threads (default: 4)")
    args = parser.parse_args(argv)

    rows = read_refs(args.refs)
    if args.outdir.exists():
        sys.exit(f"{args.outdir} exists: give a new folder (one database per build)")
    skani = shlex.split(args.skani)
    datasets = shlex.split(args.datasets)

    genome_dir = args.outdir / "genomes"
    genome_dir.mkdir(parents=True)
    # skani runs inside the output folder, so a container only needs that folder mounted
    skani_version = tool_version(skani)
    if skani_version != args.skani_version:
        sys.exit(f"skani {skani_version} found, {args.skani_version} required (the pipeline's SKANI_SEARCH)")
    datasets_version = tool_version(datasets)

    accessions = [row["accession"] for row in rows]
    download(accessions, genome_dir, datasets, args.api_key)

    # relative paths, so skani search reports genomes/<accession>.fna whatever the folder is mounted as
    fasta_list = args.outdir / "genomes.txt"
    fasta_list.write_text("".join(f"genomes/{accession}.fna\n" for accession in accessions))
    run(skani + ["sketch", "-l", "genomes.txt", "-o", "sketch", "-t", str(args.threads)], cwd=args.outdir)

    shutil.copyfile(args.refs, args.outdir / "species_refs.tsv")
    version = {
        "database": "MGAP species reference set (skani)",
        "built": datetime.date.today().isoformat(),
        "skani_version": skani_version,
        "datasets_version": datasets_version,
        "genomes": str(len(accessions)),
        "species_refs_md5": md5(args.outdir / "species_refs.tsv"),
    }
    (args.outdir / "version.txt").write_text("".join(f"{key}: {value}\n" for key, value in version.items()))
    print(f"Built {args.outdir}: {len(accessions)} genomes, skani {skani_version}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        sys.exit(f"{error.cmd[0]} failed (exit {error.returncode}); remove the partial --outdir folder and run again")
