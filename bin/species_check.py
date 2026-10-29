#!/usr/bin/env python3
"""One species decision per sample: which AMRFinderPlus organism and species tools the sample gets.

Combines the MLST scheme (mlst's automatic choice), the skani call against the curated reference set
(--species_refs, accessory_scripts/build_species_refs.py), the top species in the reads (Bracken, or
Kraken2 for ONT) and the species given in the samplesheet. Rules (D-13, D-14, D-15; docs/output.md):

- skani call: the best reference hit with ANI >= --min-ani and query aligned fraction >= --min-af.
    When the next hit of another organism is within --ani-margin ANI, the call is ambiguous.
- agree: the reference species is typed by the MLST scheme (its mlst_schemes) and maps to the same
    organism as the scheme. refine: typed by the scheme, but the species maps to another organism (or
    none), e.g. Bordetella pertussis under bordetella_3 (no organism of its own, D-14) or Staphylococcus
    argenteus under saureus. The organism follows the reference species.
- override: the scheme does not type the reference species (a K. pneumoniae typed as
    ecoli_achtman_4); routing follows the reference species and mlst is run again with its first scheme.
- reference: mlst found no scheme; routing follows the reference species (mlst again with its scheme).
- ambiguous: two references of different organisms within the margin; the MLST organism is kept when
    it is one of them, else none.
- mlst: no reference call (no --species_refs, or no hit); routing follows the MLST scheme. With a
    read classification, its genus is compared with the MLST genus (a difference is flagged).
- samplesheet: no reference call and no MLST organism; the samplesheet species gives the organism.

Writes one TSV row (with header). The pipeline fills rerun_st after the mlst re-run.
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

COLUMNS = [
    "sample",
    "seq_type",
    "mlst_scheme",
    "mlst_st",
    "mlst_organism",
    "skani_species",
    "skani_accession",
    "skani_ani",
    "skani_af_ref",
    "skani_af_query",
    "skani_next",
    "reads_species",
    "reads_percent",
    "reads_source",
    "expected_species",
    "species",
    "organism",
    "rerun_scheme",
    "rerun_st",
    "decision",
    "flagged",
    "reason",
]
REF_COLUMNS = ["accession", "species", "strain", "amrfinder_organism", "mlst_schemes"]
FASTA_SUFFIX = re.compile(r"\.(fna|fa|fasta)(\.gz)?$")


def none_if_dash(value):
    value = (value or "").strip()
    return None if value in ("", "-") else value


def read_refs(path):
    """species_refs.tsv -> {accession: {species, organism, schemes}}."""
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in REF_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing column(s) {', '.join(missing)}")
        refs = {}
        for row in reader:
            schemes = none_if_dash(row["mlst_schemes"])
            refs[row["accession"].strip()] = {
                "species": row["species"].strip(),
                "organism": none_if_dash(row["amrfinder_organism"]),
                "schemes": schemes.split(";") if schemes else [],
            }
    return refs


def read_mlst(path):
    """mlst line: file, scheme, ST, alleles; '-' (or an empty file) is no scheme / no ST."""
    with open(path) as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if not lines:
        return None, None
    fields = lines[-1]
    scheme = none_if_dash(fields[1]) if len(fields) > 1 else None
    st = none_if_dash(fields[2]) if len(fields) > 2 else None
    return scheme, st


def accession_of(ref_file):
    return FASTA_SUFFIX.sub("", Path(ref_file).name)


def read_skani(path):
    """skani search table -> hits sorted by ANI (highest first)."""
    hits = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hits.append(
                {
                    "accession": accession_of(row["Ref_file"]),
                    "ani": float(row["ANI"]),
                    "af_ref": float(row["Align_fraction_ref"]),
                    "af_query": float(row["Align_fraction_query"]),
                }
            )
    return sorted(hits, key=lambda hit: -hit["ani"])


def read_bracken(path):
    """Bracken species table -> (top species, percent of the classified reads)."""
    best = None
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            percent = float(row["fraction_total_reads"]) * 100
            if best is None or percent > best[1]:
                best = (row["name"].strip(), percent)
    return best


def read_kraken2(path):
    """Kraken2 report -> (top species-rank taxon, percent of all reads)."""
    best = None
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[-3] != "S":
                continue
            percent = float(fields[0])
            if best is None or percent > best[1]:
                best = (fields[-1].strip(), percent)
    return best


# Shigella is genomically Escherichia coli (one AMRFinderPlus organism); read classifiers report either
GENUS_ALIASES = {"Shigella": "Escherichia"}


def genus(name):
    if not name:
        return None
    first = name.replace("_", " ").split()[0]
    return GENUS_ALIASES.get(first, first)


def same_species(expected, names):
    """Samplesheet species against called species; a genus alone is compared at genus level."""
    if len(expected.split()) == 1:
        return genus(expected) in {genus(name) for name in names}
    return expected.lower() in (name.lower() for name in names)


def scheme_genera(scheme, refs, mlst_organism):
    """Genera the MLST scheme points to: of the references it types, else of its organism."""
    genera = {genus(ref["species"]) for ref in refs.values() if scheme in ref["schemes"]}
    if not genera and mlst_organism:
        genera = {genus(mlst_organism)}
    return genera


def species_call(hits, refs, min_ani, min_af, margin):
    """(call, runner-up) from the skani hits; call None when no hit passes; runner-up only if ambiguous."""
    passing = [hit for hit in hits if hit["ani"] >= min_ani and hit["af_query"] >= min_af and hit["accession"] in refs]
    if not passing:
        return None, None
    top = passing[0]
    for other in passing[1:]:
        if refs[other["accession"]]["organism"] == refs[top["accession"]]["organism"]:
            continue
        if top["ani"] - other["ani"] < margin:
            return top, other
        break
    return top, None


def decide(args):
    refs = read_refs(args.refs)
    with open(args.mlst_organisms) as handle:
        mlst_map = json.load(handle)
    scheme, st = read_mlst(args.mlst)
    mlst_organism = mlst_map.get(scheme) if scheme else None
    row = dict.fromkeys(COLUMNS, "")
    row.update(sample=args.sample, seq_type=args.seq_type, mlst_scheme=scheme or "-", mlst_st=st or "-")
    row["mlst_organism"] = mlst_organism or ""
    notes = []
    flagged = False

    # read classification (Illumina: Bracken; ONT: Kraken2)
    reads = None
    if args.taxa:
        reads = read_bracken(args.taxa) if args.taxa_format == "bracken" else read_kraken2(args.taxa)
        if reads:
            row.update(
                reads_species=reads[0],
                reads_percent=f"{reads[1]:.1f}",
                reads_source="Bracken" if args.taxa_format == "bracken" else "Kraken2",
            )

    # expected species from the samplesheet: matched to a reference species by name
    expected = (args.expected or "").strip()
    expected_ref = None
    if expected:
        row["expected_species"] = expected
        expected_ref = next((ref for ref in refs.values() if ref["species"].lower() == expected.lower()), None)

    call = runner_up = None
    if args.skani:
        hits = read_skani(args.skani)
        call, runner_up = species_call(hits, refs, args.min_ani, args.min_af, args.ani_margin)
        best = call or (hits[0] if hits else None)
        if best and best["accession"] in refs:
            row.update(
                skani_species=refs[best["accession"]]["species"],
                skani_accession=best["accession"],
                skani_ani=f"{best['ani']:.2f}",
                skani_af_ref=f"{best['af_ref']:.2f}",
                skani_af_query=f"{best['af_query']:.2f}",
            )
        if runner_up:
            row["skani_next"] = f"{refs[runner_up['accession']]['species']} ({runner_up['ani']:.2f})"

    organism = mlst_organism
    species = None
    if call and runner_up:
        candidates = {refs[call["accession"]]["organism"], refs[runner_up["accession"]]["organism"]}
        organism = mlst_organism if mlst_organism in candidates else None
        species = f"{refs[call['accession']]['species']} / {refs[runner_up['accession']]['species']}"
        decision = "ambiguous"
        flagged = True
        notes.append(
            f"skani cannot separate {species} (ANI {call['ani']:.2f} vs {runner_up['ani']:.2f}, "
            f"margin {args.ani_margin:g}); "
            + (f"MLST organism {mlst_organism} kept" if organism else "no AMRFinderPlus organism")
        )
    elif call:
        ref = refs[call["accession"]]
        species = ref["species"]
        organism = ref["organism"]
        label = f"skani {species} (ANI {call['ani']:.2f}, AF {call['af_query']:.1f} %)"
        if not scheme:
            decision = "reference"
            row["rerun_scheme"] = ref["schemes"][0] if ref["schemes"] else ""
            notes.append(f"no MLST scheme; {label}")
        elif scheme in ref["schemes"]:
            decision = "agree" if organism == mlst_organism else "refine"
            if decision == "refine":
                notes.append(
                    f"{label} is typed by {scheme} but its organism is {organism or 'none'} "
                    f"(MLST alone: {mlst_organism or 'none'})"
                )
                flagged = mlst_organism is not None
        else:
            decision = "override"
            flagged = True
            row["rerun_scheme"] = ref["schemes"][0] if ref["schemes"] else ""
            rerun = f"; mlst re-run with {row['rerun_scheme']}" if row["rerun_scheme"] else ""
            notes.append(f"MLST scheme {scheme} does not type {label}: routed by the reference{rerun}")
    else:
        decision = "mlst"
        if not args.skani:
            notes.append("no reference set (--species_refs)")
        else:
            flagged = True
            best = row["skani_species"]
            notes.append(
                "no reference hit"
                + (
                    f" (best: {best}, ANI {row['skani_ani']}, AF {row['skani_af_query']} %)"
                    if best
                    else " (no reference within the skani screen)"
                )
            )
        if organism is None and expected_ref and expected_ref["organism"]:
            organism = expected_ref["organism"]
            species = expected_ref["species"]
            decision = "samplesheet"
            notes.append(f"organism from the samplesheet species {expected_ref['species']}")

    # second opinions: read classification and samplesheet species (genus level)
    called_genus = genus(species.split(" / ")[0]) if species else None
    reference_genera = (
        {called_genus} if called_genus else scheme_genera(scheme, refs, mlst_organism) if scheme else set()
    )
    if reads and reference_genera:
        if genus(reads[0]) in reference_genera:
            if decision == "mlst":
                notes.append(f"{row['reads_source']} {reads[0]} agrees with the MLST genus")
        else:
            flagged = True
            notes.append(f"{row['reads_source']} top species {reads[0]} ({reads[1]:.1f} %) is another genus")
    if expected and decision != "samplesheet":
        compare = species.split(" / ") if species else None
        if compare:
            if not same_species(expected, compare):
                flagged = True
                notes.append(f"samplesheet species {expected} differs")
        elif reference_genera and genus(expected) not in reference_genera:
            flagged = True
            notes.append(f"samplesheet species {expected} differs from the MLST genus")
        if not expected_ref and len(expected.split()) > 1:
            notes.append(f"samplesheet species {expected} is not in the reference list")

    row.update(
        species=species or "",
        organism=organism or "",
        decision=decision,
        flagged="yes" if flagged else "no",
        reason="; ".join(notes),
    )
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--seq-type", required=True)
    parser.add_argument("--mlst", required=True, help="mlst TSV (automatic scheme)")
    parser.add_argument("--refs", required=True, help="species_refs.tsv")
    parser.add_argument("--mlst-organisms", required=True, help="JSON: MLST scheme -> AMRFinderPlus organism")
    parser.add_argument("--skani", help="skani search TSV (with --species_refs)")
    parser.add_argument("--taxa", help="Bracken species table or Kraken2 report")
    parser.add_argument("--taxa-format", choices=["bracken", "kraken2"], default="bracken")
    parser.add_argument("--expected", help="species from the samplesheet")
    parser.add_argument("--min-ani", type=float, default=95.0, help="minimum ANI for a species call (95)")
    parser.add_argument("--min-af", type=float, default=50.0, help="minimum aligned fraction of the sample (50)")
    parser.add_argument("--ani-margin", type=float, default=0.5, help="ANI margin to the next organism (0.5)")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    row = decide(args)
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    sys.exit(main())
