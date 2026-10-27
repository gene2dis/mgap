#!/usr/bin/env python3
"""Plasmid calls per contig for one MGAP sample (roadmap D-22).

Two subcommands, both reading the sample's annotation FASTA (contigs renamed contig_1..N, the
assembler's header kept as the description):

circular
    Writes the circular contigs to a FASTA for MOB-typer. A contig is circular when Flye's
    assembly_info.txt marks it (circ. = Y; matched on the first word of the description, the
    Flye contig name) or when its header carries circular=true or topology=circular (Autocycler,
    Unicycler, Bakta). With --all-if-untagged (contig mode) a sample whose headers carry no
    topology tag at all gets every contig. Contigs of --max-length bp or more (chromosomes) are
    left out. Nothing is written when no contig is selected.

consensus
    One row per contig with the evidence from each plasmid tool and a consensus call. Three
    sources vote: geNomad (contig in the plasmid summary), MOB-suite (mob_recon molecule_type
    plasmid, or a replicon or relaxase from MOB-typer on a contig known to be circular; on other
    contigs, such as every contig of an untagged contig-mode FASTA, MOB-typer adds evidence
    columns but no vote) and PlasmidFinder
    (a replicon hit). Calls:
        plasmid     2 or 3 sources
        candidate   1 source
        check       any source, on a contig of --large-contig bp or more (chromosome-sized)
        chromosome  no source, and mob_recon calls the contig chromosome
        no call     no source, and MOB-suite did not run (nothing calls it chromosome)
    Circularity (Flye, header tag, or a geNomad DTR topology) and mobility (MOB-typer for the
    contig, else mob_recon's plasmid cluster) are evidence columns, not votes. A tool that did not
    run is "not run" in its columns.
"""

import argparse
import csv
import re
import sys

NOT_RUN = "not run"
EMPTY = "-"
LARGE_CONTIG = 1_000_000

# circular=true (Autocycler, Unicycler, Plassembler), topology=circular and [topology=circular] (Bakta)
CIRCULAR_TAG = re.compile(r"(?:^|[\s\[;,])(?:circular=(?:true|yes|y)|topology=circular)(?=$|[\s\];,])", re.I)
LINEAR_TAG = re.compile(r"(?:^|[\s\[;,])(?:circular=(?:false|no|n)|topology=linear)(?=$|[\s\];,])", re.I)

GENOMAD_COLUMNS = ["genomad_plasmid", "genomad_score", "genomad_topology", "genomad_conjugation_genes"]
MOBSUITE_COLUMNS = [
    "mobsuite_molecule_type",
    "mobsuite_cluster",
    "mobsuite_replicons",
    "mobsuite_relaxases",
    "mobsuite_mpf",
    "mobsuite_orit",
    "mob_typer",
    "predicted_mobility",
    "mobility_source",
    "mash_nearest_neighbor",
    "mash_neighbor_identification",
]
COLUMNS = [
    "sample",
    "contig",
    "description",
    "length",
    "circular",
    "circular_source",
    "consensus",
    "support",
    "evidence",
    *GENOMAD_COLUMNS,
    *MOBSUITE_COLUMNS,
    "plasmidfinder_replicons",
    "plasmidfinder_identity",
]


def first_word(text):
    words = (text or "").split()
    return words[0] if words else ""


def read_fasta(path):
    """[(contig, description, length, sequence lines)] in file order."""
    records = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:]
                contig = first_word(header)
                description = header[len(contig) :].strip()
                records.append([contig, description, 0, []])
            elif records and line:
                records[-1][2] += len(line.strip())
                records[-1][3].append(line)
    return records


def read_table(path):
    """Rows of a tab-separated file with a header, or [] for an empty file."""
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def flye_circular(path):
    """Flye contig names marked circular in assembly_info.txt."""
    circular, linear = set(), set()
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) < 4:
                continue
            (circular if fields[3].strip().upper() == "Y" else linear).add(fields[0].strip())
    return circular, linear


def topology(records, flye_info=None):
    """{contig: (circular yes/no/-, source)} from Flye's assembly_info or the header tags."""
    flye = flye_circular(flye_info) if flye_info else None
    result = {}
    for contig, description, _length, _seq in records:
        if flye is not None:
            original = first_word(description) or contig
            if original in flye[0]:
                result[contig] = ("yes", "Flye")
            elif original in flye[1]:
                result[contig] = ("no", "Flye")
            else:
                result[contig] = (EMPTY, "")
        elif CIRCULAR_TAG.search(description):
            result[contig] = ("yes", "header")
        elif LINEAR_TAG.search(description):
            result[contig] = ("no", "header")
        else:
            result[contig] = (EMPTY, "")
    return result


def select_circular(args):
    records = read_fasta(args.fasta)
    shape = topology(records, args.flye_info)
    candidates = [record for record in records if record[2] < args.max_length]
    selected = [record for record in candidates if shape[record[0]][0] == "yes"]
    reason = f"circular, shorter than {args.max_length:,} bp"
    tagged = any(source for _status, source in shape.values())
    if not selected and args.all_if_untagged and not tagged:
        selected, reason = (
            candidates,
            f"no topology tag in the headers, every contig shorter than {args.max_length:,} bp",
        )
    print(f"{len(selected)} of {len(records)} contigs selected ({reason})", file=sys.stderr)
    if not selected:
        return
    with open(args.out, "w") as handle:
        for contig, description, _length, sequence in selected:
            handle.write(f">{contig} {description}".rstrip() + "\n")
            handle.write("\n".join(sequence) + "\n")


def value(row, *names):
    """The first of the named columns that has a value, '-' otherwise."""
    for name in names:
        text = (row.get(name) or "").strip()
        if text and text != "-":
            return text
    return EMPTY


def join_unique(*values):
    items = []
    for text in values:
        for item in (text or "").split(","):
            item = item.strip()
            if item and item != EMPTY and item not in items:
                items.append(item)
    return ",".join(items) or EMPTY


def by_contig(rows, column):
    table = {}
    for row in rows:
        table.setdefault(first_word(row.get(column, "")).split("|")[0], row)
    return table


def plasmidfinder_hits(path):
    hits = {}
    for row in read_table(path):
        contig = first_word(row.get("Contig", ""))
        if contig:
            hits.setdefault(contig, []).append((row.get("Plasmid", "").strip(), row.get("Identity", "").strip()))
    return hits


def recon_mobility(path):
    """{mob_recon cluster: predicted mobility} from mobtyper_results.txt."""
    mobility = {}
    for row in read_table(path):
        cluster = value(row, "primary_cluster_id")
        if cluster == EMPTY:
            cluster = (row.get("sample_id") or "").rsplit(":", 1)[-1]
        mobility[cluster] = value(row, "predicted_mobility")
    return mobility


def consensus(args):
    records = read_fasta(args.fasta)
    shape = topology(records, args.flye_info)
    genomad = by_contig(read_table(args.genomad), "seq_name") if args.genomad else None
    recon = by_contig(read_table(args.mob_contigs), "contig_id") if args.mob_contigs else None
    recon_mob = recon_mobility(args.mob_plasmids) if args.mob_plasmids else {}
    typer = by_contig(read_table(args.mob_typer), "sample_id") if args.mob_typer else {}
    finder = plasmidfinder_hits(args.plasmidfinder) if args.plasmidfinder else None
    mobsuite_run = recon is not None or bool(args.mob_typer)

    rows = []
    for contig, description, length, _seq in records:
        circular, circular_source = shape[contig]
        row = {"sample": args.sample, "contig": contig, "description": description or EMPTY, "length": length}
        votes = []

        if genomad is None:
            row.update(dict.fromkeys(GENOMAD_COLUMNS, NOT_RUN))
        else:
            hit = genomad.get(contig)
            row["genomad_plasmid"] = "yes" if hit else "no"
            row["genomad_score"] = value(hit, "plasmid_score") if hit else EMPTY
            row["genomad_topology"] = value(hit, "topology") if hit else EMPTY
            row["genomad_conjugation_genes"] = value(hit, "conjugation_genes") if hit else EMPTY
            if hit:
                votes.append("geNomad")
                # direct terminal repeats: geNomad found the contig's ends overlapping, so it is circular
                if circular == EMPTY and row["genomad_topology"] == "DTR":
                    circular, circular_source = "yes", "geNomad DTR"
        row["circular"], row["circular_source"] = circular, circular_source or EMPTY

        mob_vote = False
        if not mobsuite_run:
            row.update(dict.fromkeys(MOBSUITE_COLUMNS, NOT_RUN))
        else:
            contig_row = (recon or {}).get(contig, {})
            typed = typer.get(contig)
            molecule = value(contig_row, "molecule_type")
            cluster = value(contig_row, "primary_cluster_id")
            row["mobsuite_molecule_type"] = molecule if recon is not None else NOT_RUN
            row["mobsuite_cluster"] = cluster if recon is not None else NOT_RUN
            row["mobsuite_replicons"] = join_unique(value(contig_row, "rep_type(s)"), value(typed or {}, "rep_type(s)"))
            row["mobsuite_relaxases"] = join_unique(
                value(contig_row, "relaxase_type(s)"), value(typed or {}, "relaxase_type(s)")
            )
            row["mobsuite_mpf"] = join_unique(value(contig_row, "mpf_type"), value(typed or {}, "mpf_type"))
            row["mobsuite_orit"] = join_unique(value(contig_row, "orit_type(s)"), value(typed or {}, "orit_type(s)"))
            row["mob_typer"] = ("yes" if typed else "no") if args.mob_typer else NOT_RUN
            # MOB-typer's prediction is for the whole molecule only on a circular contig; on a fragment of a
            # mob_recon plasmid, the cluster's prediction stands
            if typed and (circular == "yes" or molecule != "plasmid"):
                row["predicted_mobility"], row["mobility_source"] = value(typed, "predicted_mobility"), "MOB-typer"
                row["mash_nearest_neighbor"] = value(typed, "mash_nearest_neighbor")
                row["mash_neighbor_identification"] = value(typed, "mash_neighbor_identification")
            else:
                mobility = recon_mob.get(cluster, EMPTY) if molecule == "plasmid" else EMPTY
                row["predicted_mobility"] = mobility
                row["mobility_source"] = "mob_recon cluster" if mobility != EMPTY else EMPTY
                row["mash_nearest_neighbor"] = value(contig_row, "mash_nearest_neighbor")
                row["mash_neighbor_identification"] = value(contig_row, "mash_neighbor_identification")
            # MOB-typer's markers vote only on a contig known to be circular (a whole molecule); on a
            # fragment of a draft, mob_recon's molecule_type decides
            typed_markers = (
                typed
                and circular == "yes"
                and (value(typed, "rep_type(s)") != EMPTY or value(typed, "relaxase_type(s)") != EMPTY)
            )
            mob_vote = molecule == "plasmid" or bool(typed_markers)
            if mob_vote:
                votes.append("MOB-suite")

        if finder is None:
            row["plasmidfinder_replicons"] = row["plasmidfinder_identity"] = NOT_RUN
        else:
            hits = finder.get(contig, [])
            row["plasmidfinder_replicons"] = ",".join(name for name, _identity in hits) or EMPTY
            row["plasmidfinder_identity"] = ",".join(identity for _name, identity in hits) or EMPTY
            if hits:
                votes.append("PlasmidFinder")

        row["support"] = len(votes)
        if votes and length >= args.large_contig:
            row["consensus"] = "check"
            row["evidence"] = ", ".join(votes) + f" (contig of {length:,} bp)"
        elif len(votes) >= 2:
            row["consensus"], row["evidence"] = "plasmid", ", ".join(votes)
        elif votes:
            row["consensus"], row["evidence"] = "candidate", f"{votes[0]} only"
        elif mobsuite_run and row["mobsuite_molecule_type"] == "chromosome":
            row["consensus"], row["evidence"] = "chromosome", "MOB-suite: chromosome"
        else:
            row["consensus"] = "no call"
            row["evidence"] = "no plasmid evidence" if (genomad is not None or finder is not None) else EMPTY
        rows.append(row)

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    calls = {}
    for row in rows:
        calls[row["consensus"]] = calls.get(row["consensus"], 0) + 1
    print(f"{args.sample}: {len(rows)} contigs, " + ", ".join(f"{n} {c}" for c, n in calls.items()), file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    circular = commands.add_parser("circular", help="write the circular contigs for MOB-typer")
    circular.add_argument("--fasta", required=True, help="annotation FASTA")
    circular.add_argument("--flye-info", help="Flye assembly_info.txt")
    circular.add_argument(
        "--all-if-untagged", action="store_true", help="every contig when no header carries a topology tag"
    )
    circular.add_argument(
        "--max-length",
        type=int,
        default=LARGE_CONTIG,
        help=f"leave out contigs this long or longer: chromosomes (default {LARGE_CONTIG})",
    )
    circular.add_argument("--out", required=True, help="FASTA to write (not written when nothing is selected)")
    circular.set_defaults(func=select_circular)

    merge = commands.add_parser("consensus", help="one row per contig with every tool's evidence and a call")
    merge.add_argument("--sample", required=True)
    merge.add_argument("--fasta", required=True, help="annotation FASTA")
    merge.add_argument("--flye-info", help="Flye assembly_info.txt")
    merge.add_argument("--genomad", help="geNomad <sample>_plasmid_summary.tsv")
    merge.add_argument("--mob-contigs", help="mob_recon contig_report.txt")
    merge.add_argument("--mob-plasmids", help="mob_recon mobtyper_results.txt")
    merge.add_argument("--mob-typer", help="MOB-typer --multi report on the circular contigs")
    merge.add_argument("--plasmidfinder", help="PlasmidFinder results_tab.tsv")
    merge.add_argument(
        "--large-contig",
        type=int,
        default=LARGE_CONTIG,
        help=f"contigs this long or longer get 'check' instead of a plasmid call (default {LARGE_CONTIG})",
    )
    merge.add_argument("--out", required=True, help="TSV to write")
    merge.set_defaults(func=consensus)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
