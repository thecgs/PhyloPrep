#!/usr/bin/env python3
"""Extract selected one-to-one or multi-copy OrthoFinder orthogroups as FASTAs."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from Bio import SeqIO

from msap_io import atomic_output, check_output_path

ORTHOFINDER_FASTA_SUFFIXES = (".fa", ".faa", ".fasta", ".fas", ".pep")


def non_ambiguous_length(sequence: str, seqtype: str) -> int:
    """Count unambiguous bases/residues without modifying the source sequence."""
    alphabet = "ACGTU" if seqtype != "prot" else "ACDEFGHIKLMNPQRSTVWY"
    return sum(char in alphabet for char in sequence.upper())


def orthofinder_safe_id(identifier: str) -> str:
    """Apply OrthoFinder's accession substitutions in result labels."""
    return identifier.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_")


def read_taxa_file(path: Path) -> list[str]:
    taxa = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            name = line.strip()
            if not name or name.startswith("#"):
                continue
            taxa.append(name)
    if not taxa:
        raise ValueError(f"Taxon list is empty: {path}")
    return taxa


def expand_taxa(values: list[str] | None) -> list[str] | None:
    if not values:
        return None
    taxa = []
    for value in values:
        path = Path(value)
        taxa.extend(read_taxa_file(path) if path.is_file() else [value])
    duplicates = sorted({taxon for taxon in taxa if taxa.count(taxon) > 1})
    if duplicates:
        raise ValueError("Taxa were specified more than once: " + ", ".join(duplicates))
    return taxa


def select_taxa(all_taxa: list[str], keep: list[str] | None, drop: list[str] | None) -> list[str]:
    named_taxa = (keep or []) + (drop or [])
    unknown = [taxon for taxon in named_taxa if taxon not in all_taxa]
    if unknown:
        raise ValueError("Taxa absent from Orthogroups.tsv: " + ", ".join(unknown))
    overlap = sorted(set(keep or []).intersection(drop or []))
    if overlap:
        raise ValueError("--keep-taxa and --drop-taxa select the same taxa: " + ", ".join(overlap))
    requested = [taxon for taxon in all_taxa if (keep is None or taxon in keep) and taxon not in (drop or [])]
    if len(requested) < 2:
        raise ValueError("At least two taxa are required after taxon filtering")
    return requested


def read_orthogroups(path: Path) -> tuple[list[str], list[tuple[str, dict[str, str]]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "Orthogroup":
            raise ValueError(f"{path} must be an OrthoFinder Orthogroups.tsv file with an 'Orthogroup' first column")
        taxa = reader.fieldnames[1:]
        if not taxa or len(set(taxa)) != len(taxa):
            raise ValueError(f"{path} has no taxa columns or contains duplicate taxon columns")
        rows, seen = [], set()
        for line_number, row in enumerate(reader, 2):
            orthogroup = (row.get("Orthogroup") or "").strip()
            if not orthogroup:
                raise ValueError(f"{path}:{line_number}: empty Orthogroup value")
            if orthogroup in seen:
                raise ValueError(f"{path}:{line_number}: duplicate Orthogroup value: {orthogroup}")
            seen.add(orthogroup)
            rows.append((orthogroup, {taxon: (row.get(taxon) or "").strip() for taxon in taxa}))
    return taxa, rows


def genes_in_cell(cell: str, orthogroup: str, taxon: str) -> list[str]:
    """Read an OrthoFinder cell as zero, one, or several comma-separated IDs."""
    if not cell:
        return []
    genes = [gene.strip() for gene in cell.split(",")]
    if any(not gene for gene in genes):
        raise ValueError(f"{orthogroup}: malformed gene list for taxon '{taxon}'")
    return genes


def detect_orthofinder_id_style(rows: list[tuple[str, dict[str, str]]], taxa: list[str],
                                requested: str = "auto") -> str:
    """Detect whether table IDs are normal OrthoFinder or OrthoFinder ``-X`` IDs."""
    if requested != "auto":
        return requested
    styles = set()
    for orthogroup, cells in rows:
        for taxon in taxa:
            prefix = orthofinder_safe_id(taxon) + "_"
            for identifier in genes_in_cell(cells[taxon], orthogroup, taxon):
                styles.add("prefixed" if orthofinder_safe_id(identifier).startswith(prefix) else "unprefixed")
    if len(styles) > 1:
        raise ValueError(
            "Orthogroups.tsv mixes taxon-prefixed and unprefixed gene IDs. "
            "Use a consistent OrthoFinder result, or set --orthofinder-id-style explicitly "
            "when the input is a nonstandard converted table."
        )
    return next(iter(styles), "unprefixed")


def table_gene_id(identifier: str, taxon: str, id_style: str) -> tuple[str, str]:
    """Return normalized table and source-FASTA IDs for the selected ID style."""
    table_id = orthofinder_safe_id(identifier)
    if id_style == "unprefixed":
        return table_id, table_id
    prefix = orthofinder_safe_id(taxon) + "_"
    if not table_id.startswith(prefix):
        raise ValueError(
            f"{taxon}: expected a '{prefix}geneID' value for --orthofinder-id-style prefixed, "
            f"found '{identifier}'"
        )
    return table_id, table_id[len(prefix):]


def occupied_taxa(cells: dict[str, str], taxa: list[str], orthogroup: str, copy_mode: str) -> list[str]:
    """Return taxa eligible under the requested single- or multi-copy mode."""
    result = []
    for taxon in taxa:
        genes = genes_in_cell(cells[taxon], orthogroup, taxon)
        if (copy_mode == "single-copy" and len(genes) == 1) or (copy_mode == "all-copies" and genes):
            result.append(taxon)
    return result


def fasta_files_by_taxon(sequence_directory: Path) -> dict[str, Path]:
    if not sequence_directory.is_dir():
        raise ValueError(f"Sequence directory does not exist or is not a directory: {sequence_directory}")
    files = {}
    for path in sorted(sequence_directory.iterdir()):
        if not path.is_file() or path.suffix.lower() not in ORTHOFINDER_FASTA_SUFFIXES:
            continue
        taxon = path.name.rsplit(".", 1)[0]
        if taxon in files:
            raise ValueError(f"More than one OrthoFinder-style FASTA file has taxon name '{taxon}': {files[taxon]}, {path}")
        files[taxon] = path
    if not files:
        raise ValueError(f"No FASTA files found in {sequence_directory} (accepted suffixes: {', '.join(ORTHOFINDER_FASTA_SUFFIXES)})")
    return files


def load_sequences(path: Path, wanted_ids: set[str] | None = None) -> dict[str, str]:
    """Scan a taxon FASTA once, retaining sequence contents only for requested IDs."""
    sequences = {}
    seen = set()
    for record in SeqIO.parse(path, "fasta"):
        identifier = orthofinder_safe_id(record.id)
        if not identifier:
            raise ValueError(f"{path}: empty FASTA ID")
        if identifier in seen:
            raise ValueError(f"{path}: duplicate OrthoFinder-normalized sequence ID: {identifier}")
        seen.add(identifier)
        if wanted_ids is None or identifier in wanted_ids:
            sequences[identifier] = str(record.seq)
    if not seen:
        raise ValueError(f"No sequences found in FASTA: {path}")
    return sequences


def extract_orthologs(orthogroups: Path, sequence_directory: Path, outdir: Path,
                      keep_taxa: list[str] | None = None, drop_taxa: list[str] | None = None,
                      min_taxa: int | None = None, no_species_prefix: bool = False,
                      map_output: Path | None = None, copy_mode: str = "single-copy",
                      pathlist_output: Path | None = None, orthofinder_id_style: str = "auto",
                      longest_copy: bool = False, seqtype: str = "codon") -> int:
    if seqtype not in ("codon", "prot", "nucl", "pseudogene"):
        raise ValueError("seqtype must be codon, prot, nucl or pseudogene")
    if longest_copy and copy_mode != "all-copies":
        raise ValueError("--longest-copy requires --copy-mode all-copies")
    all_taxa, rows = read_orthogroups(orthogroups)
    taxa = select_taxa(all_taxa, keep_taxa, drop_taxa)
    id_style = detect_orthofinder_id_style(rows, taxa, orthofinder_id_style)
    if min_taxa is None:
        min_taxa = len(taxa)
    if not 2 <= min_taxa <= len(taxa):
        raise ValueError(f"--min-taxa must be between 2 and {len(taxa)} after taxon filtering")
    sequence_files = fasta_files_by_taxon(sequence_directory)
    missing_files = [taxon for taxon in taxa if taxon not in sequence_files]
    if missing_files:
        raise ValueError("Sequence FASTA is missing for selected taxa: " + ", ".join(missing_files))
    selected = [(og, cells, occupied_taxa(cells, taxa, og, copy_mode)) for og, cells in rows]
    selected = [(og, cells, present_taxa) for og, cells, present_taxa in selected if len(present_taxa) >= min_taxa]
    if not selected:
        raise ValueError("No orthogroups meet the --min-taxa threshold after taxon filtering")
    if copy_mode == "single-copy" and min_taxa < len(taxa):
        incomplete = sum(len(present) < len(taxa) for _, _, present in selected)
        print(
            f"WARNING: --copy-mode single-copy --min-taxa {min_taxa} permits missing taxa among "
            f"{len(taxa)} selected taxa ({incomplete} retained OGs are incomplete). Multi-copy taxa are omitted, "
            "not resolved by choosing a copy. To retain incomplete loci during downstream concatenation, "
            "use phyloprep.py --missing-taxa pad-gaps (or get_supergenes.py --missing-taxa pad-gaps); "
            "the default skip-gene policy excludes loci missing a target taxon.", file=sys.stderr)
    if longest_copy:
        print(
            "WARNING: --copy-mode all-copies --longest-copy selects representatives by non-ambiguous length, "
            "not by inferred orthology. Different species may retain different paralogs; the longest "
            "sequence can also be chimeric or a fusion. This does not establish single-copy orthologs. "
            "Inspect copy_selection.tsv and use independent orthology evidence before concatenation.",
            file=sys.stderr)
    # Filter OGs before retaining sequence contents. All copies in a retained
    # cell remain candidates for --longest-copy; none may be skipped by length.
    needed_ids = {taxon: set() for taxon in taxa}
    for orthogroup, cells, present_taxa in selected:
        for taxon in present_taxa:
            for table_id in genes_in_cell(cells[taxon], orthogroup, taxon):
                _, gene_id = table_gene_id(table_id, taxon, id_style)
                needed_ids[taxon].add(gene_id)
    sequences_by_taxon = {
        taxon: load_sequences(sequence_files[taxon], needed_ids[taxon]) for taxon in taxa
    }
    outdir.mkdir(parents=True, exist_ok=True)
    if not outdir.is_dir():
        raise ValueError(f"Output path is not a directory: {outdir}")
    map_output = map_output or outdir / "map.tsv"
    pathlist_output = pathlist_output or outdir / "orthogroups.pathlist"
    inputs = [orthogroups, *[sequence_files[taxon] for taxon in taxa]]
    check_output_path(map_output, inputs=inputs)
    check_output_path(pathlist_output, inputs=inputs + [map_output])
    selection_output = outdir / "copy_selection.tsv" if longest_copy else None
    if selection_output:
        check_output_path(selection_output, inputs=inputs + [map_output, pathlist_output])
    prepared, map_rows, selection_rows = [], [], []
    mapped_taxa = {}
    # Validate every orthogroup before publishing any output, so a missing sequence
    # cannot leave a partially updated output directory.
    for orthogroup, cells, present_taxa in selected:
        output = outdir / f"{orthogroup}.fasta"
        check_output_path(output, inputs=inputs + [map_output, pathlist_output] + ([selection_output] if selection_output else []))
        output_rows, seen_ids = [], set()
        for taxon in present_taxa:
            candidates = []
            for table_id in genes_in_cell(cells[taxon], orthogroup, taxon):
                _normalized_table_id, gene_id = table_gene_id(table_id, taxon, id_style)
                if gene_id not in sequences_by_taxon[taxon]:
                    raise ValueError(
                        f"{orthogroup}: sequence ID '{gene_id}' for taxon '{taxon}' was not found in {sequence_files[taxon]}. "
                        "Check that the input FASTAs and Orthogroups.tsv use identical OrthoFinder-normalized IDs."
                    )
                candidates.append(gene_id)
            if longest_copy:
                if len(set(candidates)) != len(candidates):
                    raise ValueError(f"{orthogroup}: duplicate candidate gene IDs for taxon '{taxon}'")
                effective_lengths = {
                    gene: non_ambiguous_length(sequences_by_taxon[taxon][gene], seqtype)
                    for gene in candidates
                }
                chosen = min(candidates, key=lambda gene: (-effective_lengths[gene], gene))
                selection_rows.extend((orthogroup, taxon, gene, len(sequences_by_taxon[taxon][gene]),
                                       effective_lengths[gene], seqtype,
                                       'yes' if gene == chosen else 'no') for gene in candidates)
                candidates = [chosen]
            for gene_id in candidates:
                export_id = gene_id if no_species_prefix else orthofinder_safe_id(taxon) + "_" + gene_id
                if any(char.isspace() for char in export_id):
                    raise ValueError(f"{orthogroup}: output ID contains whitespace: {export_id!r}")
                if export_id in seen_ids:
                    raise ValueError(
                        f"{orthogroup}: output IDs are not unique after {'-X' if no_species_prefix else 'species-prefix'} processing: "
                        f"'{export_id}'. Use globally unique gene IDs or omit -X."
                    )
                seen_ids.add(export_id)
                output_rows.append((export_id, gene_id, taxon, sequences_by_taxon[taxon][gene_id]))
                if export_id in mapped_taxa and mapped_taxa[export_id] != taxon:
                    raise ValueError(
                        f"Output ID '{export_id}' maps to multiple taxa across orthogroups: "
                        f"{mapped_taxa[export_id]} and {taxon}. Omit -X or use globally unique gene IDs."
                    )
                mapped_taxa[export_id] = taxon
                map_rows.append((export_id, taxon))
        prepared.append((output, orthogroup, output_rows))
    output_fastas = [output for output, _orthogroup, _output_rows in prepared]
    check_output_path(pathlist_output, inputs=inputs + [map_output, *output_fastas])
    for output, orthogroup, output_rows in prepared:
        with atomic_output(output, inputs=inputs + [map_output, pathlist_output]) as handle:
            for export_id, _gene_id, taxon, sequence in output_rows:
                handle.write(f">{export_id} taxon={taxon} orthogroup={orthogroup}\n{sequence}\n")
    with atomic_output(map_output, inputs=inputs + output_fastas) as handle:
        for fasta_id, taxon in map_rows:
            handle.write(f"{fasta_id}\t{taxon}\n")
    with atomic_output(pathlist_output, inputs=inputs + [map_output, *output_fastas]) as handle:
        for output in output_fastas:
            handle.write(str(output.resolve()) + "\n")
    if selection_output:
        with atomic_output(selection_output, inputs=inputs + output_fastas + [map_output, pathlist_output]) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["orthogroup", "taxon", "gene_id", "length", "effective_length", "sequence_type", "selected"])
            writer.writerows(selection_rows)
    return len(prepared)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Extract selected single-copy or multi-copy sequence orthogroups from OrthoFinder Orthogroups.tsv.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  extract_orthofinder_orthogroups.py -og Results/Orthogroups/Orthogroups.tsv -i cds/ -o one_to_one
  extract_orthofinder_orthogroups.py -og Orthogroups.tsv -i cds/ --keep-taxa Human Mouse Zebrafish
  extract_orthofinder_orthogroups.py -og Orthogroups.tsv -i cds/ --drop-taxa Fruit_fly Frog
  extract_orthofinder_orthogroups.py -og Orthogroups.tsv -i cds/ --min-taxa 8 -X
  extract_orthofinder_orthogroups.py -og Orthogroups.tsv -i cds/ --copy-mode all-copies --min-taxa 8
  extract_orthofinder_orthogroups.py -og Orthogroups.tsv -i cds/ --copy-mode all-copies --longest-copy

The input sequence directory follows OrthoFinder -f conventions: one FASTA per taxon;
the filename stem must exactly equal the taxon column in Orthogroups.tsv. Accepted
suffixes are .fa, .faa, .fasta, .fas and .pep. Taxon options accept names or a
one-name-per-line text file. Default: retain all taxa in Orthogroups.tsv.

--keep-taxa and --drop-taxa may be supplied together, but overlapping taxa are
an error. --copy-mode single-copy (default) counts a taxon only when its cell
has exactly one gene and writes that gene. --copy-mode all-copies counts a taxon
when it has at least one gene and writes every listed copy; this mode preserves
paralogs for gene-tree methods such as ASTRAL-Pro3. With --longest-copy, retain
the sequence with the most non-ambiguous characters within each OG and taxon.
-st codon (default), nucl and pseudogene count A/C/G/T/U; prot counts the 20 standard amino
acids, including N. Counting is case-insensitive; ambiguity, gap, stop and other
symbols do not contribute. Effective-length ties choose the lexicographically
smallest normalized gene ID, without using raw length as a secondary criterion.
No sequence cleaning is performed. copy_selection.tsv records raw length,
effective_length, sequence_type and selection status for all candidates.
This is representative selection, not orthology inference.
--min-taxa still counts occupied taxa, not copies, and is evaluated
under the chosen mode. Its maximum is the number of taxa remaining after taxon
filtering; its default is that maximum. The default --orthofinder-id-style auto
detects taxon-prefixed table IDs (normal OrthoFinder output) or unprefixed IDs
(OrthoFinder -X output). Detection controls source FASTA lookup only.
Default exported IDs are always Taxon_geneID; -X always emits geneID without
adding the taxon prefix. An existing table prefix is stripped before the
output policy is applied, so normal prefixed tables are not double-prefixed. Use
--orthofinder-id-style prefixed|unprefixed only to override auto detection for
a nonstandard converted table. map.tsv contains one
'FASTA_ID<TAB>Taxa' row for every written sequence; its first column always
matches the exported FASTA ID. orthogroups.pathlist contains
absolute output FASTA paths and can be passed directly to phyloprep.py -i.
""")
    parser.add_argument("-og", "--orthogroups", required=True, metavar="TSV", help="OrthoFinder Orthogroups.tsv file.")
    parser.add_argument("-i", "-f", "--sequence-dir", dest="sequence_dir", required=True, metavar="DIR",
                        help="Directory of input FASTA files, one file per taxon (like OrthoFinder -f).")
    parser.add_argument("-k", "--keep-taxa", nargs="+", metavar="TAXON_OR_TXT",
                        help="Retain only these taxa; overlapping --drop-taxa entries are errors.")
    parser.add_argument('-d',"--drop-taxa", nargs="+", metavar="TAXON_OR_TXT",
                        help="Exclude these taxa; overlapping --keep-taxa entries are errors.")
    parser.add_argument('-m',"--copy-mode", choices=("single-copy", "all-copies"), default="single-copy",
                        help="single-copy: require exactly one gene per counted taxon (default); all-copies: retain every copy and count taxa with >=1 gene, for e.g. ASTRAL-Pro3.")
    parser.add_argument("--longest-copy", action="store_true",
                        help="all-copies only: retain the sequence with most non-ambiguous characters per taxon within each OG; ties use gene ID order. This does not establish orthology.")
    parser.add_argument("-st", "--seqtype", "--seq-type", choices=("codon", "prot", "nucl", "pseudogene"), default="codon",
                        help="Alphabet for --longest-copy scoring: codon (default), nucl and pseudogene count A/C/G/T/U; prot counts the 20 standard amino acids. No translation or codon grouping is performed. Does not modify sequences.")
    parser.add_argument("--min-taxa", type=int, metavar="N",
                        help="Require N occupied taxa after filtering; occupancy follows --copy-mode (default: all selected taxa).")
    parser.add_argument("-X", action="store_true", dest="no_species_prefix",
                        help="Write gene IDs without adding the taxon prefix. Without -X, always write Taxon_geneID, regardless of input table style.")
    parser.add_argument("--orthofinder-id-style", choices=("auto", "prefixed", "unprefixed"), default="auto",
                        help="Orthogroups.tsv ID style: auto-detect (default), prefixed (normal OrthoFinder), or unprefixed (OrthoFinder -X).")
    parser.add_argument("-o", "--outdir", default="one-to-one_orthologs", metavar="DIR", help="Output directory.")
    parser.add_argument('-mo',"--map-output", metavar="TSV", help="Exported-FASTA-ID-to-taxon TSV; default is OUTDIR/map.tsv.")
    parser.add_argument("--pathlist-output", metavar="TXT",
                        help="Exported OG FASTA path list; default is OUTDIR/orthogroups.pathlist.")
    parser.add_argument("-v", "--version", action="version", version="v1.0.0")
    args = parser.parse_args(argv)
    try:
        count = extract_orthologs(Path(args.orthogroups), Path(args.sequence_dir), Path(args.outdir),
                                  expand_taxa(args.keep_taxa), expand_taxa(args.drop_taxa),
                                  args.min_taxa, args.no_species_prefix,
                                  Path(args.map_output) if args.map_output else None,
                                  args.copy_mode,
                                  Path(args.pathlist_output) if args.pathlist_output else None,
                                  args.orthofinder_id_style, args.longest_copy, args.seqtype)
    except (OSError, ValueError) as error:
        print(f"extract_orthofinder_orthogroups: {error}", file=sys.stderr)
        return 2
    print(f"orthogroups written: {count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
