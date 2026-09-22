#!/usr/bin/env python3
"""Convert Proteinortho .proteinortho/.poff tables to the PhyloPrep standard."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
from orthogroup_table import parse_matrix_table, write_standard_table


def main(argv=None):
    parser = argparse.ArgumentParser(description="Convert Proteinortho group output into Orthogroup-by-taxon TSV for extract_orthofinder_orthogroups.py.")
    parser.add_argument("-i", "--input", required=True, metavar="PROTEINORTHO", help="Proteinortho tab-separated result table.")
    parser.add_argument("-o", "--output", required=True, metavar="TSV", help="Standard Orthogroup TSV.")
    args = parser.parse_args(argv)
    try:
        taxa, groups = parse_matrix_table(Path(args.input), metadata_columns=3, header_comment=True, group_prefix="POG")
        write_standard_table(Path(args.output), taxa, groups, inputs=[args.input])
    except (OSError, ValueError) as error:
        print(f"proteinortho_to_orthogroups: {error}", file=sys.stderr)
        return 2
    print(f"Proteinortho groups converted: {len(groups)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
