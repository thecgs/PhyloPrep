#!/usr/bin/env python3
"""Convert OMA Standalone OrthologousGroups.txt to the PhyloPrep TSV standard."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
from orthogroup_table import parse_prefixed_tab_groups, write_standard_table


def main(argv=None):
    parser = argparse.ArgumentParser(description="Convert OMA OrthologousGroups.txt 'GROUP: genome|gene ...' output into TSV for extract_orthofinder_orthogroups.py.")
    parser.add_argument("-i", "--input", required=True, metavar="GROUPS", help="OMA Output/OrthologousGroups.txt.")
    parser.add_argument("-o", "--output", required=True, metavar="TSV", help="Standard Orthogroup TSV.")
    parser.add_argument("--taxon-separator", default="|", help="Genome/gene separator in members (default: '|').")
    parser.add_argument("--strip-taxon-prefix", action="store_true", help="Write gene IDs without the leading 'genome<separator>' portion.")
    args = parser.parse_args(argv)
    try:
        taxa, groups = parse_prefixed_tab_groups(Path(args.input), separator=args.taxon_separator, strip_prefix=args.strip_taxon_prefix)
        write_standard_table(Path(args.output), taxa, groups, inputs=[args.input])
    except (OSError, ValueError) as error:
        print(f"oma_to_orthogroups: {error}", file=sys.stderr)
        return 2
    print(f"OMA groups converted: {len(groups)}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
