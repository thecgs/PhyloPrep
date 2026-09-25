#!/bin/sh
# Build PhyloPrep from any working directory.  Singularity resolves %files
# paths from its current directory, so always invoke it from this repository.
set -eu

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 OUTPUT.sif" >&2
    exit 2
fi

caller_dir=$(pwd -P)
case $1 in
    /*) output=$1 ;;
    *) output=$caller_dir/$1 ;;
esac

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
cd "$repo_dir"
exec singularity build "$output" "$repo_dir/Singularity.def"
