#!/usr/bin/env python
# coding: utf-8

import os
import sys
import glob
#import gzip
import argparse
import subprocess
import shutil
import tempfile
from Bio import SeqIO, Seq
from Bio.Data import CodonTable
from msap_io import atomic_output, read_fasta_alignment, normalize_dna
from sequence_audit import (preprocess_codon, preprocess_nucl, rejected_input, add_macse_options,
                            normalize_workflow, macse_align, preprocess_protein, add_protein_options,
                            merge_reports, EmptyCodonLocusError)

def remove_stop_codon(infile, outfile, genetic_code):
    preprocess_codon(infile, outfile, genetic_code, get_prefix(infile) + '.sequence_changes.tsv')


def translate_seq(infile, outfile, genetic_code):
    with atomic_output(outfile, inputs=[infile]) as out:
        for record in SeqIO.parse(infile, 'fasta'):
            record.seq = Seq.Seq(normalize_dna(record.seq))
            start_codon = str(record.seq[:3]).upper()
            if start_codon in Seq.CodonTable.unambiguous_dna_by_id[genetic_code].start_codons:
                protein_sequence = list(record.seq.translate(table=genetic_code, cds=False))
                protein_sequence[0] = "M"
                protein_sequence = ''.join(protein_sequence)
                print(f">{record.id}\n{protein_sequence}", file=out)
            else:
                print(f">{record.id}\n{record.seq.translate(table=genetic_code)}", file=out)
    return None

def get_prefix(infile):
    return os.path.splitext(os.path.basename(infile))[0]

def detect_seq_type(infile):
    seqence = ""
    lst = ["V", "L", "I", "E", "Q", "D", "N",
           "M", "S", "F", "W", "Y", "R", "H", 
           "P","K","X"]
    for n, record in enumerate(SeqIO.parse(infile, 'fasta')):
        if n > 1:
            break
        else:
            seqence += str(record.seq.upper())

    status = False
    for s in seqence:
        if s in lst:
            status = False
        else:
            if status==True:
                break
    
    if status==True:
        seq_type="PROTEIN"
    else:
        seq_type="DNA"
    return seq_type

def check_dependencies(softwares):
    #print("Check dependencies...\n")
    main_path = os.path.dirname(os.path.abspath(__file__))
    
    paths = os.environ.get('PATH', os.defpath).split(os.pathsep)
    software_path = {software:None for software in softwares}
    for i in glob.glob(os.path.join(main_path, '*')):
        if os.path.isdir(i):
            paths.append(i)
    for i in glob.glob(os.path.join(main_path, '*', 'bin')):
        if os.path.isdir(i):
            paths.append(i)
    for i in glob.glob(os.path.join(main_path, '*', 'script')):
        if os.path.isdir(i):
            paths.append(i)
            
    for software in softwares:
        for path in paths:
            candidate = os.path.join(path, software)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                software_path[software] = os.path.abspath(candidate)
                break
    
    for software in software_path:
        if software_path[software] == None:
            print(software, "is not installed in the environment. Please install it.", file=sys.stderr)
            sys.exit(1)
        #else:
        #    print(software, ":", software_path[software])
    return software_path


def validate_input_fasta(infile, seqtype, protein_stop_symbol='*'):
    """Validate FASTA records before invoking an external alignment program."""
    if not os.path.isfile(infile):
        raise ValueError(f"Input FASTA file does not exist: {infile}")

    records = list(SeqIO.parse(infile, "fasta"))
    if len(records) < 2:
        raise ValueError("Input FASTA must contain at least two sequences.")

    ids = set()
    for record in records:
        if record.id in ids:
            raise ValueError(f"Duplicate sequence ID: {record.id}")
        ids.add(record.id)
        if not record.seq:
            raise ValueError(f"{record.id}: sequence is empty")

    if seqtype in {"nucl", "codon", "pseudogene"}:
        allowed_characters = set("ACGTURYSWKMBDHVN?-")
    else:
        allowed_characters = set("ACDEFGHIKLMNPQRSTVWYBXZJUO?-") | {protein_stop_symbol}

    for record in records:
        invalid_characters = set(str(record.seq).upper()) - allowed_characters
        if invalid_characters:
            invalid = "".join(sorted(invalid_characters))
            raise ValueError(f"{record.id}: invalid {seqtype} character(s): {invalid}")

        if seqtype == "pseudogene" and ("-" in record.seq or "U" in record.seq.upper()):
            raise ValueError(f"{record.id}: pseudogene input must be ungapped DNA (use T, not U)")
        if seqtype == "codon":
            if "-" in record.seq:
                raise ValueError(f"{record.id}: codon input must not contain gaps")


def trim_alignment(
    trim_software,
    trimal_path,
    trimal_args,
    infile,
    outfile,
    seqtype,
    gap_ratio,
    n_ratio,
    x_ratio,
):
    """Trim an alignment with the configured trimming program."""
    if trim_software == "trimAlnSeq":
        command = [
            sys.executable,
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "trimAlnSeq.py"),
            "-i", infile,
            "-o", outfile,
            "-st", seqtype,
            "-G", str(gap_ratio),
        ]
        if seqtype == "prot":
            command.extend(["-X", str(x_ratio)])
        else:
            command.extend(["-N", str(n_ratio)])
    else:
        print(
            "Note: trimAl does not use MSAP's -G, -N, or -X options. "
            "Use --trimal-args to pass trimAl options.",
            file=sys.stderr,
        )
        command = [
            trimal_path,
            "-in", infile,
            "-out", outfile,
            "-fasta",
        ]
        if trimal_args:
            command.extend(trimal_args)
        else:
            command.append("-automated1")

    subprocess.run(command, check=True)
    read_fasta_alignment(outfile, codon=(seqtype == "codon"))


def trim_codon_with_trimal(trimal_path, trimal_args, protein_file, codon_file,
                          protein_output, codon_output):
    """Apply trimAl's original protein column indices to whole codons."""
    # This function owns the input, output, format and column-map contract.
    # Letting a forwarded argument replace any of them could make trimAl
    # operate on another alignment or return a back-translated product, while
    # the following code still interpreted its columns as protein positions.
    reserved_options = {"-in", "-out", "-fasta", "-backtrans", "-colnumbering"}
    conflicting = sorted(reserved_options.intersection(trimal_args))
    if conflicting:
        raise ValueError(
            "--trimal-args must not contain " + ", ".join(conflicting) +
            "; these options are managed automatically for codon trimming"
        )
    proteins = {r.id: str(r.seq) for r in read_fasta_alignment(protein_file)}
    codons = {r.id: str(r.seq) for r in read_fasta_alignment(codon_file, codon=True)}
    if proteins.keys() != codons.keys() or any(
        len(codons[name]) != 3 * len(seq) for name, seq in proteins.items()
    ):
        raise ValueError("Protein/codon alignment dimensions or IDs do not match")
    options = list(trimal_args) if trimal_args else ["-automated1"]
    if "-colnumbering" not in options:
        options.append("-colnumbering")
    # Validate both products before replacing either destination.
    with tempfile.TemporaryDirectory(prefix=".msap-trimal-", dir=".") as temporary:
        trimmed_path = os.path.join(temporary, "protein.fasta")
        result = subprocess.run(
            [trimal_path, "-in", protein_file, "-out", trimmed_path, "-fasta", *options],
            check=True, stdout=subprocess.PIPE, text=True,
        )
        maps = [line.split("\t", 1)[1] for line in result.stdout.splitlines()
                if line.startswith("#ColumnsMap\t")]
        if len(maps) != 1:
            raise ValueError("trimAl did not return one #ColumnsMap; cannot preserve codon positions")
        try:
            columns = [int(item.strip()) for item in maps[0].split(",") if item.strip()]
        except ValueError as error:
            raise ValueError("Invalid trimAl column mapping") from error
        original_length = len(next(iter(proteins.values())))
        if (not columns or columns != sorted(set(columns)) or
                columns[0] < 0 or columns[-1] >= original_length):
            raise ValueError("trimAl column mapping is empty, unordered, or out of range")
        trimmed = read_fasta_alignment(trimmed_path)
        trimmed_ids = {record.id for record in trimmed}
        if trimmed_ids != proteins.keys() or len(trimmed) != len(proteins):
            missing = sorted(proteins.keys() - trimmed_ids)
            unexpected = sorted(trimmed_ids - proteins.keys())
            details = []
            if missing:
                details.append("removed IDs: " + ", ".join(missing))
            if unexpected:
                details.append("unexpected IDs: " + ", ".join(unexpected))
            raise ValueError(
                "trimAl changed the sequence set during codon trimming (" +
                "; ".join(details) + ")"
            )
        for record in trimmed:
            if record.id not in proteins or str(record.seq).upper() != "".join(
                proteins[record.id][i] for i in columns
            ).upper():
                raise ValueError(f"trimAl column mapping does not match output for {record.id}")
        with atomic_output(codon_output, inputs=[protein_file, codon_file]) as out:
            for record in trimmed:
                seq = "".join(codons[record.id][3*i:3*i+3] for i in columns)
                print(f">{record.id}\n{seq}", file=out)
        with atomic_output(protein_output, inputs=[protein_file, codon_file]) as out:
            SeqIO.write(trimmed, out, "fasta")

def run_pseudogene(args, infile, prefix, tools):
    # MACSE accepts a narrower alphabet than the rest of the workflow.  Make
    # uncertainty explicit before invocation, then retain those audit rows.
    cleaned = prefix + '.pseudogene.cleaned.fasta'
    input_report = prefix + '.input.sequence_changes.tsv'
    preprocess_nucl(infile, cleaned, input_report)
    nt, aa = macse_align(args, cleaned, prefix, tools['java'])
    merge_reports([input_report, prefix + '.sequence_changes.tsv'], prefix + '.merged.sequence_changes.tsv')
    os.replace(prefix + '.merged.sequence_changes.tsv', prefix + '.sequence_changes.tsv')
    os.remove(cleaned)
    os.remove(input_report)
    if args.notrim:
        return
    out_nt, out_aa = prefix + '.macse.codon.trimal.aln', prefix + '.macse.prot.trimal.aln'
    if args.trim_software == 'trimal':
        trim_codon_with_trimal(tools['trimal'], args.trimal_args, aa, nt, out_aa, out_nt)
        return
    nucleotides = read_fasta_alignment(nt, codon=True)
    proteins = {r.id: str(r.seq).upper() for r in read_fasta_alignment(aa)}
    length = len(nucleotides[0].seq) // 3
    count = len(nucleotides)
    keep = []
    for i in range(length):
        codons = [str(r.seq[3*i:3*i+3]).upper() for r in nucleotides]
        residues = [proteins[r.id][i] for r in nucleotides]
        if (sum('-' in c for c in codons) / count <= args.G and
                sum('N' in c for c in codons) / count <= args.N and
                residues.count('X') / count <= args.X):
            keep.append(i)
    if not keep:
        raise ValueError('No codons remain after paired pseudogene trimming')
    with atomic_output(out_nt) as handle_nt, atomic_output(out_aa) as handle_aa:
        for r in nucleotides:
            handle_nt.write(f'>{r.id}\n' + ''.join(str(r.seq[3*i:3*i+3]) for i in keep) + '\n')
            handle_aa.write(f'>{r.id}\n' + ''.join(proteins[r.id][i] for i in keep) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Align nucleotide, protein, or CDS FASTA sequences.

Workflows:
  codon  Remove terminal stops, mask internal stops as NNN, translate and align,
         back-translate to codons, then optionally trim both alignments.
  pseudogene  MACSE alignment allowing frameshifts/stops; audited export and paired trimming.
  prot   Align protein sequences, then optionally trim the alignment.
  nucl   Align nucleotide sequences, then optionally trim the alignment.""",
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  MSAP.py -i CDS.fasta -st codon -g 1
  MSAP.py -i gene.fa -st pseudogene --macse-jar macse_v2.07.jar
  MSAP.py -i protein.fasta -st prot -s muscle -t 4
  MSAP.py -i 16S.fasta -st nucl --trim-software trimal
  MSAP.py -i 16S.fasta -st nucl --trim-software trimal \\
      --trimal-args -automated1

Notes:
  * nucl/codon accept RNA: U is normalized to T and IUPAC ambiguity codes are
    normalized to N, with one audit row per edit. Inputs remain unchanged.
  * prot removes trailing stop symbols and masks internal stops, ambiguous and
    non-standard residues as X, with an audit report. Use --protein-stop-symbol
    '.' for dot-encoded stops (default: '*').
  * Input FASTA must contain at least two non-empty sequences with unique IDs.
  * Codon input must contain ungapped nucleotide sequences. A terminal excess
    of one or two bases is trimmed with a warning and audit record.
  * clustalw2 and prank do not accept duplicated sequence IDs.
  * -g/--genetic-code uses NCBI genetic-code table IDs (default: 1).
    https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi""",
    )
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('workflow options')
    trim_aln_seq_options = parser.add_argument_group('trimAlnSeq.py options')
    trimal_options = parser.add_argument_group('trimAl options')
    required.add_argument('-i', '--input', metavar='FASTA', required=True,
                          help='Input FASTA file.')
    optional.add_argument('-t', '--thread', metavar='N', default=os.cpu_count(),
                          type=int, help=f'MAFFT/MUSCLE threads; not used by PRANK, ClustalW2 or MACSE (default: {os.cpu_count()}).')
    optional.add_argument('-s', '--align_software', '--align-software', default=None,
                          choices=["mafft", "muscle", "prank", "clustalw2", "macse"],
                          help='Alignment program (default: mafft; pseudogene requires macse).')
    optional.add_argument('-n', '--notrim', action='store_true',
                          help='Skip alignment trimming, regardless of the selected trimming software.')
    optional.add_argument('-ts', '--trim-software', default='trimAlnSeq',
                          choices=['trimAlnSeq', 'trimal'],
                          help='Trimming program (default: trimAlnSeq; trimal uses -automated1 by default).')
    trim_aln_seq_options.add_argument('-G', '--G', '--trimAlnSeq-G', dest='G', default=0.2,
                                      type=float, metavar='float',
                                      help='Maximum gap ratio per site. range 0-1. default=0.2')
    trim_aln_seq_options.add_argument('-N', '--N', '--trimAlnSeq-N', dest='N', default=0.2,
                                      type=float, metavar='float',
                                      help='Maximum N ratio per site in nucleotide/codon data. default=0.2')
    trim_aln_seq_options.add_argument('-X', '--X', '--trimAlnSeq-X', dest='X', default=0.2,
                                      type=float, metavar='float',
                                      help='Maximum X ratio per site in protein data. default=0.2')
    trimal_options.add_argument('--trimal-args', nargs=argparse.REMAINDER, default=[],
                                help='Arguments passed to trimAl; must be final and override -automated1. In codon/pseudogene mode, do not pass -in, -out, -fasta, -backtrans or -colnumbering.')
    optional.add_argument('-st', '--seqtype', '--seq-type', default='codon',
                          choices=['codon', 'prot', 'nucl', 'pseudogene'],
                          help='Input sequence type (default: codon).')
    optional.add_argument('-g', '--genetic_code', '--genetic-code', metavar='TABLE', default=1,
                          type=int, help='NCBI genetic-code table for codon workflow (default: 1).')
    optional.add_argument('-h', '--help', action='help',
                          help="Show program's help message and exit.")
    optional.add_argument('-v', '--version', action='version', version='v1.0.0',
                          help="Show program's version number and exit.")
    
    add_macse_options(optional)
    add_protein_options(optional)
    args = parser.parse_args()
    try:
        normalize_workflow(args)
    except ValueError as error:
        parser.error(str(error))
    if args.thread is None or args.thread < 1:
        parser.error("--thread must be a positive integer")
    if args.seqtype in {"codon", "pseudogene"} and args.genetic_code not in CodonTable.unambiguous_dna_by_id:
        parser.error("Unknown NCBI genetic-code table")
    infile=args.input
    thread = args.thread
    genetic_code = args.genetic_code
    align_software = args.align_software
    seqtype = args.seqtype
    notrim= args.notrim
    trim_software = args.trim_software
    trimal_args = args.trimal_args

    for option, value in (("-G", args.G), ("-N", args.N), ("-X", args.X)):
        if not 0 <= value <= 1:
            parser.error(f"{option} must be between 0 and 1")
    if not notrim and trimal_args and trim_software != "trimal":
        parser.error("--trimal-args requires --trim-software trimal")
    try:
        validate_input_fasta(infile, seqtype, args.protein_stop_symbol)
    except ValueError as error:
        if seqtype in {"codon", "pseudogene"}:
            rejected_input(infile, get_prefix(infile), seqtype, genetic_code, str(error))
        parser.error(str(error))


    softwares = ["java" if seqtype == "pseudogene" else align_software]
    if not notrim and trim_software == "trimal":
        softwares.append("trimal")
    software_path = check_dependencies(softwares=softwares)

    # run main
    prefix = get_prefix(infile)
    infile = os.path.realpath(infile)
    if seqtype == "pseudogene":
        run_pseudogene(args, infile, prefix, software_path)
        sys.exit(0)
    if seqtype == "codon":
        infile_tmp = prefix + '.CDS.fasta'
        try:
            remove_stop_codon(infile, outfile=infile_tmp, genetic_code=genetic_code)
        except EmptyCodonLocusError:
            # Exit code 3 is a deliberate locus skip, recognized by MSAP_batch.
            sys.exit(3)
        infile = infile_tmp
        translate_seq(infile, outfile=prefix + '.pep.fasta', genetic_code=genetic_code)
        align_infile = os.path.realpath(f"{prefix}.pep.fasta")
    elif seqtype == 'nucl':
        align_infile = os.path.realpath(prefix + '.DNA.fasta')
        preprocess_nucl(infile, align_infile, prefix + '.sequence_changes.tsv')
    else:
        align_infile = os.path.realpath(prefix + '.protein.cleaned.fasta')
        preprocess_protein(infile, align_infile, prefix + '.sequence_changes.tsv', args.protein_stop_symbol)
    
    ## run align
    if align_software == "mafft":
        if seqtype == "codon" or seqtype == "prot":
            align_outfile = os.path.realpath(f"{prefix}.mafft.prot.aln")
        elif seqtype == "nucl":
            align_outfile = os.path.realpath(f"{prefix}.mafft.nucl.aln")
        with atomic_output(align_outfile, inputs=[align_infile]) as out:
            subprocess.run(
                [software_path["mafft"], "--amino" if seqtype in {"prot", "codon"} else "--nuc",
                 "--thread", str(thread), "--quiet", "--auto", align_infile],
                stdout=out,
                check=True,
            )
    
    
    elif align_software == "muscle":
        if seqtype == "codon" or seqtype == "prot":
            align_outfile = os.path.realpath(f"{prefix}.muscle.prot.aln")
        elif seqtype == "nucl":
            align_outfile = os.path.realpath(f"{prefix}.muscle.nucl.aln")
        subprocess.run(
            [software_path["muscle"], "-threads", str(thread), "-align", align_infile, "-output", align_outfile],
            check=True,
        )
    
    
    elif align_software == "prank":
        subprocess.run(
            [software_path["prank"], f"-d={align_infile}", f"-o={prefix}.prank.aln", "-f=fasta",
             "-DNA" if seqtype == "nucl" else "-protein"],
            check=True,
        )
        if seqtype == "codon" or seqtype == "prot":
            os.rename(f"{prefix}.prank.aln.best.fas", f"{prefix}.prank.prot.aln")
            align_outfile = os.path.realpath(f"{prefix}.prank.prot.aln")
        elif seqtype == "nucl":
            os.rename(f"{prefix}.prank.aln.best.fas", f"{prefix}.prank.nucl.aln")
            align_outfile = os.path.realpath(f"{prefix}.prank.nucl.aln")

        
    elif align_software == "clustalw2":
        suffix = "nucl" if seqtype == "nucl" else "prot"
        clustal_type = "DNA" if seqtype == "nucl" else "PROTEIN"
        align_outfile = os.path.realpath(f"{prefix}.clustalw2.{suffix}.aln")
        # ClustalW creates guide trees beside its input. Keep all such files private.
        with tempfile.TemporaryDirectory(prefix=".msap-clustalw-", dir=".") as temporary:
            local_input = os.path.abspath(os.path.join(temporary, "input.fasta"))
            shutil.copyfile(align_infile, local_input)
            subprocess.run(
                [software_path["clustalw2"], f"-INFILE={local_input}",
                 f"-TYPE={clustal_type}", "-OUTPUT=FASTA", f"-OUTFILE={align_outfile}"],
                check=True,
            )

    else:
        print("Please input a align software. [matff, muscle, prank, clustalw2]")
        sys.exit()
    
    aligned = read_fasta_alignment(align_outfile)
    source_sequences = {r.id: str(r.seq).upper().replace('-', '') for r in SeqIO.parse(align_infile, 'fasta')}
    aligned_sequences = {r.id: str(r.seq).upper().replace('-', '') for r in aligned}
    if source_sequences != aligned_sequences:
        raise ValueError('Alignment changed taxon IDs or non-gap sequence content; results rejected')

    ## run AA2codon.py
    if seqtype == "codon":
        subprocess.run(
            [
                sys.executable,
                os.path.join(os.path.dirname(os.path.abspath(__file__)), "AA2Codon.py"),
                "-c", infile,
                "-p", f"{prefix}.{align_software}.prot.aln",
                "-o", f"{prefix}.{align_software}.codon.aln",
                "-g", str(genetic_code),
            ],
            check=True,
        )
        if not notrim and trim_software == "trimal":
            trim_codon_with_trimal(
                software_path["trimal"], trimal_args,
                f"{prefix}.{align_software}.prot.aln",
                f"{prefix}.{align_software}.codon.aln",
                f"{prefix}.{align_software}.prot.trimal.aln",
                f"{prefix}.{align_software}.codon.trimal.aln",
            )
        elif not notrim:
            trim_alignment(
                trim_software, software_path.get("trimal"), trimal_args,
                f"{prefix}.{align_software}.codon.aln",
                f"{prefix}.{align_software}.codon.trimal.aln",
                "codon", args.G, args.N, args.X,
            )
            trim_alignment(
                trim_software, software_path.get("trimal"), trimal_args,
                f"{prefix}.{align_software}.prot.aln",
                f"{prefix}.{align_software}.prot.trimal.aln",
                "prot", args.G, args.N, args.X,
            )
    else:
        ## run trimAlnSeq.py
        if notrim == False:
            trim_alignment(
                trim_software, software_path.get("trimal"), trimal_args,
                f"{prefix}.{align_software}.{seqtype}.aln",
                f"{prefix}.{align_software}.{seqtype}.trimal.aln",
                seqtype, args.G, args.N, args.X,
            )

    if seqtype == "codon":
        os.remove(infile)
    elif seqtype in {'nucl', 'prot'}:
        os.remove(align_infile)
