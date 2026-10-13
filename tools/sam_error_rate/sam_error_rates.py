#!/usr/bin/env python3
"""
Estimate substitution and indel error patterns from a reference FASTA and SAM.

Definitions:
- Substitution rate = mismatching aligned query/reference bases / aligned M, =, X columns.
- Indel-event rate = insertion + deletion events / reference-aligned bases (M, =, X, D).
- Indel-base rate = inserted + deleted bases / reference-aligned bases.
- Insertions are summarized by inserted query nucleotide; deletions by deleted reference nucleotide.
  Deletion base identities require a reference FASTA.
This is a small teaching tool, not a production-grade variant/error caller.
"""
import argparse
from collections import Counter
from pathlib import Path

BASES = "ACGT"
COMP = {b: set(BASES) - {b} for b in BASES}
TRANSITIONS = {("A","G"), ("G","A"), ("C","T"), ("T","C")}
TRANSVERSIONS = {(a,b) for a in BASES for b in BASES if a != b and (a,b) not in TRANSITIONS}

def read_fasta(path):
    seqs, name, chunks = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs

def parse_cigar(cigar):
    import re
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not ops or "".join(n+op for n,op in ops) != cigar:
        raise ValueError(f"Invalid/unsupported CIGAR: {cigar}")
    return [(int(n), op) for n,op in ops]

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sam", help="Input SAM file")
    ap.add_argument("--reference", required=True, help="Reference FASTA used for alignment")
    args = ap.parse_args()
    refs = read_fasta(args.reference)
    subs = Counter()
    ins_bases, del_bases = Counter(), Counter()
    indel_sizes = Counter()
    read_counts = Counter()
    aligned_query_cols = ref_aligned_bases = mismatch_count = indel_events = indel_bases_total = 0
    for line in Path(args.sam).read_text().splitlines():
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq = f[:10]
        flag, pos = int(flag), int(pos)
        if flag & 4 or rname == "*" or cigar == "*":
            continue
        if rname not in refs:
            raise ValueError(f"{qname}: reference {rname!r} absent from FASTA")
        seq = seq.upper()
        ref = refs[rname]
        ri, qi = pos - 1, 0
        for n, op in parse_cigar(cigar):
            if op in ("M", "=", "X"):
                for k in range(n):
                    if ri+k >= len(ref) or qi+k >= len(seq):
                        raise ValueError(f"{qname}: CIGAR extends beyond sequence/reference")
                    rb, qb = ref[ri+k], seq[qi+k]
                    if rb in BASES and qb in BASES:
                        aligned_query_cols += 1
                        ref_aligned_bases += 1
                        if rb != qb:
                            subs[rb+qb] += 1
                            mismatch_count += 1
                ri += n
                qi += n
            elif op == "I":
                inserted = seq[qi:qi+n]
                ins_bases.update(b for b in inserted if b in BASES)
                indel_sizes[n] += 1
                indel_events += 1
                indel_bases_total += n
                qi += n
            elif op == "D":
                deleted = ref[ri:ri+n]
                del_bases.update(b for b in deleted if b in BASES)
                indel_sizes[n] += 1
                indel_events += 1
                indel_bases_total += n
                ref_aligned_bases += n
                ri += n
            elif op == "N":
                ri += n
            elif op == "S":
                qi += n
            elif op in ("H", "P"):
                pass
        read_counts["mapped_reads"] += 1

    aligned_den = aligned_query_cols
    ref_den = ref_aligned_bases
    total_subs = sum(subs.values())
    print("=== Input summary ===")
    print(f"Mapped reads analyzed: {read_counts['mapped_reads']}")
    print(f"Aligned M/= /X query-reference columns: {aligned_query_cols}")
    print(f"Reference-aligned bases including deletions: {ref_aligned_bases}")
    print()
    print("=== Substitutions ===")
    print(f"Substitutions observed: {total_subs}")
    print(f"Substitution rate: {total_subs}/{aligned_den} = {total_subs/aligned_den:.6%}" if aligned_den else "Substitution rate: n/a")
    print("Type\tCount\tFraction of substitutions")
    for a in BASES:
        for b in BASES:
            if a != b:
                c = subs[a+b]
                frac = c/total_subs if total_subs else 0
                label = "transition" if (a,b) in TRANSITIONS else "transversion"
                print(f"{a}>{b}\t{c}\t{frac:.2%}\t{label}")
    print(f"Transitions total: {sum(subs[x] for x in ('AG','GA','CT','TC'))}")
    print(f"Transversions total: {sum(subs[x] for x in ('AC','AT','CA','CG','GC','GT','TA','TG'))}")
    print()
    print("=== Insertions and deletions ===")
    print(f"Indel events: {indel_events}")
    print(f"Indel-event rate: {indel_events}/{ref_den} = {indel_events/ref_den:.6%}" if ref_den else "Indel-event rate: n/a")
    print(f"Indel bases: {indel_bases_total}")
    print(f"Indel-base rate: {indel_bases_total}/{ref_den} = {indel_bases_total/ref_den:.6%}" if ref_den else "Indel-base rate: n/a")
    print("Indel size (bp)\tEvents")
    for size in sorted(indel_sizes):
        print(f"{size}\t{indel_sizes[size]}")
    print("Inserted nucleotide\tCount")
    for b in BASES:
        print(f"{b}\t{ins_bases[b]}")
    print("Deleted reference nucleotide\tCount")
    for b in BASES:
        print(f"{b}\t{del_bases[b]}")
    print()
    print("Interpretation note: insertion-base counts describe inserted sequence composition; deletion-base counts describe")
    print("which reference bases were removed. Compare proportions, not raw counts alone, when exposure opportunities differ.")

if __name__ == "__main__":
    main()
