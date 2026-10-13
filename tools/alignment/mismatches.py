import argparse
import collections
import re
import subprocess
import sys

# Counter to track (Reference_Base -> Read_Base)
mismatch_counts = collections.Counter()

# Translation table for complementing DNA bases
COMPLEMENT_TABLE = str.maketrans("ACGTN", "TGCAN")

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence string."""
    return seq.translate(COMPLEMENT_TABLE)[::-1]

def parse_md_and_match_read(cigar, md_tag, read_seq):
    """Parses the MD tag alongside the read sequence, properly accounting for soft-clipping."""
    if not cigar or cigar == "*" or not md_tag or not read_seq or read_seq == "*":
        return

    # 1. Parse CIGAR into standard tuples: (length, op)
    cigar_ops = [(int(ln), op) for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]

    # 2. Tokenize MD tag into numbers (matches), deletions (^ACT), or single mismatches (A,C,G,T)
    md_tokens = re.findall(r"(\d+|\^[A-Za-z]+|[A-Za-z])", md_tag)

    read_idx = 0
    aligned_read_bases = []

    for op_len, op_type in cigar_ops:
        if op_type in ("M", "=", "X"):
            # These bases align to the reference and are evaluated by the MD tag
            for _ in range(op_len):
                if read_idx < len(read_seq):
                    aligned_read_bases.append(read_seq[read_idx].upper())
                    read_idx += 1
        elif op_type == "I":
            # Insertions exist in read, but not reference (MD tag skips them)
            read_idx += op_len
        elif op_type == "S":
            # Soft-clipped bases exist in read string but are ignored by MD tag.
            read_idx += op_len
        elif op_type in ("D", "N"):
            # Deletions exist in reference but not read. Insert padding.
            for _ in range(op_len):
                aligned_read_bases.append(None)

    # 3. Traverse the MD tokens over our adjusted, clean array
    aligned_idx = 0
    for token in md_tokens:
        if token.isdigit():
            aligned_idx += int(token)
        elif token.startswith("^"):
            del_len = len(token) - 1
            aligned_idx += del_len
        else:
            ref_base = token.upper()

            if aligned_idx < len(aligned_read_bases):
                read_base = aligned_read_bases[aligned_idx]
                # Only log the mismatch if the bases are genuinely different
                if read_base and ref_base in "ACGT" and read_base in "ACGT" and ref_base != read_base:
                    mismatch_counts[(ref_base, read_base)] += 1

            aligned_idx += 1


def process_stream(stream):
    """Processes an open stream of SAM lines."""
    for line in stream:
        if line.startswith("@"):
            continue

        fields = line.strip().split("\t")
        if len(fields) < 11:
            continue

        try:
            flag = int(fields[1])
        except (ValueError, IndexError):
            continue

        # Skip unmapped reads (flag bit 4)
        if flag & 4:
            continue

        cigar = fields[5]
        read_seq = fields[9]

        # FIX: If flag bit 16 is set, the read is reverse stranded.
        # samtools view outputs reverse-complemented sequences for these.
        # We must flip it back to match the forward reference orientation used by the MD tag.
        if flag & 16 and read_seq != "*":
            read_seq = reverse_complement(read_seq)

        # Locate the MD tag out of the remaining optional tag arrays
        md_tag = None
        for field in fields[11:]:
            if field.startswith("MD:Z:"):
                md_tag = field[5:]
                break

        if md_tag:
            try:
                parse_md_and_match_read(cigar, md_tag, read_seq)
            except Exception:
                continue


def main():
    parser = argparse.ArgumentParser(
        description="Count mismatched nucleotide frequencies across multiple BAM files linearly."
    )
    parser.add_argument(
        "bam_files",
        metavar="BAM",
        type=str,
        nargs="*",
        help="One or more BAM files to process. If omitted, reads from standard input.",
    )

    args = parser.parse_args()

    if args.bam_files:
        for bam_path in args.bam_files:
            print(f"Processing: {bam_path}...", file=sys.stderr)
            try:
                process = subprocess.Popen(
                    ["samtools", "view", bam_path],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                process_stream(process.stdout)
                process.wait()
                
                if process.returncode != 0:
                    err_msg = process.stderr.read().strip()
                    print(f"Error processing {bam_path}: {err_msg}", file=sys.stderr)

            except FileNotFoundError:
                print("Error: 'samtools' not found in PATH.", file=sys.stderr)
                sys.exit(1)
            except Exception as e:
                print(f"Error reading {bam_path}: {e}", file=sys.stderr)
    else:
        print("No files specified. Reading from standard input...", file=sys.stderr)
        process_stream(sys.stdin)

    # --- Print Clean TSV Output ---
    print("reference_base\tread_base\tcount")
    for (ref, read), count in sorted(mismatch_counts.items()):
        print(f"{ref}\t{read}\t{count}")


if __name__ == "__main__":
    main()
