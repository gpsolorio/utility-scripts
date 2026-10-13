Toy SAM error-rate exercise
===========================

Files
- toy_reference.fa: 1 kb synthetic reference called chrToy
- toy_reads.sam: small SAM with known substitutions, insertions, deletions, and mixed-error reads
- sam_error_rates.py: Python 3 standard-library analysis tool

Run
    python3 sam_error_rates.py toy_reads.sam --reference toy_reference.fa

The SAM is synthetic and intended for learning/testing, not biological inference.
Rates depend on the denominator:
- substitution rate: mismatches / aligned M/= /X columns
- indel-event rate: insertion/deletion events / reference-aligned bases (M/= /X plus D)
- indel-base rate: inserted/deleted bases / reference-aligned bases

The tool uses CIGAR strings to locate insertions/deletions and compares aligned bases
to the supplied reference for substitutions. It does not require MD tags. It skips
unmapped reads and handles M, =, X, I, D, N, S, H, and P CIGAR operations.
