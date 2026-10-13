# Tiny ground-truth error-analysis dataset

Designed to test substitution-spectrum and indel-analysis code.

- `reference/truth_reference.fasta`: 64 bp reference
- `reads/truth_reads.fastq`: four short synthetic reads
- `alignments/truth_alignments.sam`: hand-written alignment fixture with `=`, `X`, `I`, `D` CIGAR operations
- `expected/GROUND_TRUTH.txt`: intended errors and expected totals

Optional SAMtools conversion, from this directory:

```bash
samtools faidx reference/truth_reference.fasta
samtools view -bS alignments/truth_alignments.sam | samtools sort -o alignments/truth_sorted.bam
samtools index alignments/truth_sorted.bam
samtools calmd -b alignments/truth_sorted.bam reference/truth_reference.fasta > alignments/truth_calmd.bam
samtools index alignments/truth_calmd.bam
```

This tiny fixture tests error-counting logic; it is not realistic sequencing data and does not model ONT quality scores or context-dependent error rates.
