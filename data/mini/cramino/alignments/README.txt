Tiny synthetic long-read test data
==================================
5 random DNA reads (80-150 bp) with synthetic quality strings.
Includes FASTQ, FASTA reference, and SAM alignments to matching artificial contigs.
These are tiny test fixtures, not biologically realistic long reads.
To create BAM files with samtools installed:
  samtools view -b -o tiny_unsorted.bam tiny_alignments.sam
  samtools sort -o tiny.bam tiny_unsorted.bam
  samtools index tiny.bam
