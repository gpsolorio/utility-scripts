# Make mini-bam data using the synthetic mini-sam

# 1. Create directories
mkdir -p ../mini_bam
mkdir -p ../../cramino/outputs/tiny_test


samtools view -b \
    -o ../mini_bam/tiny_unsorted.bam \
    ../mini_sam/tiny_alignments.sam

# 3. Sort bam file
samtools sort \
    -o ../mini_bam/tiny_sorted.bam \
    ../mini_bam/tiny_unsorted.bam

# 3. Index bam file
samtools index ../mini_bam/tiny_sorted.bam




### Use samtools calmd to calculate the tags

# 1. Index the synthetic reference
samtools faidx ../mini_sam/tiny_reference.fasta

# 2. Add NM and MD alignment tags
samtools calmd -b \
  ../mini_bam/tiny_sorted.bam \
  ../mini_sam/tiny_reference.fasta \
  > ../mini_bam/tiny_calmd.bam

# 3. Index the updated BAM
samtools index ../mini_bam/tiny_calmd.bam