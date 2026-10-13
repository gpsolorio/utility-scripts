#!/usr/bin/env bash

# Create output directory
mkdir -p tools/cramino/outputs/tiny_test

# Generate summary statistics
cramino --format tsv \
  mini_data/cramino_mini_data/alignments/tiny_calmd.bam \
  > tools/cramino/outputs/tiny_test/cramino_stats.tsv

# Generate read-length histogram
cramino \
  --hist-count=tools/cramino/outputs/tiny_test/read_length_histogram.tsv \
  mini_data/cramino_mini_data/alignments/tiny_calmd.bam \
  > /dev/null