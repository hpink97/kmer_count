#!/usr/bin/env bash
set -euo pipefail

# Download BAM + index
mkdir -p data
curl -L https://raw.githubusercontent.com/IARCbioinfo/data_test/master/BAM/NA06984_N.bam -o data/NA06984_N.bam
curl -L https://raw.githubusercontent.com/IARCbioinfo/data_test/master/BAM/NA06984_N.bam.bai -o data/NA06984_N.bam.bai

# Download reference FASTA (choose one)
curl -L https://raw.githubusercontent.com/IARCbioinfo/data_test/master/REF/TP53.fasta -o data/TP53.fasta

# Optional: index FASTA for tools that need it
samtools faidx data/TP53.fasta


# build if no release binary exists
if [ ! -f ./target/release/]; then
    cargo build --release
fi


# Run kmer counter
./target/release/kmer_counter \
    --bam data/NA06984_N.bam \
    --ref data/TP53.fasta \
    --kmer 10 \
    --output kmer_counts.tsv


rm -rf data