# K-mer Enrichment Analyzer

A command-line tool written in Rust to compute and analyze k-mer enrichment by comparing a reference FASTA file and a BAM alignment file.

The tool identifies k-mers that are over-represented (enriched) in the sequencing reads (BAM) compared to their frequency in the reference genome (FASTA). This can be useful for various bioinformatics analyses, such as identifying repetitive elements, sequencing artifacts, or regions with high mutation rates.

## Features

- **K-mer Counting**: Efficiently counts k-mers in both FASTA and BAM files.
- **3-bit Encoding**: Uses a memory-efficient 3-bit encoding for the 6-character alphabet {A, C, G, T, N, -}, allowing for k-mer sizes up to 21.
- **Enrichment Calculation**: For each k-mer, it calculates its proportion in both the reference and the BAM file, and computes an enrichment score (`bam_proportion / ref_prop`).
- **Filtering**: Skips unmapped, duplicate, supplementary, and low-quality reads in the BAM file.
- **JSON Output**: Outputs the results in two separate, easy-to-parse JSON files: one for canonical k-mers (`ACGT`) and one for k-mers containing `N` or `-`.
- **Sorting**: The output lists are sorted by enrichment score in descending order, bringing the most over-represented k-mers to the top.

## How It Works

1.  **Reference K-mer Counting**: The tool first reads the provided reference FASTA file (which must be indexed with `samtools faidx`). It iterates through each contig, calculating and counting all k-mers present in the reference sequence.
2.  **BAM K-mer Counting**: Next, it processes the input BAM file. For each sequencing read that passes quality filters, it calculates and counts all k-mers.
3.  **Proportion & Enrichment**: With the two sets of k-mer counts, it calculates the proportion of each unique k-mer relative to the total number of k-mers in the reference and the BAM file, respectively. The enrichment score is then calculated as the ratio of these proportions.
4.  **Output Generation**: The final list of k-mers and their associated data is partitioned into two groups:
    - `canonical.json`: Contains k-mers composed only of A, C, G, and T.
    - `with_n.json`: Contains k-mers that include 'N' (unknown base) or '-' (gap/other).

    Both files are written to the current directory.

## Prerequisites

Before building, ensure you have the **Rust toolchain** and `libhts-dev` installed.

- **Rust**: [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install)
- **HTSlib**: This is a dependency for `rust-htslib`. On Debian/Ubuntu, you can install it with:
  ```sh
  sudo apt-get update
  sudo apt-get install libhts-dev
  ```

You will also need an indexed FASTA file. If you have `my_reference.fasta`, you can index it using `samtools`:
```sh
samtools faidx my_reference.fasta
```

## Building

Clone the repository and build the project using Cargo:

```sh
git clone <repository_url>
cd <repository_name>
cargo build --release
```

The compiled binary will be located at `target/release/kmer_enrichment`.

## Usage

Run the tool from your terminal, providing paths to your reference FASTA, BAM file, a k-mer size, and an output path (note: the output path is a required argument but is currently unused, as output files are hardcoded).

### Command-Line Arguments

-   `-r, --ref <FASTA>`: Path to the reference FASTA file (must be indexed, with a `.fai` file present).
-   `-b, --bam <BAM>`: Path to the input BAM file.
-   `-k, --kmer <K>`: The k-mer size to use (an integer, e.g., 7, 15, 21). Must satisfy `3*k <= 64`.
-   `-o, --output <OUTPUT>`: A required placeholder for the output path. The tool currently writes to `canonical.json` and `with_n.json` in the working directory.

### Example

```sh
./target/release/kmer_enrichment \
    --ref /path/to/your/reference.fasta \
    --bam /path/to/your/alignments.bam \
    --kmer 21 \
    --output ./results/
```

This will produce two files in the directory where you run the command: `canonical.json` and `with_n.json`.

## Output Format

The output JSON files contain an array of objects, where each object represents a k-mer and its statistics.

### `canonical.json` / `with_n.json`

```json
[
  {
    "kmer": "ATCGATCGATCGATCGATCGA",
    "ref_prop": 0.000015,
    "bam_prop": 0.0025,
    "enrichment": 166.66
  },
  {
    "kmer": "GATTACAGATTACAGATTACA",
    "ref_prop": 0.000021,
    "bam_prop": 0.0018,
    "enrichment": 85.71
  },
  // ... more k-mers
]
```

-   `kmer`: The DNA sequence of the k-mer.
-   `ref_prop`: The proportion of this k-mer relative to the total k-mers in the reference FASTA.
-   `bam_prop`: The proportion of this k-mer relative to the total k-mers in the BAM file.
-   `enrichment`: The calculated enrichment score (`bam_prop / ref_prop`). `null` if the `ref_prop` is zero, as this would result in division by zero.