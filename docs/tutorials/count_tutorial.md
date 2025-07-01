---
layout: default
title: Count K-mers Tutorial
---

# Tutorial: Counting K-mers with `orion-kmer count`

[Back to Main Manual](../index.md)

The `orion-kmer count` subcommand is used to count all k-mers in one or more input FASTA or FASTQ files. It outputs the k-mers and their counts in a simple tab-separated format.

## Command Overview

```bash
orion-kmer count -k <KMER_SIZE> -i <INPUT_FILE1> [INPUT_FILE2 ...] -o <OUTPUT_TSV> [-m <MIN_COUNT>]
```

### Arguments:

*   `-k, --kmer-size <INT>`: **Required.** The length of the k-mer (1-32).
*   `-i, --input <FILE>...`: **Required.** One or more input FASTA/FASTQ files. These can be gzipped (e.g., `sequences.fasta.gz`).
*   `-o, --output <FILE>`: **Required.** The path to the output file where k-mer counts will be written (e.g., `counts.tsv`).
*   `-m, --min-count <INT>`: Optional. The minimum count required for a k-mer to be reported. Defaults to `1`.

## Tutorial with Toy Data

Let's use the toy data we created to count k-mers. First, ensure you have the `orion-kmer` binary accessible in your PATH or provide the full path to it. You'll also need the toy FASTA file `genome1.fasta`.

**Toy Data:** `genome1.fasta`
```fasta
>genome1_seq1
ATGCGTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCG
>genome1_seq2
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
```

### Step 1: Running `orion-kmer count`

Let's count 5-mers from `genome1.fasta`.

Navigate to the `docs/tutorials/data` directory in your terminal or adjust paths accordingly. For this example, we'll assume you are running the command from a directory where `orion-kmer` is accessible and `genome1.fasta` is in a subdirectory `data`.

```bash
# Assuming orion-kmer is in your PATH
# and you are in the 'docs/tutorials' directory
mkdir -p output # Create an output directory for results
orion-kmer count -k 5 -i data/genome1.fasta -o output/genome1_counts.tsv
```

If `orion-kmer` is not in your PATH, you would replace `orion-kmer` with the actual path to the compiled binary (e.g., `../../orion-kmer/target/release/orion-kmer`).

### Step 2: Examining the Output

The command will create a file named `genome1_counts.tsv` in the `output` directory (relative to `docs/tutorials`). The content will look something like this (order might vary, and actual k-mers/counts depend on the exact sequences and k-mer size):

**Expected `output/genome1_counts.tsv` (example for k=5):**
```tsv
ATGCG	1
TGCGT	1
GCGTA	1
CGTAG	1
GTAGC	1
TAGCA	1
AGCAT	1
GCATC	1
CATCG	2
ATCGA	8
TCGAT	8
CGATC	8
GATCG	7
// ... and so on for all 5-mers
```
*Note: K-mers are canonical by default in orion-kmer. This means the lexicographically smaller of a k-mer and its reverse complement is stored and counted. The exact output will reflect this.*

Let's manually verify a few for `k=3` from `ATGCGTAGC`:
- `ATG` (rev-comp: `CAT`) -> `ATG`
- `TGC` (rev-comp: `GCA`) -> `GCA`
- `GCG` (rev-comp: `CGC`) -> `CGC`
- `CGT` (rev-comp: `ACG`) -> `ACG`
- `GTA` (rev-comp: `TAC`) -> `TAC`
- `TAG` (rev-comp: `CTA`) -> `CTA`
- `AGC` (rev-comp: `GCT`) -> `AGC`

If you run `orion-kmer count -k 3 -i data/genome1.fasta -o output/genome1_counts_k3.tsv`, you can check these.

### Step 3: Using `--min-count`

If you only want to see k-mers that appear at least, say, 5 times:
```bash
orion-kmer count -k 5 -i data/genome1.fasta -o output/genome1_counts_min5.tsv -m 5
```
The output file `genome1_counts_min5.tsv` would then only contain k-mers that appeared 5 or more times.

## Next Steps

*   Try different k-mer sizes.
*   Use multiple input files: `orion-kmer count -k 7 -i data/genome1.fasta data/genome2.fasta -o output/combined_counts.tsv`
*   Experiment with gzipped input files if you have them.

---
[Back to Main Manual](../index.md) | [Next Tutorial: Building Databases (`build`)](./build_tutorial.md)
