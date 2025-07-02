---
theme: default
title: Orion-Kmer User Manual
---

# Orion-Kmer User Manual

Welcome to the Orion-Kmer user manual! Orion-Kmer is a blazing-fast k-mer toolkit written in Rust, designed for efficient k-mer counting, database construction, comparison, and querying.

This manual provides a comprehensive guide to installing and using Orion-Kmer.

## Table of Contents

*   [Installation](#installation)
*   [Global Options](#global-options)
*   Tutorials
    *   [Counting K-mers (`count`)](./tutorials/count_tutorial.md)
    *   [Building K-mer Databases (`build`)](./tutorials/build_tutorial.md)
    *   [Comparing K-mer Databases (`compare`)](./tutorials/compare_tutorial.md)
    *   [Querying Reads against a Database (`query`)](./tutorials/query_tutorial.md)
    *   [Classifying Sequences (`classify`)](./tutorials/classify_tutorial.md)
*   [Example Workflow](#example-workflow)

## Installation

Currently, Orion-Kmer must be built from source using Cargo:

```bash
# Clone the repository (assuming it's in a public repo later)
# git clone <repo_url>
# cd orion-kmer
cargo build --release
# The binary will be at target/release/orion-kmer
```
*Ensure Rust and Cargo are installed on your system. You can find installation instructions at [rust-lang.org](https://www.rust-lang.org/tools/install).*

## Global Options

Orion-Kmer supports the following global options that can be used with any subcommand:

*   `-t, --threads <THREADS>`: Number of threads to use. Defaults to the number of logical cores if set to 0.
    *   Example: `orion-kmer -t 4 count ...`
*   `-v, --verbose`: Increase verbosity. Use `-v` for informational messages and `-vv` for debug messages.
    *   Example: `orion-kmer -vv build ...`
*   `-h, --help`: Print help information for the main command or a specific subcommand.
    *   Example: `orion-kmer --help` or `orion-kmer count --help`
*   `-V, --version`: Print version information.

## Tutorials

The following tutorials demonstrate how to use each of Orion-Kmer's subcommands with example data.

*   **[Counting K-mers (`count`)](./tutorials/count_tutorial.md)**: Learn how to count k-mers in FASTA/FASTQ files.
*   **[Building K-mer Databases (`build`)](./tutorials/build_tutorial.md)**: Learn how to create compact k-mer databases from genome assemblies.
*   **[Comparing K-mer Databases (`compare`)](./tutorials/compare_tutorial.md)**: Learn how to compare two k-mer databases and interpret the results.
*   **[Querying Reads against a Database (`query`)](./tutorials/query_tutorial.md)**: Learn how to screen short reads against a k-mer database.

## Example Workflow

Here's a typical workflow using Orion-Kmer:

1.  **Build a database for a reference genome:**
    ```bash
    orion-kmer build -k 31 -t 8 -g reference.fasta -o ref.db
    ```

2.  **Build a database for a newly assembled genome:**
    ```bash
    orion-kmer build -k 31 -t 8 -g new_assembly.fasta -o new.db
    ```

3.  **Compare the two assemblies:**
    ```bash
    orion-kmer compare --db1 ref.db --db2 new.db -o report.json
    ```

4.  **Check which reads from a sequencing experiment belong to the reference:**
    ```bash
    orion-kmer query -d ref.db -r experiment.fastq.gz -o matched_reads.txt -c 5 -t 8
    ```

---
*For more details on command-line arguments, run `orion-kmer <subcommand> --help`.*
