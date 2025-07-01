---
layout: default
title: Build K-mer Database Tutorial
---

# Tutorial: Building K-mer Databases with `orion-kmer build`

[Back to Main Manual](../index.md)

The `orion-kmer build` subcommand scans genome assemblies (FASTA format) and creates a compact, binary database file (`.db`). This database contains the set of *unique* k-mers found in the input genome(s). These database files are then used by other subcommands like `compare` and `query`.

## Command Overview

```bash
orion-kmer build -k <KMER_SIZE> -g <GENOME_FILE1> [GENOME_FILE2 ...] -o <OUTPUT_DB>
```

### Arguments:

*   `-k, --kmer-size <INT>`: **Required.** The length of the k-mer (1-32) to be stored in the database. This should match the k-mer size used for other operations (e.g., querying).
*   `-g, --genomes <FILE>...`: **Required.** One or more input genome assembly files in FASTA format. These can be gzipped.
*   `-o, --output <FILE>`: **Required.** The path where the output binary k-mer database file (e.g., `my_genome.db`) will be saved.

## Tutorial with Toy Data

We will use our toy FASTA files `genome1.fasta` and `genome2.fasta` to create k-mer databases.

**Toy Data:**
*   `genome1.fasta`:
    ```fasta
    >genome1_seq1
    ATGCGTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCG
    >genome1_seq2
    CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
    ```
*   `genome2.fasta`:
    ```fasta
    >genome2_seq1
    TACGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
    >genome2_seq2
    GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
    ```

### Step 1: Building a Database for `genome1.fasta`

Let's build a database for `genome1.fasta` using a k-mer size of 7.

Ensure `orion-kmer` is accessible and you are in the `docs/tutorials` directory or adjust paths.
```bash
# Assuming orion-kmer is in your PATH
# and you are in the 'docs/tutorials' directory
# (output directory created in previous tutorial)
orion-kmer build -k 7 -g data/genome1.fasta -o output/genome1.db
```

This command will:
1.  Read `data/genome1.fasta`.
2.  Extract all unique 7-mers (using canonical representation).
3.  Save these unique 7-mers into a binary file named `output/genome1.db`.

### Step 2: Building a Database for `genome2.fasta`

Similarly, let's build a database for `genome2.fasta` with the same k-mer size.
```bash
orion-kmer build -k 7 -g data/genome2.fasta -o output/genome2.db
```
This creates `output/genome2.db` containing unique 7-mers from `genome2.fasta`.

### Step 3: Building a Combined Database

You can also provide multiple genome files to the `build` command to create a single database containing unique k-mers from all specified files.

```bash
orion-kmer build -k 7 -g data/genome1.fasta data/genome2.fasta -o output/combined_genomes.db
```
The file `output/combined_genomes.db` will contain all unique 7-mers present in either `genome1.fasta` or `genome2.fasta` (or both).

### What's in the `.db` file?

The `.db` file is a binary file optimized for quick loading and querying by Orion-Kmer. It's not human-readable directly. It essentially stores a set of k-mer hashes. The `compare` and `query` subcommands know how to read this format.

## Important Considerations

*   **K-mer Size Consistency:** When you plan to compare databases or query a database with reads, ensure you use the *same k-mer size* for the `build` command as you do for the `compare` or `query` commands.
*   **Canonical K-mers:** Orion-Kmer uses canonical k-mers. This means for each k-mer and its reverse complement, only the lexicographically smaller one is stored. This is handled automatically.
*   **Memory and Time:** For very large genomes, the `build` process can consume significant memory and time. The toy examples here will be very fast.

## Next Steps

Now that you have created k-mer databases, you can:
*   Compare them using `orion-kmer compare` ([Compare Tutorial](./compare_tutorial.md)).
*   Query reads against them using `orion-kmer query` ([Query Tutorial](./query_tutorial.md)).

---
[Back to Main Manual](../index.md) | [Previous Tutorial: Counting K-mers (`count`)](./count_tutorial.md) | [Next Tutorial: Comparing Databases (`compare`)](./compare_tutorial.md)
