---
layout: default
title: Query Reads Tutorial
---

# Tutorial: Querying Reads with `orion-kmer query`

[Back to Main Manual](../index.md)

The `orion-kmer query` subcommand screens a short-read file (FASTQ format) against a k-mer database (created by `orion-kmer build`). It identifies and outputs the IDs of reads that contain a certain number of k-mers also found in the database. This is useful for quickly determining if reads from a sequencing experiment likely originate from a genome or set of genomes represented in the database.

## Command Overview

```bash
orion-kmer query -d <DATABASE_DB> -r <READS_FASTQ> -o <OUTPUT_READ_IDS> [-c <MIN_HITS>]
```

### Arguments:

*   `-d, --database <FILE>`: **Required.** Path to the k-mer database file (e.g., `genome1.db`) to query against.
*   `-r, --reads <FILE>`: **Required.** Path to the short-read file in FASTQ format (can be gzipped, e.g., `my_reads.fastq.gz`).
*   `-o, --output <FILE>`: **Required.** Path to the output file where the IDs of matching reads will be written (one read ID per line).
*   `-c, --min-hits <INT>`: Optional. The minimum number of k-mers from a read that must be present in the database for that read to be considered a match and its ID reported. Defaults to `1`.

## Tutorial with Toy Data

We will use the `genome1.db` created in the [Build Tutorial](./build_tutorial.md) (using k=7) and the toy FASTQ file `reads.fastq`.

**Prerequisites:**
*   `output/genome1.db` (created from `data/genome1.fasta` with k=7)
*   `data/reads.fastq`:
    ```fastq
    @read1
    ATGCGTAGCATCGATCGATCGATCG
    +
    IIIIIIIIIIIIIIIIIIIIIIIII
    @read2
    CGATCGATCGATCGATCGATCGATC
    +
    IIIIIIIIIIIIIIIIIIIIIIIII
    @read3
    TACGCTAGCTAGCTAGCTAGCTAGC
    +
    IIIIIIIIIIIIIIIIIIIIIIIII
    @read4
    GCTAGCTAGCTAGCTAGCTAGCTAG
    +
    IIIIIIIIIIIIIIIIIIIIIIIII
    @read5
    AAAAAAAAAAAAAAAAAAAAAAAAA
    +
    IIIIIIIIIIIIIIIIIIIIIIIII
    ```

Ensure `orion-kmer` is accessible and you are in the `docs/tutorials` directory or adjust paths.

### Step 1: Querying `reads.fastq` Against `genome1.db`

Let's find reads in `data/reads.fastq` that have at least one 7-mer present in `output/genome1.db`.

```bash
# Assuming orion-kmer is in your PATH
# and you are in the 'docs/tutorials' directory
orion-kmer query -d output/genome1.db -r data/reads.fastq -o output/matching_reads_g1.txt
```

This command will:
1.  Load the k-mer set from `output/genome1.db` (which contains 7-mers from `genome1.fasta`).
2.  For each read in `data/reads.fastq`:
    a.  Extract all 7-mers from the read sequence.
    b.  Count how many of these read 7-mers are found in `genome1.db`.
3.  If this count is `>= 1` (the default for `--min-hits`), the read ID (e.g., `@read1`) is written to `output/matching_reads_g1.txt`.

### Step 2: Examining the Output

The file `output/matching_reads_g1.txt` will contain the IDs of the matching reads.
Given `genome1.fasta` starts `ATGCGTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCG` and `CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC`.
And `reads.fastq` contains:
- `@read1: ATGCGTAGCATCGATCGATCGATCG` (25 bases)
- `@read2: CGATCGATCGATCGATCGATCGATC` (25 bases)
- `@read3: TACGCTAGCTAGCTAGCTAGCTAGC` (25 bases)
- `@read4: GCTAGCTAGCTAGCTAGCTAGCTAG` (25 bases)
- `@read5: AAAAAAAAAAAAAAAAAAAAAAAAA` (25 bases)

The database `genome1.db` was built with k=7.
- **Read1:** `ATGCGTAGCATCGATCGATCGATCG`. Many 7-mers from this read (e.g., `ATGCGTA`, `TGCGTGA`, `GCGTACG`, ..., `ATCGATC`, `TCGATCG`) should be present in `genome1.db` because `read1` is a prefix of `genome1_seq1`.
- **Read2:** `CGATCGATCGATCGATCGATCGATC`. Similarly, many 7-mers from this read should be present in `genome1.db` as it's derived from the repetitive parts of `genome1.fasta`.
- **Read3:** `TACGCTAGCTAGCTAGCTAGCTAGC`. This sequence is more similar to `genome2.fasta`. It's unlikely its 7-mers (e.g., `TACGCTA`, `ACGCTAG`) are in `genome1.db` unless there are accidental overlaps or reverse complements that match. Given our toy data design, it's less likely.
- **Read4:** `GCTAGCTAGCTAGCTAGCTAGCTAG`. Similar to `read3`, less likely to match `genome1.db`.
- **Read5:** `AAAAAAAAAAAAAAAAAAAAAAAAA`. The 7-mer `AAAAAAA` (and its revcomp `TTTTTTT`) would need to be in `genome1.db`. It's not in our defined `genome1.fasta`.

So, the expected output in `output/matching_reads_g1.txt` would be:
```
@read1
@read2
```
*(The exact output depends on the k-mer content and canonicalization. Order is not guaranteed.)*

### Step 3: Using `--min-hits`

Let's try to be more stringent. Suppose we only want reads where at least 5 of their 7-mers are found in the database.
The length of reads is 25 bp. The number of 7-mers in a 25bp read is `25 - 7 + 1 = 19`.

```bash
orion-kmer query -d output/genome1.db -r data/reads.fastq -o output/matching_reads_g1_min5.txt -c 5
```
The output `output/matching_reads_g1_min5.txt` might still contain `@read1` and `@read2` if they indeed have at least 5 k-mer hits. If `read1` is almost a perfect match to the start of `genome1_seq1`, it will likely have many more than 5 hits. Same for `read2`. `read3`, `read4`, and `read5` are unlikely to meet this threshold.

The expected output in `output/matching_reads_g1_min5.txt` would likely be the same:
```
@read1
@read2
```

### Step 4: Querying Against a Different Database

If you also created `output/genome2.db` (from `genome2.fasta` which starts `TACGCTAGC...` and `GCTAGCT...`), you could query against it:
```bash
orion-kmer query -d output/genome2.db -r data/reads.fastq -o output/matching_reads_g2.txt -c 1
```
In this case, you would expect `@read3` and `@read4` to be in `output/matching_reads_g2.txt`.

Expected `output/matching_reads_g2.txt`:
```
@read3
@read4
```

## Important Considerations

*   **K-mer Size:** The k-mer size used to build the database (`-k` in `orion-kmer build`) is implicitly used for querying. The reads will be broken down into k-mers of that same size.
*   **FASTQ Format:** Ensure your read file is in valid FASTQ format.
*   **Performance:** `orion-kmer query` is designed to be fast, especially for large databases and read sets. It uses multithreading if the `-t` global option is specified.
*   **Interpreting `--min-hits`:** A higher `--min-hits` value increases specificity (fewer false positives) but may reduce sensitivity (miss reads with fewer true k-mer matches, e.g., due to sequencing errors or variations). The optimal value can depend on read length, k-mer size, and error rates.

## Conclusion

The `query` subcommand is a powerful tool for rapidly screening sequencing reads against existing genomic data. By adjusting the `--min-hits` parameter, you can control the stringency of the matching.

---
[Back to Main Manual](../index.md) | [Previous Tutorial: Comparing Databases (`compare`)](./compare_tutorial.md)
