---
layout: default
title: Compare K-mer Databases Tutorial
---

# Tutorial: Comparing K-mer Databases with `orion-kmer compare`

[Back to Main Manual](../index.md)

The `orion-kmer compare` subcommand takes two k-mer databases (created by `orion-kmer build`) and calculates similarity statistics between them. These statistics include the number of unique k-mers in each database, the size of their intersection (common k-mers), the size of their union (total unique k-mers across both), and the Jaccard index. The output is in JSON format.

## Command Overview

```bash
orion-kmer compare --db1 <DATABASE1_DB> --db2 <DATABASE2_DB> -o <OUTPUT_JSON>
```

### Arguments:

*   `--db1 <FILE>`: **Required.** Path to the first k-mer database file (e.g., `genome1.db`).
*   `--db2 <FILE>`: **Required.** Path to the second k-mer database file (e.g., `genome2.db`).
*   `-o, --output <FILE>`: **Required.** Path to the output file where the comparison statistics (in JSON format) will be written (e.g., `comparison_stats.json`).

## Tutorial with Toy Data

In the [Build K-mer Database Tutorial](./build_tutorial.md), we created `genome1.db` and `genome2.db` from our toy FASTA files. We will now compare these two databases.

**Prerequisites:**
*   `output/genome1.db` (created from `data/genome1.fasta` with k=7)
*   `output/genome2.db` (created from `data/genome2.fasta` with k=7)

Ensure `orion-kmer` is accessible and you are in the `docs/tutorials` directory or adjust paths.

### Step 1: Comparing `genome1.db` and `genome2.db`

```bash
# Assuming orion-kmer is in your PATH
# and you are in the 'docs/tutorials' directory
orion-kmer compare --db1 output/genome1.db --db2 output/genome2.db -o output/comparison_g1_g2.json
```

This command will:
1.  Load the k-mer sets from `output/genome1.db` and `output/genome2.db`.
2.  Calculate the intersection, union, and Jaccard index.
3.  Write these statistics to `output/comparison_g1_g2.json`.

### Step 2: Examining the Output JSON

The output file `output/comparison_g1_g2.json` will contain something like this (exact numbers will depend on the k-mer content of `genome1.fasta` and `genome2.fasta` with k=7):

```json
{
  "db1_path": "output/genome1.db",
  "db2_path": "output/genome2.db",
  "kmer_size": 7,
  "db1_unique_kmers": 56,   // Example value
  "db2_unique_kmers": 56,   // Example value
  "intersection_size": 0, // Example: if no k-mers are shared (due to ATGC vs TACG nature of toy data)
  "union_size": 112,        // Example: db1_unique + db2_unique - intersection
  "jaccard_index": 0.0    // Example: intersection_size / union_size
}
```

**Explanation of fields:**

*   `db1_path`: Path to the first database.
*   `db2_path`: Path to the second database.
*   `kmer_size`: The k-mer size used when building these databases. The `compare` command infers this from the database files. **It's crucial that both databases were built with the same k-mer size.**
*   `db1_unique_kmers`: Number of unique k-mers in the first database.
*   `db2_unique_kmers`: Number of unique k-mers in the second database.
*   `intersection_size`: Number of k-mers common to both databases.
*   `union_size`: Total number of unique k-mers present in either database. Calculated as `db1_unique_kmers + db2_unique_kmers - intersection_size`.
*   `jaccard_index`: A measure of similarity, calculated as `intersection_size / union_size`. Ranges from 0 (no similarity) to 1 (identical k-mer sets).

### Manual Verification (Conceptual for k=3)

Let's consider a very simplified example with k=3 for conceptual understanding.
*   `genome1.fasta` starts: `ATGCGTAGC...`
    *   Unique 3-mers (canonical): `ATG`, `CGC`, `CTA`, `GCA`, `GCG`, `GTA`, `TAG` ...
*   `genome2.fasta` starts: `TACGCTAGC...`
    *   Unique 3-mers (canonical): `ACG`, `AGC`, `CGC`, `CTA`, `GCT`, `TAG`, `TAC` ...

If these were our full genomes:
*   `db1_unique_kmers`: (count of unique canonical 3-mers from genome1)
*   `db2_unique_kmers`: (count of unique canonical 3-mers from genome2)
*   `intersection_size`: `CGC`, `CTA`, `TAG` (shared ones) -> 3
*   `union_size`: (`db1_unique` + `db2_unique` - 3)
*   `jaccard_index`: `3 / union_size`

The actual toy data is longer and uses k=7, so the numbers will be different, but the principle is the same. Our toy data `genome1.fasta` and `genome2.fasta` are designed to be mostly reverse complements of each other for the first parts of their sequences. Given canonical k-mer representation, this means many k-mers might actually be identical.

For `genome1_seq1` (`ATGCGTAGCATCGATCG...`) and `genome2_seq1` (`TACGCTAGCTAGCTA...`), if k=7:
- `genome1`: `ATGCGTA` (revcomp `TACGCAT`) -> `ATGCGTA`
- `genome2`: `TACGCTA` (revcomp `TAGCGTA`) -> `TACGCTA`
These are different.

Let's look at the repetitive parts:
`genome1_seq1`: `...CATCGATCGATCGATCGA...`
  - `CATCGAT` (revcomp `ATCGATG`)
  - `ATCGATC` (revcomp `GATCGAT`)
  - `TCGATCG` (revcomp `CGATCGA`)
  - `CGATCGA` (revcomp `TCGATCG`)
  - `GATCGAT` (revcomp `ATCGATC`)
`genome2_seq2`: `GCTAGCTAGCTAGCTAG...`
  - `GCTAGCT` (revcomp `AGCTAGC`)
  - `CTAGCTA` (revcomp `TAGCTAG`)
  - `TAGCTAG` (revcomp `CTAGCTA`)
  - `AGCTAGC` (revcomp `GCTAGCT`)

It seems for the toy data, the intersection might be small or zero for k=7, especially between `genome1_seq1` and `genome2_seq1` if one is based on `ATGC...` and the other `TACG...`. The repetitive `CGATCG...` and `GCTAGC...` parts are also not direct reverse complements.
This means a Jaccard index close to 0 for k=7 is expected for the provided toy data.

## Important Considerations

*   **Matching K-mer Sizes:** The `compare` command will fail or produce meaningless results if the two databases were built with different k-mer sizes. Orion-Kmer stores the k-mer size in the database and should warn you.
*   **Interpretation:** The Jaccard index is a common measure for genome similarity based on k-mer content. A higher Jaccard index indicates more shared k-mers and thus higher similarity.

## Next Steps

*   Try comparing `genome1.db` with `combined_genomes.db` (if you created it in the `build` tutorial). What would you expect the Jaccard index to be?
*   Proceed to the [Query Tutorial](./query_tutorial.md) to learn how to use these databases to find reads.

---
[Back to Main Manual](../index.md) | [Previous Tutorial: Building Databases (`build`)](./build_tutorial.md) | [Next Tutorial: Querying Reads (`query`)](./query_tutorial.md)
