---
layout: default
title: Classifying Sequences (orion-kmer classify)
---

# Classifying Sequences with `orion-kmer classify`

The `orion-kmer classify` subcommand compares an input FASTA/FASTQ file against one or more k-mer databases (built by `orion-kmer build`). It reports the proportion of k-mers from the input that match each database and each reference within those databases, along with coverage and depth statistics. The output is in JSON format.

This tutorial will guide you through the usage of the `classify` command.

## Prerequisites

*   Orion-Kmer installed.
*   One or more k-mer databases created using `orion-kmer build`. For example, `ref_genome1.db` and `pathogen_panel.db`.
*   An input file (FASTA or FASTQ) containing sequences to classify, for example, `my_reads.fastq`.
*   Ensure `genome1.db` (from `build_tutorial.md`) and `genome2.db` are available in `output/` if you want to follow the example directly.

## Command Usage

The basic syntax for `orion-kmer classify` is:

```bash
orion-kmer classify -i <INPUT_FASTA_OR_FASTQ> -d <DB1.db> [DB2.db ...] -o <OUTPUT_JSON> [OPTIONS]
```

### Arguments and Options:

*   `-i, --input <FILE>`: **Required.** Input genome (FASTA) or reads (FASTQ) file to classify.
*   `-d, --databases <FILE>...`: **Required.** One or more k-mer database files (`.db`) to classify against.
*   `-o, --output <FILE>`: **Required.** Output file for classification results (JSON format).
*   `--kmer-size <INT>`: Optional. If provided, this k-mer size is validated against the k-mer size stored in the databases. The command will fail if they don't match. If not provided, the k-mer size from the first database is used, and subsequent databases are validated against it. (Assumed k=7 for examples if databases from previous tutorials are used).
*   `--min-kmer-frequency <INT>`: Minimum frequency for a k-mer in the input file to be considered for matching and depth calculations. Default: `1`.
*   Global options like `-t, --threads` and `-v, --verbose` are also applicable.

## Example

Let's use `data/reads.fastq` as our input and classify it against `output/genome1.db` and `output/genome2.db` (which should have been built with k=7 from previous tutorials).

```bash
# Assuming orion-kmer is in your PATH
# and you are in the 'docs/tutorials' directory
orion-kmer classify \
    -i data/reads.fastq \
    -d output/genome1.db output/genome2.db \
    -o output/classify_reads_report.json \
    --min-kmer-frequency 1 \
    -t 4 # Optional: use threads
```

This command will:
1.  Read k-mers from `data/reads.fastq`.
2.  Consider k-mers that appear at least once (`--min-kmer-frequency 1`). The k-mer size will be inferred from `output/genome1.db` (expected to be 7).
3.  Compare these k-mers against `output/genome1.db` and `output/genome2.db`.
4.  Write the classification statistics to `output/classify_reads_report.json`.

## Output JSON Structure

The output JSON file (`output/classify_reads_report.json` in the example) will have a structure similar to this (exact numbers depend on k-mer size, input data, and database contents):

```json
{
  "input_file_path": "data/reads.fastq",
  "total_unique_kmers_in_input": 63, // Example: Actual number of unique 7-mers in reads.fastq
  "min_kmer_frequency_filter": 1,
  "databases_analyzed": [
    {
      "database_path": "output/genome1.db",
      "database_kmer_size": 7,
      "total_unique_kmers_in_db_across_references": 56, // From genome1.db (k=7)
      "overall_input_kmers_matched_in_db": 38,         // Example: k-mers from reads found in genome1.db
      "overall_sum_depth_of_matched_kmers_in_input": 38, // Assuming each k-mer appears once in reads
      "overall_avg_depth_of_matched_kmers_in_input": 1.0,
      "proportion_input_kmers_in_db_overall": 0.603,   // 38 / 63
      "proportion_db_kmers_covered_overall": 0.678,    // 38 / 56
      "references": [ // genome1.db was built from a single file data/genome1.fasta
        {
          "reference_name": "data/genome1.fasta",
          "total_kmers_in_reference": 56,
          "input_kmers_hitting_reference": 38,
          "sum_depth_of_matched_kmers_in_input": 38,
          "avg_depth_of_matched_kmers_in_input": 1.0,
          "proportion_input_kmers_hitting_reference": 0.603, // 38 / 63
          "reference_breadth_of_coverage": 0.678        // 38 / 56
        }
      ]
    },
    {
      "database_path": "output/genome2.db",
      "database_kmer_size": 7,
      "total_unique_kmers_in_db_across_references": 56, // From genome2.db (k=7)
      "overall_input_kmers_matched_in_db": 38,         // Example: k-mers from reads found in genome2.db
      "overall_sum_depth_of_matched_kmers_in_input": 38,
      "overall_avg_depth_of_matched_kmers_in_input": 1.0,
      "proportion_input_kmers_in_db_overall": 0.603,   // 38 / 63
      "proportion_db_kmers_covered_overall": 0.678,    // 38 / 56
      "references": [ // genome2.db was built from a single file data/genome2.fasta
        {
          "reference_name": "data/genome2.fasta",
          "total_kmers_in_reference": 56,
          "input_kmers_hitting_reference": 38,
          "sum_depth_of_matched_kmers_in_input": 38,
          "avg_depth_of_matched_kmers_in_input": 1.0,
          "proportion_input_kmers_hitting_reference": 0.603,
          "reference_breadth_of_coverage": 0.678
        }
      ]
    }
    // Note: The example numbers above for matches (38) for both databases would imply reads.fastq has k-mers
    // that hit both genome1 and genome2. Read1 & Read2 hit genome1.db. Read3 & Read4 hit genome2.db.
    // Read5 (AAAAAAA...) hits neither. So the actual numbers will be different.
    // A more realistic scenario:
    // genome1.db: overall_input_kmers_matched_in_db might be around 19 (from read1) + 19 (from read2) = 38 (if no overlap)
    // genome2.db: overall_input_kmers_matched_in_db might be around 19 (from read3) + 19 (from read4) = 38 (if no overlap)
    // total_unique_kmers_in_input (from reads.fastq, k=7):
    //   read1 (ATGCGTAGCATCGATCGATCGATCG, len 25, 19 kmers)
    //   read2 (CGATCGATCGATCGATCGATCGATC, len 25, 19 kmers) - some overlap with read1's end
    //   read3 (TACGCTAGCTAGCTAGCTAGCTAGC, len 25, 19 kmers)
    //   read4 (GCTAGCTAGCTAGCTAGCTAGCTAG, len 25, 19 kmers) - some overlap with read3's end
    //   read5 (AAAAAAAAAAAAAAAAAAAAAAAAA, len 25, 1 kmer type: AAAAAAA)
    // The exact unique count needs careful calculation.
  ]
}
```

### Key Metrics in the Output:

*   **`total_unique_kmers_in_input`**: The number of distinct k-mers in your input file after the `--min-kmer-frequency` filter is applied.
*   For each database (`databases_analyzed` array):
    *   **`overall_input_kmers_matched_in_db`**: How many unique k-mers from your input were found in this database.
    *   **`overall_sum_depth_of_matched_kmers_in_input`**: The total count (frequency) in the input file of all k-mers that matched this database.
    *   **`overall_avg_depth_of_matched_kmers_in_input`**: Average frequency in the input of the k-mers that matched.
    *   **`proportion_input_kmers_in_db_overall`**: (`overall_input_kmers_matched_in_db` / `total_unique_kmers_in_input`). What fraction of your input's k-mers hit this database.
    *   **`proportion_db_kmers_covered_overall`**: (`overall_input_kmers_matched_in_db` / `total_unique_kmers_in_db_across_references`). What fraction of this database's k-mers were hit by your input.
    *   The `references` array provides these metrics broken down for each individual reference sequence that was part of the database build (e.g., each FASTA file given to `orion-kmer build`).
        *   **`reference_breadth_of_coverage`**: (`input_kmers_hitting_reference` / `total_kmers_in_reference`). What fraction of a specific reference's k-mers were hit by your input.

This information helps you understand how well your input sequences are represented in each database and in specific genomes within those databases.
---

*Next, you might want to check out other tutorials or head back to the [main documentation page](../index.md).*
