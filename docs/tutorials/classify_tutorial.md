---
theme: default
title: Classifying Sequences (orion-kmer classify)
---

# Classifying Sequences with `orion-kmer classify`

The `orion-kmer classify` subcommand compares an input FASTA/FASTQ file against one or more k-mer databases (built by `orion-kmer build`). It reports the proportion of k-mers from the input that match each database and each reference within those databases, along with coverage and depth statistics. The output is in JSON format.

This tutorial will guide you through the usage of the `classify` command.

## Prerequisites

*   Orion-Kmer installed.
*   One or more k-mer databases created using `orion-kmer build`. For example, `ref_genome1.db` and `pathogen_panel.db`.
*   An input file (FASTA or FASTQ) containing sequences to classify, for example, `my_reads.fastq`.

## Command Usage

The basic syntax for `orion-kmer classify` is:

```bash
orion-kmer classify -i <INPUT_FASTA_OR_FASTQ> -d <DB1.db> [DB2.db ...] -o <OUTPUT_JSON> [OPTIONS]
```

### Arguments and Options:

*   `-i, --input <FILE>`: **Required.** Input genome (FASTA) or reads (FASTQ) file to classify.
*   `-d, --databases <FILE>...`: **Required.** One or more k-mer database files (`.db`) to classify against.
*   `-o, --output <FILE>`: **Required.** Output file for classification results (JSON format).
*   `--kmer-size <INT>`: Optional. If provided, this k-mer size is validated against the k-mer size stored in the databases. The command will fail if they don't match. If not provided, the k-mer size from the first database is used, and subsequent databases are validated against it.
*   `--min-kmer-frequency <INT>`: Minimum frequency for a k-mer in the input file to be considered for matching and depth calculations. Default: `1`.
*   Global options like `-t, --threads` and `-v, --verbose` are also applicable.

## Example

Let's say you have a set of sequencing reads in `my_reads.fastq` and you want to classify them against two databases: `ref_genome1.db` (a reference genome) and `pathogen_panel.db` (a database of pathogen k-mers).

```bash
orion-kmer classify \
    -i my_reads.fastq \
    -d ref_genome1.db pathogen_panel.db \
    -o classification_report.json \
    --min-kmer-frequency 2 \
    -t 8
```

This command will:
1.  Read k-mers from `my_reads.fastq`.
2.  Consider only those k-mers that appear at least 2 times (`--min-kmer-frequency 2`).
3.  Compare these filtered k-mers against `ref_genome1.db` and `pathogen_panel.db`.
4.  Write the classification statistics to `classification_report.json`.
5.  Use 8 threads for processing.

## Output JSON Structure

The output JSON file (`classification_report.json` in the example) will have a structure similar to this:

```json
{
  "input_file_path": "my_reads.fastq",
  "total_unique_kmers_in_input": 150000, // After min_kmer_frequency filter
  "min_kmer_frequency_filter": 2,        // Value used in the command
  "databases_analyzed": [
    {
      "database_path": "ref_genome1.db",
      "database_kmer_size": 31,
      "total_unique_kmers_in_db_across_references": 4500000,
      "overall_input_kmers_matched_in_db": 75000,
      "overall_sum_depth_of_matched_kmers_in_input": 300000,
      "overall_avg_depth_of_matched_kmers_in_input": 4.0,
      "proportion_input_kmers_in_db_overall": 0.5, // 75000 / 150000
      "proportion_db_kmers_covered_overall": 0.01666, // 75000 / 4500000
      "references": [
        {
          "reference_name": "ref_genome1_contig1.fa", // Filename used during 'build'
          "total_kmers_in_reference": 2000000,
          "input_kmers_hitting_reference": 60000,
          "sum_depth_of_matched_kmers_in_input": 250000,
          "avg_depth_of_matched_kmers_in_input": 4.1667,
          "proportion_input_kmers_hitting_reference": 0.4, // 60000 / 150000
          "reference_breadth_of_coverage": 0.03 // 60000 / 2000000
        },
        {
          "reference_name": "ref_genome1_contig2.fa",
          "total_kmers_in_reference": 2500000,
          "input_kmers_hitting_reference": 15000,
          "sum_depth_of_matched_kmers_in_input": 50000,
          "avg_depth_of_matched_kmers_in_input": 3.3333,
          "proportion_input_kmers_hitting_reference": 0.1, // 15000 / 150000
          "reference_breadth_of_coverage": 0.006 // 15000 / 2500000
        }
      ]
    },
    {
      // ... results for pathogen_panel.db would follow a similar structure ...
      "database_path": "pathogen_panel.db",
      // ...
      "references": [
        // ...
      ]
    }
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
    *   The `references` array provides these metrics broken down for each individual reference sequence that was part of the database build.
        *   **`reference_breadth_of_coverage`**: (`input_kmers_hitting_reference` / `total_kmers_in_reference`). What fraction of a specific reference's k-mers were hit by your input.

This information helps you understand how well your input sequences are represented in each database and in specific genomes within those databases.
---

*Next, you might want to check out other tutorials or head back to the [main documentation page](../index.md).*
