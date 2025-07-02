# Orion-Kmer: A Blazing-Fast K-mer Toolkit in Rust

Orion-Kmer is a command-line bioinformatics tool designed for high-performance k-mer analysis. It leverages Rust's performance and memory safety to count k-mers, build k-mer databases, compare genomic datasets, and query reads against databases with speed and efficiency.

n**For a detailed user manual, tutorials, and full documentation, please visit the [Orion-Kmer GitHub Pages site](https://motroy.github.io/orion-kmer/).**

## Features

*   **Ultra-Fast K-mer Counting**: Counts k-mers from FASTA/FASTQ files.
*   **Compact K-mer Databases**: Creates binary databases of unique k-mers.
*   **Genomic Comparison**: Computes set operations (intersection, union) and Jaccard index between k-mer databases.
*   **Short Read Querying**: Screens short-read sequencing data (FASTQ) against a k-mer database.
*   **Multithreaded**: Utilizes available CPU cores for parallel processing (primarily for the `query` subcommand currently).
*   **Efficient Memory Usage**: Uses canonical k-mer representation.

## Installation

Currently, Orion-Kmer must be built from source using Cargo:

```bash
# Clone the repository (assuming it's in a public repo later)
# git clone <repo_url>
# cd orion-kmer
cargo build --release
# The binary will be at target/release/orion-kmer
```

## Usage

Orion-Kmer uses a subcommand interface, similar to `git` or `samtools`.

```
$ orion-kmer --help
Orion-Kmer v0.1.0
A blazing-fast k-mer toolkit in Rust

USAGE:
    orion-kmer [OPTIONS] <SUBCOMMAND>

OPTIONS:
    -h, --help       Print help information
    -t, --threads    Number of threads to use (0 for all logical cores) [default: 0]
    -V, --version    Print version information
    -v, --verbose    Verbosity level (e.g., -v, -vv)

SUBCOMMANDS:
    count      Count k-mers in FASTA/FASTQ files
    build      Build a unique k-mer database from genome assemblies
    compare    Compare two k-mer databases
    query      Query short reads against a k-mer database
    classify   Classify input sequences against k-mer databases and report coverage statistics
```

### Global Options

*   `-t, --threads <THREADS>`: Number of threads to use. Defaults to the number of logical cores if set to 0.
*   `-v, --verbose`: Increase verbosity (e.g., `-v` for info, `-vv` for debug).

### Subcommands

#### 1. `count`

Counts all k-mers in the input file(s) and outputs them in a simple text format (`kmer<TAB>count`).

**Usage:**

```bash
orion-kmer count -k <KMER_SIZE> -i <INPUT_FILE1> [INPUT_FILE2 ...] -o <OUTPUT_TSV> [-m <MIN_COUNT>]
```

**Arguments:**

*   `-k, --kmer-size <INT>`: The length of the k-mer (1-32) \[required].
*   `-i, --input <FILE>...`: One or more input FASTA/FASTQ files (can be gzipped) \[required].
*   `-o, --output <FILE>`: Output file for k-mer counts \[required].
*   `-m, --min-count <INT>`: Minimum count to report a k-mer \[default: 1].

**Example:**

```bash
orion-kmer count -k 21 -i genome.fasta -o counts.tsv
```

#### 2. `build`

Scans genome assemblies (FASTA) and creates a compact, binary database file (`.db`). This database stores the set of unique k-mers found in each input genome file separately, allowing for per-reference analysis by commands like `classify`.

**Usage:**

```bash
orion-kmer build -k <KMER_SIZE> -g <GENOME_FILE1> [GENOME_FILE2 ...] -o <OUTPUT_DB>
```

**Arguments:**

*   `-k, --kmer-size <INT>`: The length of the k-mer (1-32) \[required].
*   `-g, --genomes <FILE>...`: One or more input genome assembly files (FASTA format, can be gzipped) \[required].
*   `-o, --output <FILE>`: Output path for the binary k-mer database \[required].

**Example:**

```bash
orion-kmer build -k 31 -g genome1.fasta genome2.fasta -o unique_kmers.db
```

#### 3. `compare`

Compares two k-mer databases (created by `build`) and reports similarity statistics in JSON format.

**Usage:**

```bash
orion-kmer compare --db1 <DATABASE1_DB> --db2 <DATABASE2_DB> -o <OUTPUT_JSON>
```

**Arguments:**

*   `--db1 <FILE>`: First k-mer database file \[required].
*   `--db2 <FILE>`: Second k-mer database file \[required].
*   `-o, --output <FILE>`: Output file for comparison stats (JSON format) \[required].

**Example:**

```bash
orion-kmer compare --db1 e_coli.db --db2 salmonella.db -o comparison.json
```

**Output JSON Structure:**

```json
{
  "db1_path": "e_coli.db",
  "db2_path": "salmonella.db",
  "kmer_size": 31,
  "db1_unique_kmers": 4150234,
  "db2_unique_kmers": 4398102,
  "intersection_size": 3801299,
  "union_size": 4747037,
  "jaccard_index": 0.800773
}
```

#### 4. `query`

Takes a k-mer database and a short-read file (FASTQ) and finds reads containing k-mers present in the database. Outputs the IDs of matching reads.

**Usage:**

```bash
orion-kmer query -d <DATABASE_DB> -r <READS_FASTQ> -o <OUTPUT_READ_IDS> [-c <MIN_HITS>]
```

**Arguments:**

*   `-d, --database <FILE>`: K-mer database to query against \[required].
*   `-r, --reads <FILE>`: Short-read file (FASTQ, can be gzipped) \[required].
*   `-o, --output <FILE>`: Output file for the IDs of matching reads \[required].
*   `-c, --min-hits <INT>`: Minimum number of k-mer hits in a read to report it \[default: 1].

**Example:**

```bash
orion-kmer query -d reference_genome.db -r reads.fastq.gz -o matching_reads.txt -c 5
```

#### 5. `classify`

Compares an input FASTA/FASTQ file against one or more k-mer databases (built by `orion-kmer build`). It reports the proportion of k-mers from the input that match each database and each reference within those databases, along with coverage and depth statistics. The output is in JSON format.

**Usage:**

```bash
orion-kmer classify -i <INPUT_FASTA_OR_FASTQ> -d <DB1.db> [DB2.db ...] -o <OUTPUT_JSON> [--kmer-size <KMER_SIZE_VALIDATION>] [--min-kmer-frequency <MIN_FREQ>]
```

**Arguments:**

*   `-i, --input <FILE>`: Input genome (FASTA) or reads (FASTQ) file to classify \[required].
*   `-d, --databases <FILE>...`: One or more k-mer database files (`.db`) to classify against \[required].
*   `-o, --output <FILE>`: Output file for classification results (JSON format) \[required].
*   `--kmer-size <INT>`: Optional. If provided, this k-mer size is validated against the k-mer size stored in the databases. The command will fail if they don't match. If not provided, the k-mer size from the first database is used, and subsequent databases are validated against it.
*   `--min-kmer-frequency <INT>`: Minimum frequency for a k-mer in the input file to be considered for matching and depth calculations \[default: 1].
    *   `--min-coverage <FLOAT>`: Minimum reference breadth of coverage (proportion of reference k-mers found in input) to include a reference in the output JSON and TSV reports \[default: 0.0].
    *   `--output-tsv <FILE>`: Optional. Output file path for a TSV summary of the classification results. This summary includes one row per reference that meets the `--min-coverage` threshold.

**Example:**

```bash
orion-kmer classify -i my_reads.fastq -d ref_genome1.db pathogen_panel.db -o classification_report.json --min-coverage 0.05 --output-tsv classification_summary.tsv
```

**Output JSON Structure Example:**

```json
{
  "input_file_path": "my_reads.fastq",
  "total_unique_kmers_in_input": 150000, // After min_kmer_frequency filter
  "min_kmer_frequency_filter": 1,
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
      // ... results for pathogen_panel.db ...
    }
  ]
}
```

## Example Workflow

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

5.  **Classify reads against multiple databases (e.g., the reference and the new assembly):**
    ```bash
    orion-kmer classify -i experiment.fastq.gz -d ref.db new.db -o read_classification.json -t 8
    ```

## Technical Details

*   K-mers (up to k=32) are encoded into `u64` integers.
*   Canonical k-mers (lexicographically smaller of a k-mer and its reverse complement) are used throughout.
*   Uses `needletail` for FASTA/FASTQ parsing, `clap` for CLI, `rayon` for parallelism, `dashmap` for concurrent collections, `serde` and `bincode` for serialization.
```
