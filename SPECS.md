# Orion-Kmer: A Blazing-Fast K-mer Toolkit in Rust
Orion-Kmer is a command-line bioinformatics tool designed for high-performance k-mer analysis. It leverages Rust's performance, memory safety, and first-class concurrency to count k-mers, build databases, and compare genomic datasets with maximum speed and efficiency.
## Core Features
 * Ultra-Fast K-mer Counting: Utilizes optimized data structures and multithreading to count k-mers from FASTA and FASTQ files at near I/O-bound speeds.
 * Compact K-mer Database Creation: Generates a memory-efficient, binary database of unique k-mers from one or more genome assemblies.
 * Genomic K-mer Comparison: Computes set operations (intersection, union, difference) and similarity metrics (e.g., Jaccard index) between two or more k-mer databases.
 * Short Read Querying: Rapidly screens short-read sequencing data (FASTQ) against a pre-built k-mer database to identify matching reads.
 * Multithreaded by Design: Scales its performance across all available CPU cores.
 * Efficient Memory Usage: Employs canonical k-mer representation and memory-optimized data structures.
## Under the Hood: Key Design Choices for Speed
Orion-Kmer's performance comes from several key design decisions:
 * K-mer Representation: Instead of storing k-mers as strings, they are encoded into a 64-bit unsigned integer (u64). Each DNA base (A, C, G, T) is represented by 2 bits. This allows for storing k-mers up to a length of k=32 in a single machine word, enabling extremely fast comparisons and hashing.
   * A \\rightarrow 00
   * C \\rightarrow 01
   * G \\rightarrow 10
   * T \\rightarrow 11
 * Canonical K-mers: To handle DNA's double-stranded nature, Orion-Kmer always stores the canonical k-mer. For any given k-mer, it also computes its reverse complement and stores the lexicographically smaller of the two. This ensures that AGCT and its reverse complement AGCT (from the other strand) are treated as the same entity, which is standard practice in genomics.
 * Core Data Structure: The heart of the tool is a high-performance, concurrent hash map. Instead of a standard library HashMap protected by a mutex (which causes lock contention), it will use a crate like dashmap or chashmap. These allow multiple threads to read and write to the hash map simultaneously with minimal blocking, which is critical for scaling.
 * Parallel Processing Model: Orion-Kmer uses a producer-consumer architecture.
   * I/O Thread(s): One or more dedicated threads read the FASTA/FASTQ file(s) in large chunks.
   * Worker Threads: A pool of worker threads (configurable via a command-line flag) receives these chunks. Each worker independently finds and encodes canonical k-mers from its chunk and inserts or updates their counts in the global concurrent hash map. The rayon crate is perfect for this, allowing for easy parallel iteration over sequence chunks.
 * Optimized I/O: It uses the needletail crate, the de facto standard for fast and robust FASTA/FASTQ parsing in Rust. It handles various file formats and edge cases correctly and efficiently.
## Command-Line Interface (CLI) Design
The tool would be structured with subcommands, similar to git or samtools.
```
$ orion-kmer --help
Orion-Kmer v0.1.0
A blazing-fast k-mer toolkit in Rust

USAGE:
    orion-kmer [OPTIONS] <SUBCOMMAND>

OPTIONS:
    -h, --help       Print help information
    -t, --threads    Number of threads to use [default: number of logical cores]
    -V, --version    Print version information

SUBCOMMANDS:
    count      Count k-mers in FASTA/FASTQ files
    build      Build a unique k-mer database from genome assemblies
    compare    Compare two k-mer databases
    query      Query short reads against a k-mer database
```

## 1. count Subcommand
Counts all k-mers in the input file(s) and outputs them in a simple text format.
### Usage
```
orion-kmer count -k 21 -i genome.fasta -o counts.tsv
```
### Arguments
```
-k, --kmer-size <INT>    The length of the k-mer [required]
-i, --input <FILE>...    One or more input FASTA/FASTQ files (can be gzipped)
-o, --output <FILE>      Output file for k-mer counts (kmer<TAB>count)
-m, --min-count <INT>    Minimum count to report a k-mer [default: 1]
```
## 2. build Subcommand
Scans genome assemblies and creates a compact, binary file containing only the set of unique k-mers.
### Usage
```
orion-kmer build -k 31 -g genome1.fasta genome2.fasta -o unique_kmers.db
```
### Arguments
```
-k, --kmer-size <INT>    The length of the k-mer [required]
-g, --genomes <FILE>...  One or more input genome assembly files (FASTA)
-o, --output <FILE>      Output path for the binary k-mer database [required]
```
The output .db file would be a serialized HashSet<u64> for maximum lookup speed.
## 3. compare Subcommand
Compares two k-mer databases and reports similarity statistics.
### Usage
```
orion-kmer compare --db1 e_coli.db --db2 salmonella.db -o comparison.json
```
### Arguments
```
--db1 <FILE>      First k-mer database file [required]
--db2 <FILE>      Second k-mer database file [required]
-o, --output <FILE> Output file for comparison stats (JSON format)
```
The output comparison.json would look like:
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

Where the Jaccard index is calculated as J(A, B) = \\frac{|A \\cap B|}{|A \\cup B|}.
## 4. query Subcommand
Takes a k-mer database and a short-read file (FASTQ) and efficiently finds which reads contain k-mers present in the database.
### Usage
```
orion-kmer query -d reference_genome.db -r reads.fastq.gz -o matching_reads.txt
```
### Arguments
```
-d, --database <FILE>    K-mer database to query against [required]
-r, --reads <FILE>       Short-read file (FASTQ, can be gzipped) [required]
-o, --output <FILE>      Output file for the IDs of matching reads
-c, --min-hits <INT>     Minimum number of k-mer hits to report a read [default: 1]
```
## Example Workflow
 * Build a database for a reference genome:
   orion-kmer build -k 31 -t 16 -g reference.fasta -o ref.db

 * Build a database for a newly assembled genome:
   orion-kmer build -k 31 -t 16 -g new_assembly.fasta -o new.db

 * Compare the two assemblies:
   orion-kmer compare --db1 ref.db --db2 new.db -o report.json

 * Check which reads from a sequencing experiment belong to the reference:
   orion-kmer query -d ref.db -r experiment.fastq.gz -o matched_reads.txt -c 5 -t 16

Key Rust Dependencies (Crates)
 * clap: For parsing command-line arguments.
 * needletail: For fast and reliable FASTA/FASTQ parsing.
 * rayon: For high-level, data-parallelism.
 * dashmap or chashmap: For the high-performance concurrent hash map.
 * serde: For serializing/deserializing the k-mer databases and JSON output.
 * flate2: For transparently handling gzipped input files.
This design provides a complete, robust, and exceptionally fast tool that meets all your specified requirements, using the best practices of modern systems programming in Rust.
