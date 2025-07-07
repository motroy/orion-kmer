use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
#[clap(propagate_version = true)]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Commands,

    #[clap(
        short,
        long,
        global = true,
        default_value_t = 0,
        help = "Number of threads to use (0 for all logical cores)"
    )]
    pub threads: usize,

    #[clap(short, long, global = true, action = clap::ArgAction::Count, help = "Verbosity level (e.g., -v, -vv)")]
    pub verbose: u8,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Count k-mers in FASTA/FASTQ files
    Count(CountArgs),
    /// Build a unique k-mer database from genome assemblies
    Build(BuildArgs),
    /// Compare two k-mer databases
    Compare(CompareArgs),
    /// Query short reads against a k-mer database
    Query(QueryArgs),
    /// Classify sequences against k-mer databases and report coverage statistics
    Classify(ClassifyArgs),
}

#[derive(Parser, Debug)]
pub struct CountArgs {
    #[clap(short, long, required = true, help = "The length of the k-mer")]
    pub kmer_size: u8,

    #[clap(short, long, required = true, num_args = 1.., help = "One or more input FASTA/FASTQ files. Supports .gz, .xz, .zst, .7z compression.")]
    pub input_files: Vec<PathBuf>,

    #[clap(
        short,
        long,
        required = true,
        help = "Output file for k-mer counts (kmer<TAB>count). Supports .gz, .xz, .zst, .7z compression based on extension."
    )]
    pub output_file: PathBuf,

    #[clap(
        short = 'm',
        long,
        default_value_t = 1,
        help = "Minimum count to report a k-mer"
    )]
    pub min_count: usize,
}

#[derive(Parser, Debug)]
pub struct BuildArgs {
    #[clap(short, long, required = true, help = "The length of the k-mer")]
    pub kmer_size: u8,

    #[clap(short = 'g', long = "genomes", required = true, num_args = 1.., help = "One or more input genome assembly files (FASTA). Supports .gz, .xz, .zst, .7z compression.")]
    pub genome_files: Vec<PathBuf>,

    #[clap(
        short,
        long,
        required = true,
        help = "Output path for the binary k-mer database. Supports .gz, .xz, .zst, .7z compression based on extension."
    )]
    pub output_file: PathBuf,
}

#[derive(Parser, Debug)]
pub struct CompareArgs {
    #[clap(long, required = true, help = "First k-mer database file. Supports .gz, .xz, .zst, .7z compression.")]
    pub db1: PathBuf,

    #[clap(long, required = true, help = "Second k-mer database file. Supports .gz, .xz, .zst, .7z compression.")]
    pub db2: PathBuf,

    #[clap(
        short,
        long,
        required = true,
        help = "Output file for comparison stats (JSON format). Supports .gz, .xz, .zst, .7z compression based on extension."
    )]
    pub output_file: PathBuf,
}

#[derive(Parser, Debug)]
pub struct QueryArgs {
    #[clap(
        short = 'd',
        long = "database",
        required = true,
        help = "K-mer database to query against. Supports .gz, .xz, .zst, .7z compression."
    )]
    pub database_file: PathBuf,

    #[clap(
        short = 'r',
        long = "reads",
        required = true,
        help = "Short-read file (FASTQ). Supports .gz, .xz, .zst, .7z compression."
    )]
    pub reads_file: PathBuf,

    #[clap(
        short,
        long,
        required = true,
        help = "Output file for the IDs of matching reads. Supports .gz, .xz, .zst, .7z compression based on extension."
    )]
    pub output_file: PathBuf,

    #[clap(
        short = 'c',
        long = "min-hits",
        default_value_t = 1,
        help = "Minimum number of k-mer hits to report a read"
    )]
    pub min_hits: usize,
}

#[derive(Parser, Debug)]
pub struct ClassifyArgs {
    #[clap(
        short,
        long,
        required = true,
        help = "Input genome (FASTA) or reads (FASTQ) file. Supports .gz, .xz, .zst, .7z compression."
    )]
    pub input_file: PathBuf,

    #[clap(
        short = 'd',
        long = "databases",
        required = true,
        num_args = 1..,
        help = "One or more k-mer database files (.db). Supports .gz, .xz, .zst, .7z compression."
    )]
    pub database_files: Vec<PathBuf>,

    #[clap(
        short,
        long,
        required = true,
        help = "Output file for classification results (JSON format). Supports .gz, .xz, .zst, .7z compression based on extension."
    )]
    pub output_file: PathBuf,

    #[clap(
        short,
        long,
        help = "Optional: K-mer size to validate against databases. If not provided, uses k from the first database."
    )]
    pub kmer_size: Option<u8>,

    #[clap(
        long,
        default_value_t = 1,
        help = "Minimum frequency for a k-mer in the input to be considered for depth calculation"
    )]
    pub min_kmer_frequency: usize,

    #[clap(
        long,
        default_value_t = 0.0,
        help = "Minimum reference breadth of coverage to include a reference in the output"
    )]
    pub min_coverage: f64,

    #[clap(
        long,
        help = "Optional: Output file path for a TSV summary of the classification results. Supports .gz, .xz, .zst, .7z compression based on extension."
    )]
    pub output_tsv: Option<PathBuf>,
}

pub fn parse_cli() -> Cli {
    Cli::parse()
}
