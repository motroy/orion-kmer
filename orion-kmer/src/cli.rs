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

    #[clap(short, long, required = true, num_args = 1.., help = "One or more input FASTA/FASTQ files (can be gzipped)")]
    pub input_files: Vec<PathBuf>,

    #[clap(
        short,
        long,
        required = true,
        help = "Output file for k-mer counts (kmer<TAB>count)"
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

    #[clap(short = 'g', long = "genomes", required = true, num_args = 1.., help = "One or more input genome assembly files (FASTA)")]
    pub genome_files: Vec<PathBuf>,

    #[clap(
        short,
        long,
        required = true,
        help = "Output path for the binary k-mer database"
    )]
    pub output_file: PathBuf,
}

#[derive(Parser, Debug)]
pub struct CompareArgs {
    #[clap(long, required = true, help = "First k-mer database file")]
    pub db1: PathBuf,

    #[clap(long, required = true, help = "Second k-mer database file")]
    pub db2: PathBuf,

    #[clap(
        short,
        long,
        required = true,
        help = "Output file for comparison stats (JSON format)"
    )]
    pub output_file: PathBuf,
}

#[derive(Parser, Debug)]
pub struct QueryArgs {
    #[clap(
        short = 'd',
        long = "database",
        required = true,
        help = "K-mer database to query against"
    )]
    pub database_file: PathBuf,

    #[clap(
        short = 'r',
        long = "reads",
        required = true,
        help = "Short-read file (FASTQ, can be gzipped)"
    )]
    pub reads_file: PathBuf,

    #[clap(
        short,
        long,
        required = true,
        help = "Output file for the IDs of matching reads"
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
        help = "Input genome (FASTA) or reads (FASTQ) file"
    )]
    pub input_file: PathBuf,

    #[clap(
        short = 'd',
        long = "databases",
        required = true,
        num_args = 1..,
        help = "One or more k-mer database files (.db)"
    )]
    pub database_files: Vec<PathBuf>,

    #[clap(
        short,
        long,
        required = true,
        help = "Output file for classification results (JSON format)"
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
}


pub fn parse_cli() -> Cli {
    Cli::parse()
}
