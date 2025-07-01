use anyhow::{Context, Result};
use log::info;
use serde::Serialize;
use std::{fs::File, io::BufReader, path::PathBuf};

use crate::{
    cli::CompareArgs,
    commands::build::KmerDb, // Assuming KmerDb is accessible here
    errors::OrionKmerError,
};

#[derive(Serialize, Debug)]
struct ComparisonOutput {
    db1_path: String,
    db2_path: String,
    kmer_size: u8,
    db1_unique_kmers: usize,
    db2_unique_kmers: usize,
    intersection_size: usize,
    union_size: usize,
    jaccard_index: f64,
}

fn load_kmer_db(path: &PathBuf) -> Result<KmerDb> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open k-mer database file: {:?}", path))?;
    let reader = BufReader::new(file);
    let kmer_db: KmerDb = bincode::deserialize_from(reader)
        .with_context(|| format!("Failed to deserialize k-mer database from {:?}", path))?;
    info!(
        "Successfully loaded k-mer database from {:?} (k={}, {} kmers)",
        path,
        kmer_db.k,
        kmer_db.kmers.len()
    );
    Ok(kmer_db)
}

pub fn run_compare(args: CompareArgs) -> Result<()> {
    info!("Starting compare command with args: {:?}", args);

    let db1 = load_kmer_db(&args.db1)?;
    let db2 = load_kmer_db(&args.db2)?;

    if db1.k != db2.k {
        return Err(OrionKmerError::KmerSizeMismatch(db1.k, db2.k).into());
    }
    let kmer_size = db1.k;

    let db1_unique_kmers = db1.kmers.len();
    let db2_unique_kmers = db2.kmers.len();

    let intersection_size = db1.kmers.intersection(&db2.kmers).count();

    // |A U B| = |A| + |B| - |A & B|
    let union_size = db1_unique_kmers + db2_unique_kmers - intersection_size;

    let jaccard_index = if union_size == 0 {
        0.0 // Avoid division by zero if both sets are empty
    } else {
        intersection_size as f64 / union_size as f64
    };

    let output_data = ComparisonOutput {
        db1_path: args.db1.to_string_lossy().into_owned(),
        db2_path: args.db2.to_string_lossy().into_owned(),
        kmer_size,
        db1_unique_kmers,
        db2_unique_kmers,
        intersection_size,
        union_size,
        jaccard_index,
    };

    info!("Comparison results: {:?}", output_data);
    // Outputting results to JSON file
    let output_file = File::create(&args.output_file)
        .with_context(|| format!("Failed to create output JSON file: {:?}", args.output_file))?;

    serde_json::to_writer_pretty(output_file, &output_data)
        .with_context(|| format!("Failed to write comparison JSON to {:?}", args.output_file))?;

    info!(
        "Successfully wrote comparison statistics to {:?}",
        args.output_file
    );

    Ok(())
}
