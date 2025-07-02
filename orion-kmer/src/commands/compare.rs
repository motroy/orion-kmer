use anyhow::{Context, Result};
use log::info;
use serde::Serialize;
use std::fs::File; // Removed BufReader, PathBuf as load_kmer_db_v2 takes &Path

use crate::{
    cli::CompareArgs,
    // KmerDbV2 is loaded by load_kmer_db_v2, direct import not needed here
    // db_types::KmerDbV2,
    errors::OrionKmerError,
    utils::load_kmer_db_v2, // Import the new loading function
};

#[derive(Serialize, Debug)]
struct ComparisonOutput {
    db1_path: String,
    db2_path: String,
    kmer_size: u8,
    db1_total_unique_kmers_across_references: usize, // Name changed for clarity
    db2_total_unique_kmers_across_references: usize, // Name changed for clarity
    intersection_size: usize,
    union_size: usize,
    jaccard_index: f64,
}

// Removed local load_kmer_db function, will use utils::load_kmer_db_v2

pub fn run_compare(args: CompareArgs) -> Result<()> {
    info!("Starting compare command with args: {:?}", args);

    // Load KmerDbV2 instances
    let db1_v2 = load_kmer_db_v2(&args.db1)?;
    let db2_v2 = load_kmer_db_v2(&args.db2)?;

    if db1_v2.k != db2_v2.k {
        return Err(OrionKmerError::KmerSizeMismatch(db1_v2.k, db2_v2.k).into());
    }
    let kmer_size = db1_v2.k;

    // Get the unified set of k-mers for each database
    let db1_all_kmers = db1_v2.get_all_kmers_unified();
    let db2_all_kmers = db2_v2.get_all_kmers_unified();

    let db1_unique_kmers_count = db1_all_kmers.len();
    let db2_unique_kmers_count = db2_all_kmers.len();

    let intersection_size = db1_all_kmers.intersection(&db2_all_kmers).count();

    // |A U B| = |A| + |B| - |A & B|
    let union_size = db1_unique_kmers_count + db2_unique_kmers_count - intersection_size;

    let jaccard_index = if union_size == 0 {
        0.0 // Avoid division by zero if both sets are empty
    } else {
        intersection_size as f64 / union_size as f64
    };

    let output_data = ComparisonOutput {
        db1_path: args.db1.to_string_lossy().into_owned(),
        db2_path: args.db2.to_string_lossy().into_owned(),
        kmer_size,
        db1_total_unique_kmers_across_references: db1_unique_kmers_count,
        db2_total_unique_kmers_across_references: db2_unique_kmers_count,
        intersection_size,
        union_size,
        jaccard_index,
    };

    info!("Comparison results: {:?}", output_data);

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
