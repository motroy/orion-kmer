use anyhow::{Context, Result};
use log::{debug, info}; // Removed warn
use serde::Serialize;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::BufWriter,
    // path::PathBuf, // Removed PathBuf, as it's not directly used as a type here
};

use crate::{
    cli::ClassifyArgs,
    db_types::KmerDbV2,
    errors::OrionKmerError,
    kmer::{canonical_u64, seq_to_u64},
    utils::{get_input_reader, get_output_writer, load_kmer_db_v2, track_progress_and_resources}, // Import the wrapper & I/O helpers
};
use csv;
use needletail::{parse_fastx_reader, Sequence}; // Changed to parse_fastx_reader
// use indicatif::ProgressBar; // Not strictly needed for the closure signature if pb is not used inside

// --- Output Structures ---

#[derive(Serialize, Debug)]
struct ReferenceClassificationResult {
    reference_name: String,
    total_kmers_in_reference: usize,
    input_kmers_hitting_reference: usize,
    sum_depth_of_matched_kmers_in_input: usize,
    avg_depth_of_matched_kmers_in_input: f64, // (sum_depth / input_kmers_hitting_reference)
    proportion_input_kmers_hitting_reference: f64, // (input_kmers_hitting_reference / total_unique_input_kmers)
    reference_breadth_of_coverage: f64, // (input_kmers_hitting_reference / total_kmers_in_reference)
}

#[derive(Serialize, Debug)]
struct DatabaseClassificationResult {
    database_path: String,
    database_kmer_size: u8,
    total_unique_kmers_in_db_across_references: usize,
    overall_input_kmers_matched_in_db: usize,
    overall_sum_depth_of_matched_kmers_in_input: usize,
    overall_avg_depth_of_matched_kmers_in_input: f64,
    proportion_input_kmers_in_db_overall: f64, // (overall_input_kmers_matched_in_db / total_unique_input_kmers)
    proportion_db_kmers_covered_overall: f64, // (overall_input_kmers_matched_in_db / total_unique_kmers_in_db_across_references)
    references: Vec<ReferenceClassificationResult>,
}

#[derive(Serialize, Debug)]
struct ClassificationOutput {
    input_file_path: String,
    total_unique_kmers_in_input: usize,
    min_kmer_frequency_filter: usize,
    databases_analyzed: Vec<DatabaseClassificationResult>,
}

// --- Main Logic ---

pub fn run_classify(args: ClassifyArgs) -> Result<()> {
    eprintln!(
        "DEBUG: Entered run_classify. Input file: {:?}, Num DBs: {}, Output: {:?}",
        args.input_file,
        args.database_files.len(),
        args.output_file
    ); // Basic entry print
    info!("Starting classify command with args: {:?}", args);

    // --- 1. Load databases and determine/validate k ---
    let mut loaded_databases: Vec<KmerDbV2> = Vec::new();
    let mut final_k: Option<u8> = None;

    if let Some(user_k) = args.kmer_size {
        if user_k == 0 || user_k > 32 {
            return Err(OrionKmerError::InvalidKmerSize(user_k).into());
        }
        final_k = Some(user_k);
        info!("User specified k-mer size for validation: {}", user_k);
    }

    for (_i, db_path) in args.database_files.iter().enumerate() {
        // Changed i to _i
        let kmer_db = load_kmer_db_v2(db_path)
            .with_context(|| format!("Failed to load database: {:?}", db_path))?;

        if let Some(current_k) = final_k {
            if kmer_db.k != current_k {
                if args.kmer_size.is_some() {
                    // Mismatch against user-provided k
                    return Err(OrionKmerError::KmerSizeMismatchValidation(
                        current_k,
                        kmer_db.k,
                        db_path.clone(),
                    )
                    .into());
                } else {
                    // Mismatch against k from the first database
                    return Err(OrionKmerError::KmerSizeMismatchBetweenDatabases(
                        current_k,
                        kmer_db.k,
                        db_path.clone(),
                    )
                    .into());
                }
            }
        } else {
            // This is the first database and user did not provide k
            if kmer_db.k == 0 || kmer_db.k > 32 {
                // Validate k from first DB
                return Err(OrionKmerError::InvalidKmerSize(kmer_db.k).into());
            }
            final_k = Some(kmer_db.k);
            info!(
                "Effective k-mer size set to {} from first database {:?}",
                kmer_db.k, db_path
            );
        }
        loaded_databases.push(kmer_db);
    }

    let k = match final_k {
        Some(k_val) => {
            if k_val == 0 || k_val > 32 {
                return Err(OrionKmerError::InvalidKmerSize(k_val).into());
            }
            k_val
        }
        None => {
            // Should not happen if args.database_files is required and not empty by clap
            return Err(OrionKmerError::Generic(
                "No databases provided to determine k-mer size.".to_string(),
            )
            .into());
        }
    };
    info!("Processing with effective k-mer size: {}", k);

    // --- 2. Process input file: count k-mers ---
    let mut input_kmer_counts: HashMap<u64, usize> = HashMap::new();
    let input_file_path_str = args.input_file.to_string_lossy().into_owned();

    track_progress_and_resources(
        &format!("Processing input file: {}", input_file_path_str),
        0, // 0 for indeterminate progress bar (spinner style) as we don't know total records easily
        |pb_input| {
            // Use get_input_reader to handle potential compression
            let input_buf_reader = get_input_reader(&args.input_file).with_context(|| {
                format!(
                    "Failed to get input reader for file: {:?}",
                    args.input_file
                )
            })?;
            // Pass the BufRead to parse_fastx_reader
            let mut reader = parse_fastx_reader(input_buf_reader).with_context(|| {
                format!(
                    "Failed to parse FASTA/Q content from: {:?}",
                    args.input_file
                )
            })?;

            let mut processed_records = 0;
            while let Some(record) = reader.next() {
                let record = record.with_context(|| {
                    format!(
                        "Error reading record from input file: {:?}",
                        args.input_file
                    )
                })?;
                let norm_seq = record.normalize(false);

                if norm_seq.len() >= k as usize {
                    for window in norm_seq.windows(k as usize) {
                        if let Some(kmer_val) = seq_to_u64(window, k) {
                            let canonical_kmer = canonical_u64(kmer_val, k);
                            *input_kmer_counts.entry(canonical_kmer).or_insert(0) += 1;
                        }
                    }
                }
                processed_records += 1;
                if processed_records % 100_000 == 0 {
                    // Update progress bar message periodically
                    pb_input.set_message(format!("Processed {} records...", processed_records));
                }
                // We don't call pb_input.inc() here as length is 0 (spinner)
            }
            pb_input.set_message(format!(
                "Processed {} total records from input file.",
                processed_records
            ));
            Ok(())
        },
    )?;

    info!(
        "Finished processing input file. Found {} unique k-mers with total occurrences before frequency filtering.",
        input_kmer_counts.len()
    );

    // Filter input_kmer_counts by min_kmer_frequency
    let filtered_input_kmer_counts: HashMap<u64, usize> = input_kmer_counts
        .into_iter()
        .filter(|&(_, count)| count >= args.min_kmer_frequency)
        .collect();

    let total_unique_input_kmers_after_filter = filtered_input_kmer_counts.len();
    info!(
        "After applying min_kmer_frequency filter (>= {}), {} unique k-mers remain in input.",
        args.min_kmer_frequency, total_unique_input_kmers_after_filter
    );

    // --- 3. Perform classification ---
    let mut db_results: Vec<DatabaseClassificationResult> = Vec::new();
    let num_databases = loaded_databases.len() as u64;

    track_progress_and_resources(
        "Classifying against databases",
        num_databases,
        |pb_classify| {
            for (idx, kmer_db_v2) in loaded_databases.iter().enumerate() {
                let db_path_str = args.database_files[idx].to_string_lossy().into_owned();
                info!("Classifying against database: {}", db_path_str);
                pb_classify.set_message(format!("Classifying against: {}", db_path_str));

                let mut overall_matched_kmers_in_db_set: HashSet<u64> = HashSet::new();
                let overall_sum_depth_for_db: usize;
                let mut per_reference_results: Vec<ReferenceClassificationResult> = Vec::new();

                for (ref_name, ref_kmers_set) in &kmer_db_v2.references {
                    debug!("Processing reference: {} from {}", ref_name, db_path_str);
                    // ---- END DEBUG PRINT ----
                    let mut matched_kmers_for_ref_set: HashSet<u64> = HashSet::new();
                    let mut sum_depth_for_ref: usize = 0;

                    for (input_kmer, input_count) in &filtered_input_kmer_counts {
                        if ref_kmers_set.contains(input_kmer) {
                            matched_kmers_for_ref_set.insert(*input_kmer);
                            sum_depth_for_ref += input_count;
                            overall_matched_kmers_in_db_set.insert(*input_kmer); // Add to overall set for the DB
                        }
                    }

                    let num_matched_for_ref = matched_kmers_for_ref_set.len();
                    let total_kmers_in_ref = ref_kmers_set.len();

                    let reference_breadth_of_coverage = if total_kmers_in_ref > 0 {
                        num_matched_for_ref as f64 / total_kmers_in_ref as f64
                    } else {
                        0.0
                    };

                    if reference_breadth_of_coverage >= args.min_coverage {
                        per_reference_results.push(ReferenceClassificationResult {
                            reference_name: ref_name.clone(),
                            total_kmers_in_reference: total_kmers_in_ref,
                            input_kmers_hitting_reference: num_matched_for_ref,
                            sum_depth_of_matched_kmers_in_input: sum_depth_for_ref,
                            avg_depth_of_matched_kmers_in_input: if num_matched_for_ref > 0 {
                                sum_depth_for_ref as f64 / num_matched_for_ref as f64
                            } else {
                                0.0
                            },
                            proportion_input_kmers_hitting_reference:
                                if total_unique_input_kmers_after_filter > 0 {
                                    num_matched_for_ref as f64
                                        / total_unique_input_kmers_after_filter as f64
                                } else {
                                    0.0
                                },
                            reference_breadth_of_coverage,
                        });
                    }
                }

                // Calculate sum of depths for overall_matched_kmers_in_db_set
                // This needs to iterate input_kmer_counts again, specifically for k-mers in overall_matched_kmers_in_db_set
                overall_sum_depth_for_db = overall_matched_kmers_in_db_set
                    .iter()
                    .map(|kmer| filtered_input_kmer_counts.get(kmer).copied().unwrap_or(0))
                    .sum();

                let num_overall_matched_kmers = overall_matched_kmers_in_db_set.len();
                // ---- END DEBUG PRINT ----
                let total_kmers_in_db_union = kmer_db_v2.total_unique_kmers();

                db_results.push(DatabaseClassificationResult {
                    database_path: db_path_str,
                    database_kmer_size: kmer_db_v2.k,
                    total_unique_kmers_in_db_across_references: total_kmers_in_db_union,
                    overall_input_kmers_matched_in_db: num_overall_matched_kmers,
                    overall_sum_depth_of_matched_kmers_in_input: overall_sum_depth_for_db,
                    overall_avg_depth_of_matched_kmers_in_input: if num_overall_matched_kmers > 0 {
                        overall_sum_depth_for_db as f64 / num_overall_matched_kmers as f64
                    } else {
                        0.0
                    },
                    proportion_input_kmers_in_db_overall: if total_unique_input_kmers_after_filter
                        > 0
                    {
                        num_overall_matched_kmers as f64
                            / total_unique_input_kmers_after_filter as f64
                    } else {
                        0.0
                    },
                    proportion_db_kmers_covered_overall: if total_kmers_in_db_union > 0 {
                        num_overall_matched_kmers as f64 / total_kmers_in_db_union as f64
                    } else {
                        0.0
                    },
                    references: per_reference_results,
                });
                pb_classify.inc(1); // Increment after processing each database
            }
            Ok(()) // Return Ok from the closure
        },
    )?;

    // --- 4. Write output ---
    let final_output = ClassificationOutput {
        input_file_path: args.input_file.to_string_lossy().into_owned(),
        total_unique_kmers_in_input: total_unique_input_kmers_after_filter,
        min_kmer_frequency_filter: args.min_kmer_frequency,
        databases_analyzed: db_results,
    };

    info!("Writing classification results to: {:?}", args.output_file);
    // Use get_output_writer for the main JSON output
    let mut writer = get_output_writer(&args.output_file).with_context(|| {
        format!(
            "Failed to get output writer for JSON file: {:?}",
            args.output_file
        )
    })?;
    serde_json::to_writer_pretty(&mut writer, &final_output).with_context(|| {
        format!(
            "Failed to write classification JSON to {:?}",
            args.output_file
        )
    })?;
    writer.flush().context("Failed to flush JSON output writer")?;


    // --- 5. Optionally write TSV output ---
    if let Some(tsv_path) = &args.output_tsv {
        info!("Writing classification summary TSV to: {:?}", tsv_path);
        // Use get_output_writer for the TSV output
        let tsv_writer_boxed = get_output_writer(tsv_path).with_context(|| {
            format!("Failed to get output writer for TSV file: {:?}", tsv_path)
        })?;
        let mut tsv_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(tsv_writer_boxed); // from_writer expects W: Write

        // Write header
        tsv_writer.write_record(&[
            "InputFile",
            "Database",
            "Reference",
            "TotalKmersInReference",
            "InputKmersHittingReference",
            "SumDepthMatchedKmers",
            "AvgDepthMatchedKmers",
            "ProportionInputKmersHittingReference",
            "ReferenceBreadthOfCoverage",
        ])?;

        // Write data rows
        for db_res in &final_output.databases_analyzed {
            for ref_res in &db_res.references {
                // Note: ref_res here are already filtered by min_coverage
                tsv_writer.write_record(&[
                    final_output.input_file_path.clone(),
                    db_res.database_path.clone(),
                    ref_res.reference_name.clone(),
                    ref_res.total_kmers_in_reference.to_string(),
                    ref_res.input_kmers_hitting_reference.to_string(),
                    ref_res.sum_depth_of_matched_kmers_in_input.to_string(),
                    format!("{:.4}", ref_res.avg_depth_of_matched_kmers_in_input),
                    format!("{:.4}", ref_res.proportion_input_kmers_hitting_reference),
                    format!("{:.4}", ref_res.reference_breadth_of_coverage),
                ])?;
            }
        }
        tsv_writer.flush()?;
        info!("TSV summary successfully written to {:?}", tsv_path);
    }

    info!("Classification successfully completed.");
    Ok(())
}
