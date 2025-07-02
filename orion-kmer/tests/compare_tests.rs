use assert_cmd::prelude::*;
use predicates::prelude::*;
use serde_json::Value as JsonValue;
use std::{
    // collections::HashSet, // No longer needed after debug code removal
    fs::{self, File},
    io::Write,
    path::PathBuf,
    process::Command,
};
use tempfile::{NamedTempFile, TempDir};
// use orion_kmer::commands::build::KmerDb; // No longer needed after debug code removal
// use orion_kmer::kmer; // kmer::u64_to_seq was only used in debug prints, now removed.

// Helper function to run 'build' and return the path to the created .db file
fn run_build_for_test(
    k: u8,
    input_files_content: Vec<(&str, &str)>, // Vec of (filename, content)
    db_output_dir: &TempDir,                // Directory to store the .db file
    db_name_prefix: &str,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let temp_input_dir = TempDir::new()?; // Temporary directory for input FASTA files
    let mut input_file_paths: Vec<PathBuf> = Vec::new();

    for (name, content) in &input_files_content {
        let file_path = temp_input_dir.path().join(name);
        fs::create_dir_all(file_path.parent().unwrap())?;
        let mut file = File::create(&file_path)?;
        file.write_all(content.as_bytes())?; // Use write_all to avoid extra newlines
        input_file_paths.push(file_path);
    }

    let string_input_paths: Vec<String> = input_file_paths
        .iter()
        .map(|p| p.to_str().unwrap().to_string())
        .collect();

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);

    // Create a unique name for the db file to avoid collisions if tests run in parallel
    // and NamedTempFile is not suitable for 'keepable' artifacts.
    let db_file_name = format!("{}_{}.db", db_name_prefix, uuid::Uuid::new_v4().simple());
    let output_db_path = db_output_dir.path().join(db_file_name);

    cmd.arg("build")
        .arg("-k")
        .arg(k.to_string())
        .arg("-o")
        .arg(&output_db_path);

    for input_path_str in &string_input_paths {
        cmd.arg("-g").arg(input_path_str);
    }

    cmd.assert().success();

    Ok(output_db_path)
}

const FASTA_DB1: &str = ">seqA\nACGTACGT\n>seqB\nTTTTGGGG"; // k=4: ACGT,CGTA,GTAC,TACG,TTTT,TTTG,TTGG,TGGG (canonicals for these)
// ACGT, CGTA, GTAC, CGTA (TACG->CGTA), TTTT, AAAC(TTTG), CCAA(TTGG), CCCA(TGGG)
// Unique: ACGT, CGTA, GTAC, TTTT, AAAC, CCAA, CCCA (7)
const FASTA_DB2: &str = ">seqC\nACGTACGG\n>seqD\nAAAACCCC"; // k=4: ACGT,CGTA,GTAC,TACG,AAAA,AAAC,AACC,ACCC
// ACGT, CGTA, GTAC, CGTA (TACG->CGTA), AAAA, AAAC, GGTT(AACC), GGGTT(ACCC)
// Unique: ACGT, CGTA, GTAC, AAAA, AAAC, GGTT, GGGTT (7)
// Intersection with DB1: ACGT, CGTA, GTAC, AAAC (4)

#[test]
fn test_compare_basic() -> Result<(), Box<dyn std::error::Error>> {
    let k = 4;
    let temp_db_dir = TempDir::new()?; // To hold the .db files

    let db1_path = run_build_for_test(k, vec![("db1.fa", FASTA_DB1)], &temp_db_dir, "db1")?;
    let db2_path = run_build_for_test(k, vec![("db2.fa", FASTA_DB2)], &temp_db_dir, "db2")?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);

    let output_json_file = NamedTempFile::new()?;
    let output_json_path_str = output_json_file.path().to_str().unwrap();

    cmd.arg("compare")
        .arg("--db1")
        .arg(db1_path.to_str().unwrap())
        .arg("--db2")
        .arg(db2_path.to_str().unwrap())
        .arg("-o")
        .arg(output_json_path_str);

    cmd.assert().success();

    let json_output_str = fs::read_to_string(output_json_path_str)?;
    let json_data: JsonValue = serde_json::from_str(&json_output_str)?;

    assert_eq!(json_data["kmer_size"], k);
    assert_eq!(json_data["db1_total_unique_kmers_across_references"], 8);
    assert_eq!(json_data["db2_total_unique_kmers_across_references"], 9);

    let intersection = 5;
    let union_val = 8 + 9 - intersection;
    assert_eq!(json_data["intersection_size"], intersection);
    assert_eq!(json_data["union_size"], union_val); // Should be 12

    let expected_jaccard = intersection as f64 / union_val as f64; // 5.0 / 12.0
    assert!((json_data["jaccard_index"].as_f64().unwrap() - expected_jaccard).abs() < 1e-6);

    Ok(())
}

#[test]
fn test_compare_identical_databases() -> Result<(), Box<dyn std::error::Error>> {
    let k = 3;
    let temp_db_dir = TempDir::new()?;
    let db_content_str = ">s1\nACGTACGTACGT"; // Canonical k-mers (k=3): {ACG, GTA} (2 unique)

    // Build the database once
    let db_path = run_build_for_test(
        k,
        vec![("identical.fa", db_content_str)],
        &temp_db_dir,
        "db_identical",
    )?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let output_json_file = NamedTempFile::new()?;
    cmd.arg("compare")
        .arg("--db1")
        .arg(&db_path)
        .arg("--db2")
        .arg(&db_path) // Comparing with itself
        .arg("-o")
        .arg(output_json_file.path());
    cmd.assert().success();

    let json_data: JsonValue = serde_json::from_reader(File::open(output_json_file.path())?)?;
    assert_eq!(json_data["kmer_size"], k);
    assert_eq!(json_data["db1_total_unique_kmers_across_references"], 2);
    assert_eq!(json_data["db2_total_unique_kmers_across_references"], 2);
    assert_eq!(json_data["intersection_size"], 2);
    assert_eq!(json_data["union_size"], 2);
    assert!((json_data["jaccard_index"].as_f64().unwrap() - 1.0).abs() < 1e-6);
    Ok(())
}

#[test]
fn test_compare_no_overlap() -> Result<(), Box<dyn std::error::Error>> {
    let k = 5;
    let temp_db_dir = TempDir::new()?;
    let db1_content = ">s1\nAAAAACCCCC"; // Unique k=5: AAAAA, AAAAC, AAACC, AACCC, ACCCC, CCCCC (6)
    let db2_content = ">s2\nTTTTTGGGGG"; // Unique k=5: TTTTT(->AAAAA), TTTTG(->CAAAA), TTTGG(->CCAAA), TTGGG(->CCCAA), TGGGG(->CCCCA), GGGGG(->CCCCC) (6 unique, but AAAAA and CCCCC are common)

    let db1_path = run_build_for_test(
        k,
        vec![("no_overlap1.fa", db1_content)],
        &temp_db_dir,
        "db_nooverlap1",
    )?;
    let db2_path = run_build_for_test(
        k,
        vec![("no_overlap2.fa", db2_content)],
        &temp_db_dir,
        "db_nooverlap2",
    )?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let output_json_file = NamedTempFile::new()?;
    cmd.arg("compare")
        .arg("--db1")
        .arg(&db1_path)
        .arg("--db2")
        .arg(&db2_path)
        .arg("-o")
        .arg(output_json_file.path());
    cmd.assert().success();

    let json_data: JsonValue = serde_json::from_reader(File::open(output_json_file.path())?)?;
    assert_eq!(json_data["kmer_size"], k);
    assert_eq!(json_data["db1_total_unique_kmers_across_references"], 6);
    assert_eq!(json_data["db2_total_unique_kmers_across_references"], 6);

    let expected_intersection = 2; // AAAAA and CCCCC
    let expected_union = 6 + 6 - expected_intersection; // 10
    let expected_jaccard = expected_intersection as f64 / expected_union as f64; // 2.0 / 10.0 = 0.2

    assert_eq!(json_data["intersection_size"], expected_intersection);
    assert_eq!(json_data["union_size"], expected_union);
    assert!((json_data["jaccard_index"].as_f64().unwrap() - expected_jaccard).abs() < 1e-6);
    Ok(())
}

#[test]
fn test_compare_kmer_size_mismatch() -> Result<(), Box<dyn std::error::Error>> {
    let temp_db_dir = TempDir::new()?;
    let db1_k3_path = run_build_for_test(3, vec![("k3.fa", FASTA_DB1)], &temp_db_dir, "db_k3")?;
    let db2_k4_path = run_build_for_test(4, vec![("k4.fa", FASTA_DB2)], &temp_db_dir, "db_k4")?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let output_json_file = NamedTempFile::new()?; // Will not be written
    cmd.arg("compare")
        .arg("--db1")
        .arg(&db1_k3_path)
        .arg("--db2")
        .arg(&db2_k4_path)
        .arg("-o")
        .arg(output_json_file.path());

    cmd.assert().failure().stderr(predicate::str::contains(
        "K-mer databases have incompatible k-mer sizes (overall comparison): 3 vs 4",
    ));
    Ok(())
}
