use assert_cmd::prelude::*;
use orion_kmer::db_types::KmerDbV2; // Updated import
use orion_kmer::kmer::{canonical_u64, seq_to_u64}; // For test data generation
use predicates::prelude::*;
use std::{
    collections::HashSet,
    fs::{self, File},
    io::Write,
    path::PathBuf,
    process::Command,
};
use tempfile::{NamedTempFile, TempDir};

use std::io::Read; // For MultiGzDecoder

// Helper to run build with actual files and load the resulting KmerDbV2 database
fn run_build_with_files_and_load_db(
    k: u8,
    input_file_paths_and_expected_names: Vec<(PathBuf, String)>, // Vec of (path_to_actual_file, expected_ref_name_in_db)
    output_is_compressed: bool, // True if the output DB file should have a .gz extension
) -> Result<KmerDbV2, Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(&project_root);

    let output_dir = TempDir::new()?;
    let mut output_db_path = output_dir.path().join("test_db.db");
    if output_is_compressed {
        output_db_path.set_extension("db.gz");
    }

    cmd.arg("build")
        .arg("-k")
        .arg(k.to_string())
        .arg("-o")
        .arg(&output_db_path);

    for (input_path, _expected_name) in &input_file_paths_and_expected_names {
        cmd.arg("-g").arg(input_path);
    }

    // println!("Running command: {:?}", cmd); // For debugging test setup
    cmd.assert().success();

    let db_bytes = if output_is_compressed {
        let file = File::open(&output_db_path)?;
        let mut decoder = flate2::read::MultiGzDecoder::new(file);
        let mut bytes = Vec::new();
        decoder.read_to_end(&mut bytes)?;
        bytes
    } else {
        fs::read(&output_db_path)?
    };

    let kmer_db_v2: KmerDbV2 = bincode::deserialize(&db_bytes)?;

    // Verify reference names if needed (important if original filenames had extensions)
    // For simplicity in this helper, we assume the caller will verify the KmerDbV2 content including ref names.
    // The `input_file_paths_and_expected_names` provides the expected name for such checks.

    Ok(kmer_db_v2)
}


// Original helper to run build and load the resulting KmerDbV2 database (using content strings)
fn run_build_and_load_db_v2(
    k: u8,
    input_files_content: Vec<(&str, &str)>, // Vec of (filename, content)
) -> Result<KmerDbV2, Box<dyn std::error::Error>> {
    let temp_dir = TempDir::new()?;
    let mut input_file_paths: Vec<PathBuf> = Vec::new();
    // let mut input_filenames: Vec<String> = Vec::new(); // No longer needed here

    for (name, content) in &input_files_content {
        let file_path = temp_dir.path().join(name);
        if let Some(parent_dir) = file_path.parent() {
            fs::create_dir_all(parent_dir)?;
        }
        let mut file = File::create(&file_path)?;
        writeln!(file, "{}", content)?;
        input_file_paths.push(file_path);
        // input_filenames.push(name.to_string()); // No longer needed here
    }

    let string_input_paths: Vec<String> = input_file_paths
        .iter()
        .map(|p| p.to_str().unwrap().to_string())
        .collect();

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);

    let output_db_file = NamedTempFile::new()?;
    let output_db_path_str = output_db_file.path().to_str().unwrap();

    cmd.arg("build")
        .arg("-k")
        .arg(k.to_string())
        .arg("-o")
        .arg(output_db_path_str);

    for input_path_str in &string_input_paths {
        cmd.arg("-g").arg(input_path_str);
    }

    cmd.assert().success();

    let db_bytes = fs::read(output_db_path_str)?;
    let kmer_db_v2: KmerDbV2 = bincode::deserialize(&db_bytes)?;

    Ok(kmer_db_v2)
}

// Helper function to create a HashSet of u64 k-mers from string representations
fn kmers_from_strings(strs: &[&str], k: u8) -> HashSet<u64> {
    strs.iter()
        .map(|s| canonical_u64(seq_to_u64(s.as_bytes(), k).unwrap(), k))
        .collect()
}

const SAMPLE1_FASTA_CONTENT: &str =
    ">seq1\nACGTACGTACGT\n>seq2\nTTTTCCCCGGGGAAAA\n>seq3\nAgCtAgCtNaCcGgTt";
const MINI_FASTA_CONTENT: &str = ">s1\nACGT\n>s2\nACGT"; // Same k-mers

#[test]
fn test_build_simple_fasta_k3() -> Result<(), Box<dyn std::error::Error>> {
    let kmer_db_v2 = run_build_and_load_db_v2(3, vec![("sample1.fasta", SAMPLE1_FASTA_CONTENT)])?;

    assert_eq!(kmer_db_v2.k, 3);
    assert_eq!(kmer_db_v2.references.len(), 1);
    assert!(kmer_db_v2.references.contains_key("sample1.fasta"));

    // Expected unique canonical k-mers for k=3 from sample1.fasta
    // AAA, AAC, ACC, ACG, AGC, CCC, CCG, CTA, GAA, GGA, GTA (11 k-mers)
    let expected_kmers_sample1 = kmers_from_strings(
        &[
            "AAA", "AAC", "ACC", "ACG", "AGC", "CCC", "CCG", "CTA", "GAA", "GGA", "GTA",
        ],
        3,
    );

    assert_eq!(
        kmer_db_v2.references["sample1.fasta"],
        expected_kmers_sample1
    );
    assert_eq!(kmer_db_v2.total_unique_kmers(), 11);

    Ok(())
}

// --- 7z specific tests ---

#[test]
fn test_build_fasta_7z_input_k7() -> Result<(), Box<dyn std::error::Error>> {
    let input_file_path = get_test_data_path("test_input1.fasta.7z");
    if !input_file_path.exists() {
        eprintln!("Skipping 7z build test, input file not found: {:?}", input_file_path);
        return Ok(());
    }
    let expected_ref_name = "test_input1.fasta.7z".to_string();
    let kmer_db = run_build_with_files_and_load_db(7, vec![(input_file_path, expected_ref_name.clone())], false)?;

    assert_eq!(kmer_db.k, 7);
    assert_eq!(kmer_db.references.len(), 1);
    assert_eq!(kmer_db.references[&expected_ref_name], get_expected_kmers_input1_k7());
    assert_eq!(kmer_db.total_unique_kmers(), get_expected_kmers_input1_k7().len());
    Ok(())
}

#[test]
fn test_build_fastq_7z_input_db_output_7z_k6() -> Result<(), Box<dyn std::error::Error>> {
    let input_file_path = get_test_data_path("test_input2.fastq.7z");
     if !input_file_path.exists() {
        eprintln!("Skipping 7z build test (input), input file not found: {:?}", input_file_path);
        return Ok(());
    }
    let expected_ref_name = "test_input2.fastq.7z".to_string();

    // Test 1: 7z input, plain DB output
    let kmer_db_plain_out = run_build_with_files_and_load_db(6, vec![(input_file_path.clone(), expected_ref_name.clone())], false)?;
    assert_eq!(kmer_db_plain_out.k, 6);
    assert_eq!(kmer_db_plain_out.references.len(), 1);
    assert_eq!(kmer_db_plain_out.references[&expected_ref_name], get_expected_kmers_input2_k6());

    // Test 2: 7z input, 7z DB output
    // The helper `run_build_with_files_and_load_db` needs to be aware of .7z output to load it back.
    // It currently only handles .gz or uncompressed.
    // So, similar to count_tests, we will call the command directly and then load the .7z DB using utils::load_kmer_db_v2

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(&project_root);
    let output_dir = TempDir::new()?;
    let output_db_path = output_dir.path().join("test_db.db.7z");

    cmd.arg("build")
        .arg("-k")
        .arg("6")
        .arg("-g")
        .arg(input_file_path)
        .arg("-o")
        .arg(&output_db_path);

    cmd.assert().success();
    assert!(output_db_path.exists());
    assert!(output_db_path.metadata()?.len() > 0);

    // Load the .7z KmerDB using our utility function (which uses get_input_reader)
    let kmer_db_7z_out = orion_kmer::utils::load_kmer_db_v2(&output_db_path)?;
    assert_eq!(kmer_db_7z_out.k, 6);
    assert_eq!(kmer_db_7z_out.references.len(), 1);
    assert_eq!(kmer_db_7z_out.references[&expected_ref_name], get_expected_kmers_input2_k6());

    Ok(())
}

#[test]
fn test_build_duplicate_kmers_k4() -> Result<(), Box<dyn std::error::Error>> {
    let kmer_db_v2 = run_build_and_load_db_v2(4, vec![("mini.fasta", MINI_FASTA_CONTENT)])?;
    assert_eq!(kmer_db_v2.k, 4);
    assert_eq!(kmer_db_v2.references.len(), 1);
    assert!(kmer_db_v2.references.contains_key("mini.fasta"));

    // ACGT (k=4)
    let expected_kmers_mini = kmers_from_strings(&["ACGT"], 4);
    assert_eq!(kmer_db_v2.references["mini.fasta"], expected_kmers_mini);
    assert_eq!(kmer_db_v2.total_unique_kmers(), 1);
    Ok(())
}

#[test]
fn test_build_multiple_files_k4_v2() -> Result<(), Box<dyn std::error::Error>> {
    let kmer_db_v2 = run_build_and_load_db_v2(
        4,
        vec![
            ("s1.fa", ">s1\nACGTACGT"), // K-mers: ACGT, CGTA, GTAC, TACG (CGTA), ACGT
            // Canonical unique for s1.fa: ACGT, CGTA, GTAC
            ("s2.fa", ">s2\nTACGTACG"), // K-mers: TACG (CGTA), ACGT, CGTA, GTAC, TACG (CGTA)
            // Canonical unique for s2.fa: ACGT, CGTA, GTAC
            ("s3.fa", ">s3\nGGGATCCC"), // K-mers: GGGA, GGAT, GATC, ATCC, TCCC
                                        // Canonical: CCC(GGGA), ATCC(GGAT), GATC, ATCC, GGG(TCCC)
                                        // Unique for s3.fa: CCC, ATCC, GATC
        ],
    )?;
    assert_eq!(kmer_db_v2.k, 4);
    assert_eq!(kmer_db_v2.references.len(), 3);
    assert!(kmer_db_v2.references.contains_key("s1.fa"));
    assert!(kmer_db_v2.references.contains_key("s2.fa"));
    assert!(kmer_db_v2.references.contains_key("s3.fa"));

    let expected_s1_kmers = kmers_from_strings(&["ACGT", "CGTA", "GTAC"], 4);
    let expected_s2_kmers = kmers_from_strings(&["ACGT", "CGTA", "GTAC"], 4); // Same as s1
    let expected_s3_kmers = kmers_from_strings(&["GGGA", "GGAT", "GATC", "ATCC", "TCCC"], 4);

    assert_eq!(kmer_db_v2.references["s1.fa"], expected_s1_kmers);
    assert_eq!(kmer_db_v2.references["s2.fa"], expected_s2_kmers);
    assert_eq!(kmer_db_v2.references["s3.fa"], expected_s3_kmers);

    // Total unique k-mers: (ACGT, CGTA, GTAC) U (CCC, ATCC, GATC) = 6
    let all_expected_unique_kmers = expected_s1_kmers
        .union(&expected_s3_kmers)
        .cloned()
        .collect::<HashSet<u64>>();
    assert_eq!(
        kmer_db_v2.total_unique_kmers(),
        all_expected_unique_kmers.len()
    );
    assert_eq!(
        kmer_db_v2.get_all_kmers_unified(),
        all_expected_unique_kmers
    );

    Ok(())
}

#[test]
fn test_build_0_byte_empty_file() -> Result<(), Box<dyn std::error::Error>> {
    // A 0-byte file is unparseable by needletail and should result in a build failure.
    let temp_dir = TempDir::new()?;
    let empty_file_path = temp_dir.path().join("empty.fa");
    // Create a 0-byte file
    File::create(&empty_file_path)?.set_len(0)?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let output_db_file = NamedTempFile::new()?;

    cmd.arg("build")
        .arg("-k")
        .arg("5")
        .arg("-g")
        .arg(empty_file_path)
        .arg("-o")
        .arg(output_db_file.path());

    cmd.assert().failure().stderr(predicate::str::contains(
        "Failed to open or parse FASTA/Q file", // Error from process_sequences_for_file
    ));

    Ok(())
}

#[test]
fn test_build_fasta_with_no_sequences() -> Result<(), Box<dyn std::error::Error>> {
    // A FASTA file with only headers or comments but no actual sequence data
    // should result in a successful build and an empty k-mer set for that reference.
    let kmer_db_v2 = run_build_and_load_db_v2(5, vec![("no_seq.fa", ">header1\n>header2\n")])?;

    assert_eq!(kmer_db_v2.k, 5);
    assert_eq!(kmer_db_v2.references.len(), 1);
    assert!(kmer_db_v2.references.contains_key("no_seq.fa"));
    assert!(kmer_db_v2.references["no_seq.fa"].is_empty());
    assert_eq!(kmer_db_v2.total_unique_kmers(), 0);
    Ok(())
}

#[test]
fn test_build_malformed_fasta_file() -> Result<(), Box<dyn std::error::Error>> {
    // A file that is not valid FASTA (e.g. binary, or bad header) should cause an error during parsing.
    let temp_dir = TempDir::new()?;
    let malformed_file_path = temp_dir.path().join("malformed.fa");
    // Not starting with ">"
    fs::write(&malformed_file_path, "This is not fasta content\nACGT")?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let output_db_file = NamedTempFile::new()?;

    cmd.arg("build")
        .arg("-k")
        .arg("3")
        .arg("-g")
        .arg(malformed_file_path)
        .arg("-o")
        .arg(output_db_file.path());

    // Error message comes from needletail via process_sequences_for_file
    cmd.assert().failure().stderr(predicate::str::contains(
        "Failed to open or parse FASTA/Q file",
    ));
    Ok(())
}

#[test]
fn test_build_invalid_k_too_large() {
    let mut cmd = Command::cargo_bin("orion-kmer").unwrap();
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let dummy_output = NamedTempFile::new().unwrap();
    let temp_dir = TempDir::new().unwrap();
    let dummy_input_path = temp_dir.path().join("dummy.fa");
    fs::write(&dummy_input_path, ">dummy\nACGT").unwrap();

    cmd.arg("build")
        .arg("-k")
        .arg("33")
        .arg("-g")
        .arg(&dummy_input_path)
        .arg("-o")
        .arg(dummy_output.path());
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Invalid K-mer size: 33"));
}

#[test]
fn test_build_file_not_found() {
    let mut cmd = Command::cargo_bin("orion-kmer").unwrap();
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let dummy_output = NamedTempFile::new().unwrap();

    cmd.arg("build")
        .arg("-k")
        .arg("5")
        .arg("-g")
        .arg("nonexistent_file.fasta")
        .arg("-o")
        .arg(dummy_output.path());
    cmd.assert().failure().stderr(predicate::str::contains(
        "Failed to open or parse FASTA/Q file: nonexistent_file.fasta",
    ));
}

// --- Tests for Compressed I/O ---

fn get_test_data_path(file_name: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("tests/data");
    path.push(file_name);
    path
}

// Expected k-mers for test_input1.fasta, k=7
// ACGTACG, CGTACGT, GTACGTA, TACGTAC, GATTACA (twice)
// Canonical: ACGTACG (from ACGTACG & CGTACGT), GTACGTA, TACGTAC, GATTACA
fn get_expected_kmers_input1_k7() -> HashSet<u64> {
    kmers_from_strings(&["ACGTACG", "GTACGTA", "TACGTAC", "GATTACA"], 7)
}

// Expected k-mers for test_input2.fastq, k=6
// CGTACG, GTACGT, TACGTA, GCATGC, CATGCA, ATGCAT, TGCATG, GATTAC
// Canonical: CGTACG (from CGTACG & GTACGT), TACGTA, GCATGC, CATGCA, ATGCAT, TGCATG, GATTAC
fn get_expected_kmers_input2_k6() -> HashSet<u64> {
    kmers_from_strings(
        &[
            "CGTACG", "TACGTA", "GCATGC", "CATGCA", "ATGCAT", "TGCATG", "GATTAC",
        ],
        6,
    )
}


#[test]
fn test_build_fasta_gz_input_k7() -> Result<(), Box<dyn std::error::Error>> {
    let input_file_path = get_test_data_path("test_input1.fasta.gz");
    let expected_ref_name = "test_input1.fasta.gz".to_string();
    let kmer_db = run_build_with_files_and_load_db(7, vec![(input_file_path, expected_ref_name.clone())], false)?;

    assert_eq!(kmer_db.k, 7);
    assert_eq!(kmer_db.references.len(), 1);
    assert_eq!(kmer_db.references[&expected_ref_name], get_expected_kmers_input1_k7());
    assert_eq!(kmer_db.total_unique_kmers(), get_expected_kmers_input1_k7().len());
    Ok(())
}

#[test]
fn test_build_fasta_xz_input_k7() -> Result<(), Box<dyn std::error::Error>> {
    let input_file_path = get_test_data_path("test_input1.fasta.xz");
    let expected_ref_name = "test_input1.fasta.xz".to_string();
    let kmer_db = run_build_with_files_and_load_db(7, vec![(input_file_path, expected_ref_name.clone())], false)?;

    assert_eq!(kmer_db.k, 7);
    assert_eq!(kmer_db.references.len(), 1);
    assert_eq!(kmer_db.references[&expected_ref_name], get_expected_kmers_input1_k7());
    Ok(())
}

#[test]
fn test_build_fasta_zst_input_k7() -> Result<(), Box<dyn std::error::Error>> {
    let input_file_path = get_test_data_path("test_input1.fasta.zst");
    if !input_file_path.exists() {
        eprintln!("Skipping Zstandard build test, input file not found: {:?}", input_file_path);
        return Ok(());
    }
    let expected_ref_name = "test_input1.fasta.zst".to_string();
    let kmer_db = run_build_with_files_and_load_db(7, vec![(input_file_path, expected_ref_name.clone())], false)?;

    assert_eq!(kmer_db.k, 7);
    assert_eq!(kmer_db.references.len(), 1);
    assert_eq!(kmer_db.references[&expected_ref_name], get_expected_kmers_input1_k7());
    Ok(())
}

#[test]
fn test_build_uncompressed_input_gz_output_k7() -> Result<(), Box<dyn std::error::Error>> {
    let input_file_path = get_test_data_path("test_input1.fasta");
    let expected_ref_name = "test_input1.fasta".to_string();
    // Pass true for output_is_compressed
    let kmer_db = run_build_with_files_and_load_db(7, vec![(input_file_path, expected_ref_name.clone())], true)?;

    assert_eq!(kmer_db.k, 7);
    assert_eq!(kmer_db.references.len(), 1);
    assert_eq!(kmer_db.references[&expected_ref_name], get_expected_kmers_input1_k7());
    Ok(())
}

#[test]
fn test_build_fastq_gz_input_gz_output_k6() -> Result<(), Box<dyn std::error::Error>> {
    let input_file_path = get_test_data_path("test_input2.fastq.gz");
    let expected_ref_name = "test_input2.fastq.gz".to_string();
    let kmer_db = run_build_with_files_and_load_db(6, vec![(input_file_path, expected_ref_name.clone())], true)?;

    assert_eq!(kmer_db.k, 6);
    assert_eq!(kmer_db.references.len(), 1);
    assert_eq!(kmer_db.references[&expected_ref_name], get_expected_kmers_input2_k6());
    assert_eq!(kmer_db.total_unique_kmers(), get_expected_kmers_input2_k6().len());
    Ok(())
}

#[test]
fn test_build_multiple_compressed_inputs_k5() -> Result<(), Box<dyn std::error::Error>> {
    let input_file1_path = get_test_data_path("test_input1.fasta.xz");
    let expected_ref_name1 = "test_input1.fasta.xz".to_string();
    let input_file2_path = get_test_data_path("test_input2.fastq.zst");
    let expected_ref_name2 = "test_input2.fastq.zst".to_string();

    if !input_file2_path.exists() {
        eprintln!("Skipping multi-compressed build test, zst input file not found: {:?}", input_file2_path);
        return Ok(());
    }

    let k = 5;
    let kmer_db = run_build_with_files_and_load_db(
        k,
        vec![
            (input_file1_path, expected_ref_name1.clone()),
            (input_file2_path, expected_ref_name2.clone())
        ],
        false
    )?;

    assert_eq!(kmer_db.k, k);
    assert_eq!(kmer_db.references.len(), 2);

    // Expected k-mers for test_input1.fasta, k=5
    // ACGTA, CGTAC, GTACG, TACGT, GATTAC, ATTACA
    // Canonical: ACGTA, CGTAC, GTACG, TACGT, GATTAC, ATTACA
    let expected_kmers1_k5 = kmers_from_strings(&["ACGTA", "CGTAC", "GTACG", "TACGT", "GATTAC", "ATTACA"], k);
    assert_eq!(kmer_db.references[&expected_ref_name1], expected_kmers1_k5);

    // Expected k-mers for test_input2.fastq, k=5
    // CGTAC, GTACG, TACGT, ACGTA, GCATG, CATGC, ATGCA, TGCAT, GATTAC
    // Canonical: CGTAC, GTACG, TACGT, ACGTA, GCATG, CATGC, ATGCA, TGCAT, GATTAC
    let expected_kmers2_k5 = kmers_from_strings(&["CGTAC", "GTACG", "TACGT", "ACGTA", "GCATG", "CATGC", "ATGCA", "TGCAT", "GATTAC"], k);
    assert_eq!(kmer_db.references[&expected_ref_name2], expected_kmers2_k5);

    let all_expected_kmers = expected_kmers1_k5.union(&expected_kmers2_k5).cloned().collect::<HashSet<u64>>();
    assert_eq!(kmer_db.total_unique_kmers(), all_expected_kmers.len());

    Ok(())
}
