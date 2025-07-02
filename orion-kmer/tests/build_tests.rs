use assert_cmd::prelude::*;
use orion_kmer::commands::build::KmerDb;
use predicates::prelude::*;
use std::{
    collections::HashSet,
    fs::{self, File},
    io::Write,
    path::PathBuf,
    process::Command,
};
use tempfile::{NamedTempFile, TempDir}; // Adjust path if KmerDb is moved/made public differently

// Helper to run build and load the resulting database
fn run_build_and_load_db(
    k: u8,
    input_files_content: Vec<(&str, &str)>, // Vec of (filename, content)
) -> Result<KmerDb, Box<dyn std::error::Error>> {
    let temp_dir = TempDir::new()?;
    let mut input_file_paths: Vec<PathBuf> = Vec::new();

    for (name, content) in &input_files_content {
        let file_path = temp_dir.path().join(name);
        fs::create_dir_all(file_path.parent().unwrap())?; // Ensure parent dir exists
        let mut file = File::create(&file_path)?;
        writeln!(file, "{}", content)?;
        input_file_paths.push(file_path);
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
        cmd.arg("-g").arg(input_path_str); // -g for genomes in build
    }

    cmd.assert().success();

    let db_bytes = fs::read(output_db_path_str)?;
    let kmer_db: KmerDb = bincode::deserialize(&db_bytes)?;

    Ok(kmer_db)
}

const SAMPLE1_FASTA_CONTENT: &str =
    ">seq1\nACGTACGTACGT\n>seq2\nTTTTCCCCGGGGAAAA\n>seq3\nAgCtAgCtNaCcGgTt";
const MINI_FASTA_CONTENT: &str = ">s1\nACGT\n>s2\nACGT"; // Same k-mers

#[test]
fn test_build_simple_fasta_k3() -> Result<(), Box<dyn std::error::Error>> {
    let kmer_db = run_build_and_load_db(3, vec![("sample1.fasta", SAMPLE1_FASTA_CONTENT)])?;

    assert_eq!(kmer_db.k, 3);

    // Expected unique canonical k-mers for k=3 from sample1.fasta
    // Based on the corrected counts from count_tests:
    // ACG, CGT, GTA, TAC (from seq1)
    // AAA, GAA, GGA, CCC, CCG (from seq2) (GGG->CCC)
    // AGC, CTA, ACC, AAC (from seq3)
    // Based on the output of the 'count' command which passed its tests,
    // the unique k-mers (canonical string form) are:
    // AAA, AAC, ACC, ACG, AGC, CCC, CCG, CTA, GAA, GGA, GTA
    // There are 11 unique k-mers.
    let mut expected_kmers = HashSet::new();
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"AAA", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"AAC", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"ACC", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"ACG", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"AGC", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"CCC", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"CCG", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"CTA", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"GAA", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"GGA", 3).unwrap(),
        3,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"GTA", 3).unwrap(),
        3,
    ));

    assert_eq!(kmer_db.kmers.len(), 11);
    assert_eq!(kmer_db.kmers, expected_kmers);

    Ok(())
}

#[test]
fn test_build_duplicate_kmers_k4() -> Result<(), Box<dyn std::error::Error>> {
    let kmer_db = run_build_and_load_db(4, vec![("mini.fasta", MINI_FASTA_CONTENT)])?;
    assert_eq!(kmer_db.k, 4);
    // ACGT (k=4) -> ACGT. Canonical: ACGT
    let mut expected_kmers = HashSet::new();
    expected_kmers.insert(orion_kmer::kmer::seq_to_u64(b"ACGT", 4).unwrap());

    assert_eq!(kmer_db.kmers.len(), 1);
    assert_eq!(kmer_db.kmers, expected_kmers);
    Ok(())
}

#[test]
fn test_build_multiple_files_k4() -> Result<(), Box<dyn std::error::Error>> {
    let kmer_db = run_build_and_load_db(
        4,
        vec![
            ("s1.fa", ">s1\nACGTACGT"), // ACGT, CGTA, GTAC, TACG, ACGT
            // Canon: ACGT, CGTA, GTAC, ACGT (TACG -> CGTA)
            // Unique: ACGT, CGTA, GTAC
            ("s2.fa", ">s2\nTACGTACG"), // TACG, ACGT, CGTA, GTAC, TACG
                                        // Canon: CGTA (TACG), ACGT, CGTA, GTAC, CGTA (TACG)
                                        // Unique: ACGT, CGTA, GTAC
        ],
    )?;
    assert_eq!(kmer_db.k, 4);

    let mut expected_kmers = HashSet::new();
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"ACGT", 4).unwrap(),
        4,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"CGTA", 4).unwrap(),
        4,
    ));
    expected_kmers.insert(orion_kmer::kmer::canonical_u64(
        orion_kmer::kmer::seq_to_u64(b"GTAC", 4).unwrap(),
        4,
    ));
    // TACG -> revcomp CGTA. CGTA (01101100=108), TACG (11000110=198). Canonical is CGTA.

    assert_eq!(kmer_db.kmers.len(), 3);
    assert_eq!(kmer_db.kmers, expected_kmers);
    Ok(())
}

#[test]
fn test_build_empty_input_file() -> Result<(), Box<dyn std::error::Error>> {
    // Expect failure because needletail cannot parse an empty file.
    let temp_dir = TempDir::new()?;
    let empty_file_path = temp_dir.path().join("empty.fa");
    fs::write(&empty_file_path, "")?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let output_db_file = NamedTempFile::new()?;

    cmd.arg("build")
        .arg("-k")
        .arg("5")
        .arg("-g")
        .arg(empty_file_path) // -g for genomes
        .arg("-o")
        .arg(output_db_file.path());

    cmd.assert().failure().stderr(predicate::str::contains(
        "Failed to open or parse FASTA file",
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
        "Failed to open or parse FASTA file: nonexistent_file.fasta",
    ));
}
