use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::{
    collections::HashSet,
    fs::{self, File},
    io::Write,
    path::PathBuf,
    process::Command,
};
use tempfile::{NamedTempFile, TempDir};

// Helper function to run 'build' and return the path to the created .db file
// (Copied and adapted from compare_tests.rs - consider a shared test utils module later)
fn run_build_for_query_test(
    k: u8,
    input_files_content: Vec<(&str, &str)>,
    db_output_dir: &TempDir,
    db_name_prefix: &str,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let temp_input_dir = TempDir::new()?;
    let mut input_file_paths: Vec<PathBuf> = Vec::new();

    for (name, content) in &input_files_content {
        let file_path = temp_input_dir.path().join(name);
        fs::create_dir_all(file_path.parent().unwrap())?;
        let mut file = File::create(&file_path)?;
        file.write_all(content.as_bytes())?;
        input_file_paths.push(file_path);
    }

    let string_input_paths: Vec<String> = input_file_paths
        .iter()
        .map(|p| p.to_str().unwrap().to_string())
        .collect();

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);

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

// Helper function to run query and return the set of matching read IDs
fn run_query_and_get_ids(
    db_path: &PathBuf,
    reads_content: &str,
    min_hits: Option<usize>,
) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    let temp_reads_dir = TempDir::new()?;
    let reads_file_path = temp_reads_dir.path().join("query_reads.fastq");
    let mut reads_file = File::create(&reads_file_path)?;
    write!(reads_file, "{}", reads_content)?; // Use write! for direct content

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);

    let output_ids_file = NamedTempFile::new()?;
    let output_ids_path_str = output_ids_file.path().to_str().unwrap();

    cmd.arg("query")
        .arg("-d")
        .arg(db_path)
        .arg("-r")
        .arg(&reads_file_path)
        .arg("-o")
        .arg(output_ids_path_str);

    if let Some(mh) = min_hits {
        cmd.arg("-c").arg(mh.to_string());
    }

    cmd.assert().success();

    let ids_output_str = fs::read_to_string(output_ids_path_str)?;
    let ids_set: HashSet<String> = ids_output_str.lines().map(String::from).collect();
    Ok(ids_set)
}

const DB_FASTA_CONTENT: &str = ">ref_genome_segment\nACGTACGTTTGCATC";
// k=4, Canonicals from ACGTACGTTTGCATC:
// ACGT, CGTA, GTAC, TACG(->CGTA), ACGT(dup), CGTT, GTTT, TTGC, TGCA, GCAT
// Unique: {ACGT, CGTA, GTAC, CGTT, GTTT, TTGC, TGCA, GCAT} (8 k-mers)

const QUERY_FASTQ_CONTENT: &str = "\
@read1_match_many
ACGTACGTTT
+
!!!!!!!!!!
@read2_match_one
TTGCXXXXXX
+
!!!!!!!!!!
@read3_no_match
CCCCCCCCCC
+
!!!!!!!!!!
@read4_match_kmer_short_read
ACG
+
!!!
@read5_match_multiple_hits_but_one_kmer
ACGTACGTACGT
+
!!!!!!!!!!!!
";
// Read1 (ACGTACGTTT) k=4: ACGT, CGTA, GTAC, TACG(->CGTA), ACGT, CGTT, GTTT. Hits: ACGT,CGTA,GTAC,CGTT,GTTT (7 hits)
// Read2 (TTGCXXXXXX) k=4: TTGC. Hits: TTGC (1 hit)
// Read3 (CCCCCCCCCC) k=4: CCCC. Hits: None (0 hits)
// Read4 (ACG) k=4: Too short. (0 hits)
// Read5 (ACGTACGTACGT) k=4: ACGT,CGTA,GTAC,TACG(->CGTA),ACGT,CGTA,GTAC,TACG(->CGTA),ACGT. Hits: ACGT,CGTA,GTAC (9 hits)

#[test]
fn test_query_basic_matches() -> Result<(), Box<dyn std::error::Error>> {
    let k = 4;
    let temp_db_storage_dir = TempDir::new()?; // To hold the .db file for this test
    let db_path = run_build_for_query_test(
        k,
        vec![("db.fa", DB_FASTA_CONTENT)],
        &temp_db_storage_dir,
        "querydb",
    )?;

    let matched_ids = run_query_and_get_ids(&db_path, QUERY_FASTQ_CONTENT, None)?; // Default min_hits = 1

    let mut expected_ids = HashSet::new();
    expected_ids.insert("read1_match_many".to_string()); // Expect ID without '@'
    expected_ids.insert("read2_match_one".to_string()); // Expect ID without '@'
    expected_ids.insert("read5_match_multiple_hits_but_one_kmer".to_string()); // Expect ID without '@'
    // read4 is too short for k=4

    assert_eq!(matched_ids, expected_ids);
    Ok(())
}

#[test]
fn test_query_min_hits_filter() -> Result<(), Box<dyn std::error::Error>> {
    let k = 4;
    let temp_db_storage_dir = TempDir::new()?;
    let db_path = run_build_for_query_test(
        k,
        vec![("db.fa", DB_FASTA_CONTENT)],
        &temp_db_storage_dir,
        "querydb_minhits",
    )?;

    // Read1: 7 hits. Read2: 1 hit. Read5: 9 hits.
    let matched_ids_min2 = run_query_and_get_ids(&db_path, QUERY_FASTQ_CONTENT, Some(2))?;
    let mut expected_ids_min2 = HashSet::new();
    expected_ids_min2.insert("read1_match_many".to_string());
    expected_ids_min2.insert("read5_match_multiple_hits_but_one_kmer".to_string());
    assert_eq!(matched_ids_min2, expected_ids_min2);

    let matched_ids_min8 = run_query_and_get_ids(&db_path, QUERY_FASTQ_CONTENT, Some(8))?;
    let mut expected_ids_min8 = HashSet::new();
    expected_ids_min8.insert("read5_match_multiple_hits_but_one_kmer".to_string()); // Only read5 has >= 8 hits
    assert_eq!(matched_ids_min8, expected_ids_min8);

    let matched_ids_min10 = run_query_and_get_ids(&db_path, QUERY_FASTQ_CONTENT, Some(10))?;
    assert!(matched_ids_min10.is_empty()); // No reads should have 10 hits

    Ok(())
}

#[test]
fn test_query_empty_reads_file() -> Result<(), Box<dyn std::error::Error>> {
    let k = 4;
    let temp_db_storage_dir = TempDir::new()?;
    let db_path = run_build_for_query_test(
        k,
        vec![("db.fa", DB_FASTA_CONTENT)],
        &temp_db_storage_dir,
        "querydb_emptyreads",
    )?;

    // Create an empty temporary FASTQ file
    let temp_reads_dir = TempDir::new()?;
    let reads_file_path = temp_reads_dir.path().join("empty_query_reads.fastq");
    fs::write(&reads_file_path, "")?; // Write empty content

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let output_ids_file = NamedTempFile::new()?;

    cmd.arg("query")
        .arg("-d")
        .arg(db_path)
        .arg("-r")
        .arg(reads_file_path) // Use path to the empty file
        .arg("-o")
        .arg(output_ids_file.path());

    // Expect failure because needletail cannot parse an empty FASTQ file
    cmd.assert().failure().stderr(predicate::str::contains(
        "Failed to open or parse FASTQ file",
    ));
    Ok(())
}

#[test]
fn test_query_db_file_not_found() {
    let mut cmd = Command::cargo_bin("orion-kmer").unwrap();
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);

    let dummy_reads_file = NamedTempFile::new().unwrap();
    let dummy_output_file = NamedTempFile::new().unwrap();

    cmd.arg("query")
        .arg("-d")
        .arg("nonexistent.db")
        .arg("-r")
        .arg(dummy_reads_file.path())
        .arg("-o")
        .arg(dummy_output_file.path());

    cmd.assert().failure().stderr(predicate::str::contains(
        "Failed to open k-mer database file: \"nonexistent.db\"",
    ));
}

#[test]
fn test_query_reads_file_not_found() -> Result<(), Box<dyn std::error::Error>> {
    let k = 4;
    let temp_db_storage_dir = TempDir::new()?;
    let db_path = run_build_for_query_test(
        k,
        vec![("db.fa", DB_FASTA_CONTENT)],
        &temp_db_storage_dir,
        "querydb_noreadsfile",
    )?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let dummy_output_file = NamedTempFile::new().unwrap();

    cmd.arg("query")
        .arg("-d")
        .arg(db_path)
        .arg("-r")
        .arg("nonexistent.fastq")
        .arg("-o")
        .arg(dummy_output_file.path());

    cmd.assert().failure().stderr(predicate::str::contains(
        "Failed to open or parse FASTQ file: \"nonexistent.fastq\"",
    ));
    Ok(())
}
