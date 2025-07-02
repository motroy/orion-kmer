use assert_cmd::prelude::*;
use csv;
use predicates::prelude::*;
use serde_json::Value as JsonValue;
use std::{
    fs::{self, File},
    io::Write,
    path::PathBuf, // Path is not directly used
    process::Command,
};
use tempfile::{NamedTempFile, TempDir};

// Helper to build a database for use in classify tests
// Similar to run_build_for_test in compare_tests.rs
fn build_db_for_classify(
    k: u8,
    input_files_content: Vec<(&str, &str)>, // Vec of (filename_in_db, content)
    db_output_dir: &TempDir,
    db_name_prefix: &str,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let temp_input_dir = TempDir::new()?;
    let mut input_file_paths_for_build: Vec<PathBuf> = Vec::new();

    for (name, content) in &input_files_content {
        let file_path = temp_input_dir.path().join(name);
        if let Some(parent_dir) = file_path.parent() {
            fs::create_dir_all(parent_dir)?;
        }
        let mut file = File::create(&file_path)?;
        file.write_all(content.as_bytes())?;
        input_file_paths_for_build.push(file_path);
    }

    let mut cmd_build = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd_build.current_dir(&project_root); // Use a reference for current_dir

    let db_file_name = format!("{}_{}.db", db_name_prefix, uuid::Uuid::new_v4().simple());
    let output_db_path = db_output_dir.path().join(db_file_name);

    cmd_build
        .arg("build")
        .arg("-k")
        .arg(k.to_string())
        .arg("-o")
        .arg(&output_db_path);

    for input_path_str_obj in &input_file_paths_for_build {
        cmd_build
            .arg("-g")
            .arg(input_path_str_obj.to_str().unwrap());
    }

    cmd_build.assert().success();
    Ok(output_db_path)
}

// Helper to run the classify command and return the parsed JSON output
fn run_classify_get_json(
    input_content: &str,  // Content of the main input file for classification
    input_filename: &str, // Filename for the main input file
    db_paths: &[PathBuf],
    k_value_user: Option<u8>, // Optional k value provided by user
    min_kmer_freq: Option<usize>,
    min_coverage: Option<f64>,
    output_tsv_path_option: Option<PathBuf>,
) -> Result<JsonValue, Box<dyn std::error::Error>> {
    let temp_dir = TempDir::new()?;

    // Create the main input file
    let main_input_file_path = temp_dir.path().join(input_filename);
    let mut main_input_file = File::create(&main_input_file_path)?;
    writeln!(main_input_file, "{}", input_content)?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);

    let output_json_file = NamedTempFile::new()?;
    let output_json_path_str = output_json_file.path().to_str().unwrap();

    cmd.arg("classify")
        .arg("-i")
        .arg(main_input_file_path.to_str().unwrap())
        .arg("-o")
        .arg(output_json_path_str);

    for db_path in db_paths {
        cmd.arg("-d").arg(db_path.to_str().unwrap());
    }

    if let Some(k) = k_value_user {
        cmd.arg("--kmer-size").arg(k.to_string());
    }
    if let Some(freq) = min_kmer_freq {
        cmd.arg("--min-kmer-frequency").arg(freq.to_string());
    }
    if let Some(cov) = min_coverage {
        cmd.arg("--min-coverage").arg(cov.to_string());
    }
    if let Some(tsv_path) = &output_tsv_path_option {
        cmd.arg("--output-tsv").arg(tsv_path);
    }

    cmd.assert().success();

    let json_output_str = fs::read_to_string(output_json_path_str)?;
    let json_data: JsonValue = serde_json::from_str(&json_output_str)?;
    Ok(json_data)
}

// --- Test Cases ---

const INPUT_FASTA_BASIC: &str =
    ">input_seq1\nACGTACGT\n>input_seq2\nACGTACGT\n>input_seq3\nTTTTGGGG";
// For k=4:
// input_seq1: ACGT (c), CGTA (c), GTAC (c), TACG (c->CGTA) => ACGT, CGTA, GTAC. Counts: ACGT:1, CGTA:2, GTAC:1
// input_seq2: ACGT (c), CGTA (c), GTAC (c), TACG (c->CGTA) => ACGT, CGTA, GTAC. Counts: ACGT:1, CGTA:2, GTAC:1
// input_seq3: TTTT (c), TTTG (c->AAAC), TTGG (c->CCAA), TGGG (c->CCCA) => TTTT, AAAC, CCAA, CCCA
// Total unique input k-mers: {ACGT, CGTA, GTAC, TTTT, AAAC, CCAA, CCCA} (7 unique)
// Counts in input: ACGT:2, CGTA:4, GTAC:2, TTTT:1, AAAC:1, CCAA:1, CCCA:1

// DB1 will contain two references
const DB1_REF1_FASTA: &str = ">db1_refA\nACGTACGTACGT"; // k=4: ACGT, CGTA, GTAC
const DB1_REF2_FASTA: &str = ">db1_refB\nGGGAAAAATTTT"; // k=4: GGGA(CCCA), GGAA(TTCC), GAAA(TTTC), AAAA, AAAT(ATTT), AATT, ATTT, TTTT

// DB2 will contain one reference
const DB2_REF1_FASTA: &str = ">db2_refC\nACGTTACGTT"; // k=4: ACGT, CGTT, GTTT(AAAC), TTAC(GTAA), TACG(CGTA)
// Unique: ACGT, CGTT, AAAC, GTAA, CGTA

#[test]
fn test_classify_basic_fasta_input() -> Result<(), Box<dyn std::error::Error>> {
    let k = 4;
    let temp_db_storage = TempDir::new()?;

    // Build DB1
    let db1_path = build_db_for_classify(
        k,
        vec![
            ("db1_refA.fa", DB1_REF1_FASTA),
            ("db1_refB.fa", DB1_REF2_FASTA),
        ],
        &temp_db_storage,
        "db1",
    )?;
    // Build DB2
    let db2_path = build_db_for_classify(
        k,
        vec![("db2_refC.fa", DB2_REF1_FASTA)],
        &temp_db_storage,
        "db2",
    )?;

    let results = run_classify_get_json(
        INPUT_FASTA_BASIC,
        "input.fa",
        &[db1_path.clone(), db2_path.clone()], // Use clone if paths are needed later
        Some(k),
        None, // Default min_kmer_frequency = 1
        None, // Default min_coverage = 0.0
        None, // No TSV output for this test
    )?;

    assert!(
        results["input_file_path"]
            .as_str()
            .unwrap()
            .ends_with("input.fa")
    );
    assert_eq!(results["total_unique_kmers_in_input"], 8); // Corrected from 7 to 8
    assert_eq!(results["min_kmer_frequency_filter"], 1);
    assert_eq!(results["databases_analyzed"].as_array().unwrap().len(), 2);

    // --- Assertions for DB1 ---
    let db1_results = &results["databases_analyzed"][0];
    assert_eq!(db1_results["database_path"], db1_path.to_str().unwrap());
    assert_eq!(db1_results["database_kmer_size"], k);
    // DB1_REF1_FASTA (ACGT, CGTA, GTAC) = 3
    // DB1_REF2_FASTA (CCCA, TTCC, TTTC, AAAA, ATTT, AATT) = 6. (TTTT from content becomes AAAA)
    // Union size = 3 + 6 = 9
    assert_eq!(db1_results["total_unique_kmers_in_db_across_references"], 9);

    // Input unique k-mers (k=4, min_freq=1): {ACGT:4, CGTA:4, GTAC:2, AAAA:1, CAAA:1, CCAA:1, CCCA:1, CCCC:1} (8 unique)
    // DB1 references k-mers:
    //  refA: {ACGT, CGTA, GTAC}
    //  refB: {GGGA(canon GGGA), AAAA(canon AAAA), TTCC, TTTC, ATTT, AATT} (6 unique)
    // Overall matches for DB1 (Input vs DB1): {ACGT, CGTA, GTAC, AAAA} (4 unique)
    // Input CCCA (canon CCCA) does NOT match DB GGGA (canon GGGA).
    assert_eq!(db1_results["overall_input_kmers_matched_in_db"], 4); // Corrected from 5
    // Depths from input (corrected): ACGT(4), CGTA(4), GTAC(2), AAAA(1) -> Sum = 4+4+2+1 = 11
    assert_eq!(
        db1_results["overall_sum_depth_of_matched_kmers_in_input"],
        11
    ); // Corrected from 12
    assert!(
        (db1_results["overall_avg_depth_of_matched_kmers_in_input"]
            .as_f64()
            .unwrap()
            - (11.0 / 4.0))
            .abs()
            < 1e-6
    ); // Corrected
    assert!(
        (db1_results["proportion_input_kmers_in_db_overall"]
            .as_f64()
            .unwrap()
            - (4.0 / 8.0))
            .abs()
            < 1e-6
    ); // Corrected
    assert!(
        (db1_results["proportion_db_kmers_covered_overall"]
            .as_f64()
            .unwrap()
            - (4.0 / 9.0))
            .abs()
            < 1e-6
    ); // Corrected

    assert_eq!(db1_results["references"].as_array().unwrap().len(), 2);

    // DB1 Ref A ("db1_refA.fa") - k-mers {ACGT, CGTA, GTAC}
    // Matches from input: {ACGT, CGTA, GTAC} (3 unique)
    let db1_refa_res = db1_results["references"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["reference_name"] == "db1_refA.fa")
        .unwrap();
    assert_eq!(db1_refa_res["total_kmers_in_reference"], 3);
    assert_eq!(db1_refa_res["input_kmers_hitting_reference"], 3);
    // Depths for these from input: ACGT(4) + CGTA(4) + GTAC(2) = 10
    assert_eq!(db1_refa_res["sum_depth_of_matched_kmers_in_input"], 10);
    assert!(
        (db1_refa_res["avg_depth_of_matched_kmers_in_input"]
            .as_f64()
            .unwrap()
            - (10.0 / 3.0))
            .abs()
            < 1e-6
    );
    assert!(
        (db1_refa_res["proportion_input_kmers_hitting_reference"]
            .as_f64()
            .unwrap()
            - (3.0 / 8.0))
            .abs()
            < 1e-6
    );
    assert!(
        (db1_refa_res["reference_breadth_of_coverage"]
            .as_f64()
            .unwrap()
            - (3.0 / 3.0))
            .abs()
            < 1e-6
    );

    // DB1 Ref B ("db1_refB.fa") - k-mers {GGGA(canon GGGA), AAAA(canon AAAA), TTCC, TTTC, ATTT, AATT}
    // Matches from input: {AAAA} (1 unique)
    // Input CCCA (canon CCCA) does not match DB GGGA (canon GGGA).
    let db1_refb_res = db1_results["references"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["reference_name"] == "db1_refB.fa")
        .unwrap();
    assert_eq!(db1_refb_res["total_kmers_in_reference"], 6);
    assert_eq!(db1_refb_res["input_kmers_hitting_reference"], 1); // Corrected from 2
    // Depths for these from input: AAAA(1) = 1
    assert_eq!(db1_refb_res["sum_depth_of_matched_kmers_in_input"], 1); // Corrected from 2
    assert!(
        (db1_refb_res["avg_depth_of_matched_kmers_in_input"]
            .as_f64()
            .unwrap()
            - (1.0 / 1.0))
            .abs()
            < 1e-6
    ); // Corrected
    assert!(
        (db1_refb_res["proportion_input_kmers_hitting_reference"]
            .as_f64()
            .unwrap()
            - (1.0 / 8.0))
            .abs()
            < 1e-6
    ); // Corrected
    assert!(
        (db1_refb_res["reference_breadth_of_coverage"]
            .as_f64()
            .unwrap()
            - (1.0 / 6.0))
            .abs()
            < 1e-6
    ); // Corrected

    // --- Assertions for DB2 ---
    // DB2_REF1_FASTA (`ACGTTACGTT` -> k=4: {ACGT, CGTT, AAAC (from GTTT), GTAA (from TTAC), CGTA (from TACG)}) = 5 unique k-mers
    // Input unique k-mers (k=4, min_freq=1): {ACGT:4, CGTA:4, GTAC:2, AAAA:1, CAAA:1, CCAA:1, CCCA:1, CCCC:1} (8 unique)
    // Input k-mer CAAA is AAAC when canonicalized if needed by string comparison, but as u64, CAAA != AAAC typically.
    // Let's use actual canon u64 values.
    // Input k-mer from TTTG is CAAA (canon). DB2 k-mer from GTTT is AAAC (canon). These are different.
    // My comment "Input has AAAC (from TTTG)" was wrong. TTTG -> CAAA.
    // DB2_REF1_FASTA: ACGT, CGTT, GTTT(AAAC), TTAC(GTAA), TACG(CGTA). Kmer set: {ACGT, CGTT, AAAC, GTAA, CGTA}
    // Input k-mers: ACGT, CGTA, GTAC, AAAA, CAAA, CCAA, CCCA, CCCC
    // Matches for DB2: ACGT, CGTA. (2 unique)
    let db2_results = &results["databases_analyzed"][1];
    assert_eq!(db2_results["database_path"], db2_path.to_str().unwrap());
    assert_eq!(db2_results["total_unique_kmers_in_db_across_references"], 5);
    assert_eq!(db2_results["overall_input_kmers_matched_in_db"], 2); // Corrected from 3
    // Depths: ACGT(4), CGTA(4) -> Sum = 4+4 = 8
    assert_eq!(
        db2_results["overall_sum_depth_of_matched_kmers_in_input"],
        8
    ); // Corrected from 7
    assert!(
        (db2_results["overall_avg_depth_of_matched_kmers_in_input"]
            .as_f64()
            .unwrap()
            - (8.0 / 2.0))
            .abs()
            < 1e-6
    ); // Corrected
    assert!(
        (db2_results["proportion_input_kmers_in_db_overall"]
            .as_f64()
            .unwrap()
            - (2.0 / 8.0))
            .abs()
            < 1e-6
    ); // Corrected
    assert!(
        (db2_results["proportion_db_kmers_covered_overall"]
            .as_f64()
            .unwrap()
            - (2.0 / 5.0))
            .abs()
            < 1e-6
    ); // Corrected

    assert_eq!(db2_results["references"].as_array().unwrap().len(), 1);
    let db2_refc_res = &db2_results["references"][0];
    assert_eq!(db2_refc_res["reference_name"], "db2_refC.fa");
    assert_eq!(db2_refc_res["total_kmers_in_reference"], 5);
    assert_eq!(db2_refc_res["input_kmers_hitting_reference"], 2); // Corrected from 3
    assert_eq!(db2_refc_res["sum_depth_of_matched_kmers_in_input"], 8); // Corrected from 7
    assert!(
        (db2_refc_res["avg_depth_of_matched_kmers_in_input"]
            .as_f64()
            .unwrap()
            - (8.0 / 2.0))
            .abs()
            < 1e-6
    ); // Corrected
    assert!(
        (db2_refc_res["proportion_input_kmers_hitting_reference"]
            .as_f64()
            .unwrap()
            - (2.0 / 8.0))
            .abs()
            < 1e-6
    ); // Corrected
    assert!(
        (db2_refc_res["reference_breadth_of_coverage"]
            .as_f64()
            .unwrap()
            - (2.0 / 5.0))
            .abs()
            < 1e-6
    ); // Corrected

    Ok(())
}

#[test]
fn test_classify_min_kmer_frequency_filter() -> Result<(), Box<dyn std::error::Error>> {
    let k = 4;
    let temp_db_storage = TempDir::new()?;

    // DB contains: ACGT, CGTA, GTAC
    let db_path = build_db_for_classify(
        k,
        vec![("db_ref.fa", DB1_REF1_FASTA)],
        &temp_db_storage,
        "db_minfreq",
    )?;

    // New simpler input:
    // >S1\nACGTACGT -> ACGT (1), CGTA (2), GTAC (1)
    // >S2\nACGTGGGG -> ACGT (1), GGGG(CCCC), GGGT(ACCA), GGTG(CACC), GTGG(CCAC)
    // Combined counts: ACGT:2, CGTA:1, GTAC:1, CCCC:1, ACCA:1, CACC:1, CCAC:1
    let simpler_input_content = ">S1\nACGTACGT\n>S2\nACGTGGGG";

    // With min_kmer_frequency = 2:
    // Filtered input k-mers should be: {ACGT: 2}. (1 unique k-mer)
    // This k-mer (ACGT) is in the DB.
    // Expected sum_depth_for_ref for "db_ref.fa" should be 2.

    let results = run_classify_get_json(
        simpler_input_content,
        "input_simple_minfreq.fa",
        &[db_path.clone()],
        Some(k),
        Some(2), // Set min_kmer_frequency to 2
        None,    // Default min_coverage
        None,    // No TSV output
    )?;

    assert_eq!(
        results["total_unique_kmers_in_input"], 2,
        "Total unique input k-mers after filter: {{ACGT:3, CGTA:2}}"
    );
    assert_eq!(results["min_kmer_frequency_filter"], 2);

    let db_res = &results["databases_analyzed"][0];
    assert_eq!(
        db_res["total_unique_kmers_in_db_across_references"], 3,
        "K-mers in DB: {{ACGT, CGTA, GTAC}}"
    );
    assert_eq!(
        db_res["overall_input_kmers_matched_in_db"], 2,
        "Overall matched k-mers: {{ACGT, CGTA}}"
    );

    let actual_sum_depth = db_res["overall_sum_depth_of_matched_kmers_in_input"]
        .as_u64()
        .unwrap_or(0);
    let expected_sum_depth = (3 + 2) as u64; // ACGT count 3, CGTA count 2
    if actual_sum_depth != expected_sum_depth {
        eprintln!("Debug info for test_classify_min_kmer_frequency_filter (simplified):");
        eprintln!(
            "Actual overall_sum_depth_of_matched_kmers_in_input: {}",
            actual_sum_depth
        );
        eprintln!(
            "Expected overall_sum_depth_of_matched_kmers_in_input: {}",
            expected_sum_depth
        );
        eprintln!(
            "Full results JSON: {}",
            serde_json::to_string_pretty(&results)
                .unwrap_or_else(|_| "Failed to serialize results".to_string())
        );
    }
    assert_eq!(
        actual_sum_depth, expected_sum_depth,
        "Mismatch in overall_sum_depth_of_matched_kmers_in_input"
    );

    let ref_res = &db_res["references"][0]; // Assuming "db_ref.fa" is the only one
    assert_eq!(
        ref_res["input_kmers_hitting_reference"], 2,
        "Input k-mers hitting reference: {{ACGT, CGTA}}"
    );
    assert_eq!(
        ref_res["sum_depth_of_matched_kmers_in_input"]
            .as_u64()
            .unwrap(),
        expected_sum_depth,
        "Per-reference sum of depths"
    );
    assert!(
        (ref_res["proportion_input_kmers_hitting_reference"]
            .as_f64()
            .unwrap()
            - (2.0 / 2.0))
            .abs()
            < 1e-6,
        "Proportion input hitting ref"
    ); // 2 matched / 2 unique in input
    assert!(
        (ref_res["reference_breadth_of_coverage"].as_f64().unwrap() - (2.0 / 3.0)).abs() < 1e-6,
        "Reference breadth"
    ); // 2 matched / 3 in ref

    Ok(())
}

#[test]
fn test_classify_k_validation_error() -> Result<(), Box<dyn std::error::Error>> {
    let temp_db_storage = TempDir::new()?;
    let db_k4_path = build_db_for_classify(
        4,
        vec![("dbk4.fa", DB1_REF1_FASTA)],
        &temp_db_storage,
        "dbk4_kvalid",
    )?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let output_json_file = NamedTempFile::new()?; // Will not be written to

    cmd.arg("classify")
        .arg("-i")
        .arg("dummy_input.fa") // Content doesn't matter as it should fail before parsing it
        .arg("-d")
        .arg(db_k4_path.to_str().unwrap())
        .arg("--kmer-size") // User provides k=3
        .arg("3")
        .arg("-o")
        .arg(output_json_file.path());

    cmd.assert().failure().stderr(predicate::str::contains(
        "User-provided k-mer size 3 does not match k-mer size 4 from database",
    ));
    Ok(())
}

#[test]
fn test_classify_k_mismatch_between_databases() -> Result<(), Box<dyn std::error::Error>> {
    let temp_db_storage = TempDir::new()?;
    let db_k4_path = build_db_for_classify(
        4,
        vec![("dbk4.fa", DB1_REF1_FASTA)],
        &temp_db_storage,
        "dbk4_mismatch",
    )?;
    let db_k3_path = build_db_for_classify(
        3,
        vec![("dbk3.fa", ">seq\nACG")],
        &temp_db_storage,
        "dbk3_mismatch",
    )?;

    let mut cmd = Command::cargo_bin("orion-kmer")?;
    let project_root = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    cmd.current_dir(project_root);
    let output_json_file = NamedTempFile::new()?;

    cmd.arg("classify")
        .arg("-i")
        .arg("dummy_input.fa")
        .arg("-d")
        .arg(db_k4_path.to_str().unwrap())
        .arg("-d")
        .arg(db_k3_path.to_str().unwrap())
        // No user k, k from first db (4) will be used, second db (3) will mismatch
        .arg("-o")
        .arg(output_json_file.path());

    cmd.assert().failure().stderr(predicate::str::contains(
        "Effective k-mer size 4 (from first database) does not match k-mer size 3 from database",
    ));
    Ok(())
}

// TODO: Add tests for FASTQ input
// TODO: Add test for empty input file (results in 0s for most things)
// TODO: Add test for no matches
// TODO: Add test for input sequences shorter than k
// TODO: Add test for empty database or reference within a database (how build handles this - it should create an empty set)

#[test]
fn test_classify_min_coverage_filter() -> Result<(), Box<dyn std::error::Error>> {
    let k = 4;
    let temp_db_storage = TempDir::new()?;

    // DB1_REF1_FASTA (ACGT, CGTA, GTAC) - 3 k-mers
    // DB1_REF2_FASTA (GGGA, GGAA, GAAA, AAAA, AAAT, AATT, ATTT, TTTT) -> (CCCA, TTCC, TTTC, AAAA, ATTT, AATT) - 6 unique k-mers
    let db_path = build_db_for_classify(
        k,
        vec![
            ("db_refA.fa", DB1_REF1_FASTA),
            ("db_refB.fa", DB1_REF2_FASTA),
        ],
        &temp_db_storage,
        "db_mincov",
    )?;

    // INPUT_FASTA_BASIC:
    // Unique k-mers: {ACGT, CGTA, GTAC, TTTT, AAAC, CCAA, CCCA, GGGG} - 8 unique
    // Counts: ACGT:2, CGTA:4, GTAC:2, TTTT:1, AAAC:1, CCAA:1, CCCA:1, GGGG:1 (Note: Original test data was slightly off, this is based on re-eval)
    // Corrected: Input unique k-mers (k=4, min_freq=1): {ACGT:4, CGTA:4, GTAC:2, AAAA:1, CAAA:1, CCAA:1, CCCA:1, CCCC:1} (8 unique)

    // Coverage for db_refA (3 k-mers: ACGT, CGTA, GTAC):
    // Matches: ACGT, CGTA, GTAC (3 matches)
    // Breadth: 3/3 = 1.0
    // Coverage for db_refB (6 k-mers: GGGA, AAAA, TTCC, TTTC, ATTT, AATT):
    // Matches: AAAA (1 match)
    // Breadth: 1/6 = 0.1666...

    // Test with min_coverage = 0.5. Expect db_refA to be present, db_refB to be filtered out.
    let results = run_classify_get_json(
        INPUT_FASTA_BASIC,
        "input_mincov.fa",
        &[db_path.clone()],
        Some(k),
        None,      // Default min_kmer_frequency
        Some(0.5), // min_coverage = 0.5
        None,      // No TSV output
    )?;

    assert_eq!(results["databases_analyzed"].as_array().unwrap().len(), 1);
    let db_results = &results["databases_analyzed"][0];
    assert_eq!(
        db_results["references"].as_array().unwrap().len(),
        1,
        "Only one reference should pass min_coverage 0.5"
    );

    let ref_res = &db_results["references"][0];
    assert_eq!(ref_res["reference_name"], "db_refA.fa");
    assert!((ref_res["reference_breadth_of_coverage"].as_f64().unwrap() - 1.0).abs() < 1e-6);

    // Test with min_coverage = 0.1. Expect both references.
    let results_low_cov = run_classify_get_json(
        INPUT_FASTA_BASIC,
        "input_mincov_low.fa",
        &[db_path.clone()],
        Some(k),
        None,
        Some(0.1), // min_coverage = 0.1
        None,
    )?;
    assert_eq!(
        results_low_cov["databases_analyzed"][0]["references"]
            .as_array()
            .unwrap()
            .len(),
        2,
        "Both references should pass min_coverage 0.1"
    );

    Ok(())
}

#[test]
fn test_classify_output_tsv() -> Result<(), Box<dyn std::error::Error>> {
    let k = 4;
    let temp_db_storage = TempDir::new()?;
    let db_path = build_db_for_classify(
        k,
        vec![
            ("db_refA.fa", DB1_REF1_FASTA),
            ("db_refB.fa", DB1_REF2_FASTA),
        ],
        &temp_db_storage,
        "db_tsv",
    )?;

    let temp_output_dir = TempDir::new()?;
    let tsv_output_path = temp_output_dir.path().join("output.tsv");

    let _json_results = run_classify_get_json(
        INPUT_FASTA_BASIC,
        "input_tsv.fa",
        &[db_path.clone()],
        Some(k),
        None,      // Default min_kmer_frequency
        Some(0.5), // min_coverage = 0.5 (so only db_refA appears)
        Some(tsv_output_path.clone()),
    )?;

    assert!(tsv_output_path.exists(), "TSV file was not created");

    let tsv_content = fs::read_to_string(&tsv_output_path)?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(tsv_content.as_bytes());

    let headers = reader.headers()?.clone();
    assert_eq!(
        headers,
        vec![
            "InputFile",
            "Database",
            "Reference",
            "TotalKmersInReference",
            "InputKmersHittingReference",
            "SumDepthMatchedKmers",
            "AvgDepthMatchedKmers",
            "ProportionInputKmersHittingReference",
            "ReferenceBreadthOfCoverage"
        ]
    );

    let records: Vec<csv::StringRecord> = reader.records().collect::<Result<_, _>>()?;
    assert_eq!(
        records.len(),
        1,
        "TSV should contain one data record due to min_coverage filter"
    );

    let record = &records[0];
    assert!(record[0].ends_with("input_tsv.fa"));
    assert!(record[1].contains("db_tsv")); // Check if db_path substring is present
    assert_eq!(&record[2], "db_refA.fa");
    assert_eq!(&record[3], "3"); // TotalKmersInReference for db_refA
    assert_eq!(&record[4], "3"); // InputKmersHittingReference for db_refA
    // SumDepth: ACGT(4) + CGTA(4) + GTAC(2) = 10
    assert_eq!(&record[5], "10"); // SumDepthMatchedKmers for db_refA
    // AvgDepth: 10/3 = 3.3333
    assert_eq!(&record[6], "3.3333");
    // ProportionInputKmersHittingReference: 3 / 8 unique input kmers = 0.375
    assert_eq!(&record[7], "0.3750");
    // ReferenceBreadthOfCoverage: 3 / 3 = 1.0
    assert_eq!(&record[8], "1.0000");

    // Test with no min_coverage (default 0.0), so both references should appear
    let tsv_output_path_all_refs = temp_output_dir.path().join("output_all_refs.tsv");
    let _json_results_all_refs = run_classify_get_json(
        INPUT_FASTA_BASIC,
        "input_tsv_all.fa",
        &[db_path.clone()],
        Some(k),
        None,
        None, // Default min_coverage (0.0)
        Some(tsv_output_path_all_refs.clone()),
    )?;

    let tsv_content_all_refs = fs::read_to_string(&tsv_output_path_all_refs)?;
    let mut reader_all_refs = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(tsv_content_all_refs.as_bytes());
    let records_all_refs: Vec<csv::StringRecord> =
        reader_all_refs.records().collect::<Result<_, _>>()?;
    assert_eq!(
        records_all_refs.len(),
        2,
        "TSV should contain two data records when min_coverage is default"
    );

    // Check second record (db_refB.fa)
    // It matches 1 kmer (AAAA) out of 6 in reference. Depth for AAAA is 1.
    // Total unique input is 8.
    let record_b = records_all_refs
        .iter()
        .find(|r| r[2] == "db_refB.fa")
        .expect("db_refB.fa not found in TSV");
    assert_eq!(&record_b[3], "6"); // TotalKmersInReference for db_refB
    assert_eq!(&record_b[4], "1"); // InputKmersHittingReference for db_refB
    assert_eq!(&record_b[5], "1"); // SumDepthMatchedKmers for db_refB (AAAA depth 1)
    assert_eq!(&record_b[6], "1.0000"); // AvgDepth
    assert_eq!(&record_b[7], format!("{:.4}", 1.0 / 8.0)); // ProportionInputKmersHittingReference
    assert_eq!(&record_b[8], format!("{:.4}", 1.0 / 6.0)); // ReferenceBreadthOfCoverage

    Ok(())
}
