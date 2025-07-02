use std::path::PathBuf; // Added for KmerSizeMismatchValidation
use thiserror::Error;

#[derive(Error, Debug)]
pub enum OrionKmerError {
    #[error("Invalid K-mer size: {0}. Must be between 1 and 32.")]
    InvalidKmerSize(u8),

    #[error("File not found: {0}")]
    FileNotFound(String),

    #[error("Failed to parse input file: {0}")]
    FileParsingError(String),

    #[error("I/O error")]
    IoError(#[from] std::io::Error),

    #[error("Serialization error: {0}")]
    SerializationError(String),

    #[error("Deserialization error: {0}")]
    DeserializationError(String),

    #[error("K-mer databases have incompatible k-mer sizes (overall comparison): {0} vs {1}")]
    KmerSizeMismatch(u8, u8),

    #[error("User-provided k-mer size {0} does not match k-mer size {1} from database: {2:?}")]
    KmerSizeMismatchValidation(u8, u8, PathBuf),

    #[error(
        "Effective k-mer size {0} (from first database) does not match k-mer size {1} from database: {2:?}"
    )]
    KmerSizeMismatchBetweenDatabases(u8, u8, PathBuf), // Specific for classify

    #[error("Generic error: {0}")]
    Generic(String),

    #[error("An unknown error occurred")]
    Unknown,
}
