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

    #[error("K-mer databases have incompatible k-mer sizes: {0} vs {1}")]
    KmerSizeMismatch(u8, u8),

    #[error("An unknown error occurred")]
    Unknown,
}
