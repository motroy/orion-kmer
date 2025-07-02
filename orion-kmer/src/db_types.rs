use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Represents a k-mer database.
///
/// The database stores k-mers associated with reference identifiers (e.g., filenames).
#[derive(Serialize, Deserialize, Debug, Clone)] // Added Clone
pub struct KmerDbV2 {
    /// The k-mer size used in this database.
    pub k: u8,
    /// A map where keys are reference identifiers (e.g., filenames from which k-mers were derived)
    /// and values are sets of unique k-mers (encoded as u64) found in that reference.
    pub references: HashMap<String, HashSet<u64>>,
}

// For compatibility, and for commands that might operate on a "flat" DB view
// or when loading older DB formats (though explicit versioning/migration is better for old formats)
// For now, let's assume build.rs will define its old KmerDb struct for loading if needed,
// or we introduce proper versioning.
// This KmerDb can be used by compare/query if they operate on a union of kmers.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct KmerDb {
    pub k: u8,
    pub kmers: HashSet<u64>,
}

impl KmerDbV2 {
    /// Creates a new, empty KmerDbV2 with a specified k-mer size.
    pub fn new(k: u8) -> Self {
        KmerDbV2 {
            k,
            references: HashMap::new(),
        }
    }

    /// Adds a reference and its set of k-mers to the database.
    /// If the reference name already exists, its k-mer set will be overwritten.
    pub fn add_reference(&mut self, name: String, kmers: HashSet<u64>) {
        self.references.insert(name, kmers);
    }

    /// Returns a unified set of all unique k-mers from all references in the database.
    pub fn get_all_kmers_unified(&self) -> HashSet<u64> {
        self.references
            .values()
            .flat_map(|kmer_set| kmer_set.iter().copied())
            .collect()
    }

    /// Returns the total number of unique k-mers across all references.
    pub fn total_unique_kmers(&self) -> usize {
        self.get_all_kmers_unified().len()
    }

    /// Returns the number of references stored in the database.
    pub fn num_references(&self) -> usize {
        self.references.len()
    }
}
