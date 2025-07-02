// K-mer processing logic

const BITS_PER_BASE: u8 = 2;

/// Encodes a single DNA base into its 2-bit representation.
/// A -> 00 (0)
/// C -> 01 (1)
/// G -> 10 (2)
/// T -> 11 (3)
/// Returns None if the base is not A, C, G, or T.
#[inline]
fn dna_base_to_u64(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0b00),
        b'C' | b'c' => Some(0b01),
        b'G' | b'g' => Some(0b10),
        b'T' | b't' => Some(0b11),
        _ => None,
    }
}

/// Decodes a 2-bit representation back to a DNA base character.
#[inline]
fn u64_to_dna_base(val: u64) -> u8 {
    match val {
        0b00 => b'A',
        0b01 => b'C',
        0b10 => b'G',
        0b11 => b'T',
        _ => panic!("Invalid 2-bit value for DNA base"), // Should not happen if encoding is correct
    }
}

/// Encodes a DNA sequence slice into a u64.
/// K-mer length `k` must be between 1 and 32, inclusive.
/// Returns `None` if the sequence contains non-ACGT characters or if k is invalid.
pub fn seq_to_u64(seq: &[u8], k: u8) -> Option<u64> {
    if k == 0 || k > 32 {
        return None; // k-mer length out of bounds for u64 representation
    }
    if seq.len() != k as usize {
        return None; // Sequence length does not match k
    }

    let mut kmer_val: u64 = 0;
    for (i, &base) in seq.iter().enumerate() {
        match dna_base_to_u64(base) {
            Some(base_val) => {
                // Shift previous bits left by 2 and add new base's bits
                // The first base will be at the most significant bits of the relevant part
                kmer_val |= base_val << (BITS_PER_BASE * (k - 1 - i as u8));
            }
            None => return None, // Invalid character in sequence
        }
    }
    Some(kmer_val)
}

/// Decodes a u64 k-mer representation back to a DNA sequence (Vec<u8>).
/// `k` specifies the length of the k-mer.
pub fn u64_to_seq(kmer_val: u64, k: u8) -> Vec<u8> {
    if k == 0 || k > 32 {
        panic!("Invalid k-mer length for decoding: {}", k);
    }
    let mut seq = Vec::with_capacity(k as usize);
    // Create a mask to extract 2 bits at a time
    let mask = 0b11;
    for i in 0..k {
        // Extract the i-th base from the most significant part
        let shift = BITS_PER_BASE * (k - 1 - i);
        let base_val = (kmer_val >> shift) & mask;
        seq.push(u64_to_dna_base(base_val));
    }
    seq
}

/// Computes the reverse complement of an encoded k-mer.
/// `k` is the length of the k-mer.
pub fn reverse_complement_u64(kmer_val: u64, k: u8) -> u64 {
    if k == 0 || k > 32 {
        panic!("Invalid k-mer length for reverse complement: {}", k);
    }
    let mut rc_val: u64 = 0;
    let mask = 0b11; // To get 2 bits
    for i in 0..k {
        // Get the i-th base from the original k-mer (from LSB side for easier processing here)
        let base_val = (kmer_val >> (BITS_PER_BASE * i)) & mask;
        // Complement: A(00) <-> T(11), C(01) <-> G(10). This is equivalent to XOR with 0b11 (3).
        let complemented_base = base_val ^ 0b11;
        // Place the complemented base in the reverse position
        rc_val |= complemented_base << (BITS_PER_BASE * (k - 1 - i));
    }
    rc_val
}

/// Returns the canonical representation of a k-mer.
/// The canonical k-mer is the lexicographically smaller of the k-mer and its reverse complement.
/// `k` is the length of the k-mer.
pub fn canonical_u64(kmer_val: u64, k: u8) -> u64 {
    let rc_kmer_val = reverse_complement_u64(kmer_val, k);
    if kmer_val < rc_kmer_val {
        kmer_val
    } else {
        rc_kmer_val
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_problematic_classify_kmers() {
        let k = 4;
        // Input sequence TTTTGGGG
        let tttt_val = seq_to_u64(b"TTTT", k).unwrap(); // 255
        let tttt_canon = canonical_u64(tttt_val, k);   // Should be AAAA (0)
        let aaaa_val = seq_to_u64(b"AAAA", k).unwrap();   // 0
        assert_eq!(tttt_canon, aaaa_val, "Canonical of TTTT should be AAAA value (0)");

        let tggg_val = seq_to_u64(b"TGGG", k).unwrap(); // 234 (11101010)
        let tggg_canon = canonical_u64(tggg_val, k);   // RC is CCCA (01010100 = 84). So canon is CCCA (84).
        let ccca_val = seq_to_u64(b"CCCA", k).unwrap();   // 84
        assert_eq!(tggg_canon, ccca_val, "Canonical of TGGG should be CCCA value (84)");

        // DB sequence GGGAAAAATTTT
        let ggga_val = seq_to_u64(b"GGGA", k).unwrap(); // 168 (10101000)
        let ggga_canon = canonical_u64(ggga_val, k);   // RC is TCCC (11010101 = 213). So canon is GGGA (168).
        assert_eq!(ggga_canon, ggga_val, "Canonical of GGGA (168) should be GGGA (168)");
        // This was the line that failed: assert_eq!(ggga_canon, ccca_val) which is 168 == 84 (false)

        // Critical checks based on correct canonicals:
        // 1. Input TTTT (canon AAAA) vs DB AAAA (canon AAAA)
        let db_aaaa_from_db_seq_val = seq_to_u64(b"AAAA", k).unwrap();
        let db_aaaa_from_db_seq_canon = canonical_u64(db_aaaa_from_db_seq_val, k);
        assert_eq!(tttt_canon, db_aaaa_from_db_seq_canon, "Canon(TTTT) from input should match Canon(AAAA) from DB");

        // 2. Input TGGG (canon CCCA=84) vs DB GGGA (canon GGGA=168)
        // These should NOT match.
        assert_ne!(tggg_canon, ggga_canon, "Canon(TGGG) from input should NOT match Canon(GGGA) from DB");
    }


    #[test]
    fn test_dna_base_to_u64() {
        assert_eq!(dna_base_to_u64(b'A'), Some(0b00));
        assert_eq!(dna_base_to_u64(b'c'), Some(0b01));
        assert_eq!(dna_base_to_u64(b'G'), Some(0b10));
        assert_eq!(dna_base_to_u64(b't'), Some(0b11));
        assert_eq!(dna_base_to_u64(b'N'), None);
        assert_eq!(dna_base_to_u64(b'X'), None);
    }

    #[test]
    fn test_u64_to_dna_base() {
        assert_eq!(u64_to_dna_base(0b00), b'A');
        assert_eq!(u64_to_dna_base(0b01), b'C');
        assert_eq!(u64_to_dna_base(0b10), b'G');
        assert_eq!(u64_to_dna_base(0b11), b'T');
    }

    #[test]
    #[should_panic]
    fn test_u64_to_dna_base_panic() {
        u64_to_dna_base(0b100); // Invalid input
    }

    #[test]
    fn test_seq_to_u64_valid() {
        assert_eq!(seq_to_u64(b"A", 1), Some(0b00)); // A
        assert_eq!(seq_to_u64(b"C", 1), Some(0b01)); // C
        assert_eq!(seq_to_u64(b"G", 1), Some(0b10)); // G
        assert_eq!(seq_to_u64(b"T", 1), Some(0b11)); // T

        // K=3, ACG = 000110
        assert_eq!(seq_to_u64(b"ACG", 3), Some(0b000110));
        // K=4, ACGT = 00011011
        assert_eq!(seq_to_u64(b"ACGT", 4), Some(0b00011011));
        // K=5, TTTTT = 1111111111
        assert_eq!(seq_to_u64(b"TTTTT", 5), Some(0b1111111111));
        // K=32, AAAA... (32 times)
        let k32_a = vec![b'A'; 32];
        assert_eq!(seq_to_u64(&k32_a, 32), Some(0));
        // K=32, TTTT... (32 times) -> all 64 bits set to 1
        let k32_t = vec![b'T'; 32];
        assert_eq!(seq_to_u64(&k32_t, 32), Some(u64::MAX));

        // Check case-insensitivity
        assert_eq!(seq_to_u64(b"acg", 3), Some(0b000110));
    }

    #[test]
    fn test_seq_to_u64_invalid_char() {
        assert_eq!(seq_to_u64(b"ACN", 3), None);
        assert_eq!(seq_to_u64(b"NA", 2), None);
        assert_eq!(seq_to_u64(b"X", 1), None);
    }

    #[test]
    fn test_seq_to_u64_invalid_k() {
        assert_eq!(seq_to_u64(b"A", 0), None); // k=0
        assert_eq!(seq_to_u64(b"A", 33), None); // k=33
        assert_eq!(seq_to_u64(b"ACG", 2), None); // seq.len != k
        assert_eq!(seq_to_u64(b"A", 2), None); // seq.len != k
    }

    #[test]
    fn test_u64_to_seq_valid() {
        assert_eq!(u64_to_seq(0b000110, 3), b"ACG".to_vec());
        assert_eq!(u64_to_seq(0b00011011, 4), b"ACGT".to_vec());
        assert_eq!(u64_to_seq(0b1111111111, 5), b"TTTTT".to_vec());
        assert_eq!(u64_to_seq(0, 1), b"A".to_vec());
        assert_eq!(u64_to_seq(0, 32), vec![b'A'; 32]);
        assert_eq!(u64_to_seq(u64::MAX, 32), vec![b'T'; 32]);
    }

    #[test]
    #[should_panic]
    fn test_u64_to_seq_invalid_k_zero() {
        u64_to_seq(0, 0);
    }

    #[test]
    #[should_panic]
    fn test_u64_to_seq_invalid_k_large() {
        u64_to_seq(0, 33);
    }

    #[test]
    fn test_reverse_complement_u64() {
        // K=1
        assert_eq!(
            reverse_complement_u64(seq_to_u64(b"A", 1).unwrap(), 1),
            seq_to_u64(b"T", 1).unwrap()
        ); // A -> T
        assert_eq!(
            reverse_complement_u64(seq_to_u64(b"T", 1).unwrap(), 1),
            seq_to_u64(b"A", 1).unwrap()
        ); // T -> A
        assert_eq!(
            reverse_complement_u64(seq_to_u64(b"C", 1).unwrap(), 1),
            seq_to_u64(b"G", 1).unwrap()
        ); // C -> G
        assert_eq!(
            reverse_complement_u64(seq_to_u64(b"G", 1).unwrap(), 1),
            seq_to_u64(b"C", 1).unwrap()
        ); // G -> C

        // K=3, ACG (000110) -> CGT (011011)
        // A (00) -> T (11)
        // C (01) -> G (10)
        // G (10) -> C (01)
        // ACG -> TGC (rc string), encoding of CGT
        // CGT: C(01) G(10) T(11) -> 011011
        assert_eq!(
            reverse_complement_u64(seq_to_u64(b"ACG", 3).unwrap(), 3),
            seq_to_u64(b"CGT", 3).unwrap()
        );

        // K=4, ATGC (00111001) -> GCAT (10010011)
        assert_eq!(
            reverse_complement_u64(seq_to_u64(b"ATGC", 4).unwrap(), 4),
            seq_to_u64(b"GCAT", 4).unwrap()
        );

        // Palindrome K=4, ATTA (00111100) -> TAAT (11000011)
        // This is not a self-complement palindrome.
        // ATTA -> TAAT
        assert_eq!(
            reverse_complement_u64(seq_to_u64(b"ATTA", 4).unwrap(), 4),
            seq_to_u64(b"TAAT", 4).unwrap()
        );

        // Self-complement palindrome K=4, GTAC (10110001)
        // GTAC -> GTAC
        assert_eq!(
            reverse_complement_u64(seq_to_u64(b"GTAC", 4).unwrap(), 4),
            seq_to_u64(b"GTAC", 4).unwrap()
        );

        // K=5, AAAAA (0000000000) -> TTTTT (1111111111)
        assert_eq!(
            reverse_complement_u64(seq_to_u64(b"AAAAA", 5).unwrap(), 5),
            seq_to_u64(b"TTTTT", 5).unwrap()
        );
    }

    #[test]
    #[should_panic]
    fn test_reverse_complement_invalid_k_zero() {
        reverse_complement_u64(0, 0);
    }

    #[test]
    #[should_panic]
    fn test_reverse_complement_invalid_k_large() {
        reverse_complement_u64(0, 33);
    }

    #[test]
    fn test_canonical_u64() {
        let k = 3;
        // ACG (000110) vs CGT (011011). ACG is smaller.
        let kmer_acg = seq_to_u64(b"ACG", k).unwrap();
        assert_eq!(canonical_u64(kmer_acg, k), kmer_acg);

        // TGT (111011) vs ACA (000100). ACA is smaller.
        let kmer_tgt = seq_to_u64(b"TGT", k).unwrap();
        let kmer_aca = seq_to_u64(b"ACA", k).unwrap();
        assert_eq!(canonical_u64(kmer_tgt, k), kmer_aca);

        // Self-complement palindrome k=4
        // GTAC (10110001) -> RC is GTAC. Canonical is GTAC.
        let k_pal = 4;
        let kmer_gtac = seq_to_u64(b"GTAC", k_pal).unwrap();
        assert_eq!(canonical_u64(kmer_gtac, k_pal), kmer_gtac);

        // Check another case
        // K=5, GATTC (1000111101) vs GAATC (1000001101)
        // GATTC -> reverse_complement is GAATC
        // GAATC is smaller.
        let k_5 = 5;
        let kmer_gattc = seq_to_u64(b"GATTC", k_5).unwrap();
        let kmer_gaatc = seq_to_u64(b"GAATC", k_5).unwrap();
        assert_eq!(canonical_u64(kmer_gattc, k_5), kmer_gaatc);
    }
}
