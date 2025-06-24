use clap::{Arg, Command};
use rust_htslib::{
    bam::{self, Read},
    faidx,
};
use serde::Serialize;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;

/// We store up to 6 bases (A, C, G, T, N, '-') in 3 bits. That means
/// k <= floor(64 / 3) = 21 if we pack them into a 64-bit integer.
const ALPHABET_SIZE: usize = 6; // A, C, G, T, N, '-'

/// Encodes a base into 3 bits, using the 6-character alphabet {A,C,G,T,N,-}.
/// Everything outside ACGTN -> '-'
fn encode_base_6(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        b'N' | b'n' => 4,
        _ => 5, // '-'
    }
}

/// Decodes a 64-bit integer k-mer (with 3 bits/base) back to a String.
fn decode_kmer_6(mut encoded: u64, k: usize) -> String {
    // codes 0..5 => [A, C, G, T, N, '-']
    let decode_table = ['A', 'C', 'G', 'T', 'N', '-', '?', '?'];
    let mut out = vec![b'?'; k];

    for i in (0..k).rev() {
        let code = (encoded & 0b111) as usize;
        encoded >>= 3;
        out[i] = decode_table[code] as u8;
    }
    String::from_utf8(out).unwrap()
}

/// Increments k-mer counts for all k-mers in `seq`.
/// Uses a rolling 3-bit representation for A,C,G,T,N,-.
fn increment_kmers_6(seq: &[u8], k: usize, counts: &mut HashMap<u64, u64>) {
    if seq.len() < k {
        return;
    }

    // Build initial k-mer
    let mut kmer: u64 = 0;
    for &base in &seq[..k] {
        let code = encode_base_6(base) as u64;
        kmer = (kmer << 3) | code;
    }
    *counts.entry(kmer).or_insert(0) += 1;

    // Mask to keep lower 3*k bits
    let mask: u64 = (1 << (3 * k)) - 1;

    // Slide over the sequence
    for &base in &seq[k..] {
        let code = encode_base_6(base) as u64;
        kmer = ((kmer << 3) & mask) | code;
        *counts.entry(kmer).or_insert(0) += 1;
    }
}

/// Process the reference FASTA using faidx, counting k-mers for each contig.
/// We read each contig in full, fetch its sequence, then increment k-mer counts.
fn count_kmers_in_fasta(
    fasta_path: &str,
    k: usize,
) -> Result<HashMap<u64, u64>, Box<dyn std::error::Error>> {
    let faidx = faidx::Reader::from_path(fasta_path)?;
    let nseqs = faidx.n_seqs() as i32;

    let mut ref_counts = HashMap::new();
    for i in 0..nseqs {
        let name = faidx.seq_name(i)?;
        let seqlen = faidx.fetch_seq_len(&name);
        // fetch the entire contig
        let seq = faidx.fetch_seq(&name, 0, seqlen as usize)?;
        // increment counts
        increment_kmers_6(&seq, k, &mut ref_counts);
    }
    Ok(ref_counts)
}

/// Process a BAM file, counting k-mers from each read.
fn count_kmers_in_bam(
    bam_path: &str,
    k: usize,
) -> Result<HashMap<u64, u64>, Box<dyn std::error::Error>> {
    let mut bam_reader = bam::Reader::from_path(bam_path)?;

    let mut bam_counts = HashMap::new();
    for record_result in bam_reader.records() {
        let record = record_result?;
        let flag = record.flags();
        // Skip unmapped, duplicates, supplementary, low mapQ, etc.
        if (flag & 0x4) != 0 || (flag & 0x400) != 0 || (flag & 0x800) != 0 {
            continue;
        }
        if record.mapq() == 0 {
            continue;
        }
        // convert seq to bytes
        let seq = record.seq().as_bytes();
        increment_kmers_6(&seq, k, &mut bam_counts);
    }
    Ok(bam_counts)
}

/// This struct holds the final proportion and enrichment data for each k-mer.
#[derive(Serialize)]
struct KmerProportion {
    kmer: String,
    ref_prop: f64,
    bam_prop: f64,
    enrichment: Option<f64>,
}

/// Generate a vector of KmerProportion, merging the reference and BAM HashMaps
/// and computing proportions & enrichment.
///
/// - ref_prop = ref_count / total_ref_kmers
/// - bam_prop = bam_count / total_bam_kmers
/// - enrichment = bam_prop / ref_prop   (if ref_prop = 0 => handle infinite or 0)
fn generate_kmer_proportions(
    ref_counts: HashMap<u64, u64>,
    bam_counts: HashMap<u64, u64>,
    k: usize,
) -> Vec<KmerProportion> {
    let total_ref_kmers: u64 = ref_counts.values().sum();
    let total_bam_kmers: u64 = bam_counts.values().sum();

    // Collect all unique k-mers from ref_counts and bam_counts
    let mut all_keys: Vec<u64> = ref_counts
        .keys()
        .chain(bam_counts.keys())
        .copied()
        .collect();
    all_keys.sort_unstable();
    all_keys.dedup();

    let mut output = Vec::with_capacity(all_keys.len());

    for key in all_keys {
        let ref_count = *ref_counts.get(&key).unwrap_or(&0);
        let bam_count = *bam_counts.get(&key).unwrap_or(&0);

        let ref_prop = if total_ref_kmers > 0 {
            ref_count as f64 / total_ref_kmers as f64
        } else {
            0.0
        };
        let bam_prop = if total_bam_kmers > 0 {
            bam_count as f64 / total_bam_kmers as f64
        } else {
            0.0
        };

        // If `ref_prop` is zero but `bam_prop > 0`, we might say it's infinity
        // or we can skip it. For demonstration, let's do "infinity" if bam_prop>0:
        let enrichment = if ref_prop > 0.0 {
            Some(bam_prop / ref_prop)
        } else if bam_prop > 0.0 {
            None
        } else {
            None
        };

        output.push(KmerProportion {
            kmer: decode_kmer_6(key, k),
            ref_prop,
            bam_prop,
            enrichment,
        });
    }

    // Sort by descending enrichment, with `None` last
    output.sort_by(|a, b| {
        match (a.enrichment, b.enrichment) {
            // both Some: do a numeric descending comparison
            (Some(ae), Some(be)) => be.partial_cmp(&ae).unwrap_or(Ordering::Equal),
            // a=Some, b=None => a should come first => Ordering::Less
            (Some(_), None) => Ordering::Less,
            // a=None, b=Some => b should come first => Ordering::Greater
            (None, Some(_)) => Ordering::Greater,
            // both None => equal
            (None, None) => Ordering::Equal,
        }
    });

    output
}

/// Write the vector of `KmerProportion` as JSON to `output_path`.
fn write_kmer_proportions_json(
    kmer_props: &[KmerProportion],
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let json = serde_json::to_string_pretty(&kmer_props)?;
    let mut file = File::create(output_path)?;
    file.write_all(json.as_bytes())?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("kmer_enrichment")
        .version("0.1.0")
        .author("Harry Pink")
        .about("Compute k-mer enrichment using a reference FASTA (via faidx) and a BAM")
        .arg(
            Arg::new("ref")
                .long("ref")
                .short('r')
                .value_name("FASTA")
                .help("Path to reference FASTA file (must be indexed: .fai)")
                .required(true),
        )
        .arg(
            Arg::new("bam")
                .long("bam")
                .short('b')
                .value_name("BAM")
                .help("Path to input BAM file")
                .required(true),
        )
        .arg(
            Arg::new("k")
                .long("kmer")
                .short('k')
                .value_name("K")
                .help("k-mer size (3*k <= 64, so typically k <= 21)")
                .required(true),
        )
        .arg(
            Arg::new("output")
                .long("output")
                .short('o')
                .value_name("OUTPUT")
                .help("Where to write JSON output")
                .required(true),
        )
        .get_matches();

    let ref_path = matches.get_one::<String>("ref").unwrap();
    let bam_path = matches.get_one::<String>("bam").unwrap();
    let k: usize = matches
        .get_one::<String>("k")
        .unwrap()
        .parse()
        .expect("Failed to parse k");
    let output_path = matches.get_one::<String>("output").unwrap();

    if 3 * k > 64 {
        eprintln!(
            "ERROR: k={} is too large for 6-base encoding in 64 bits (max ~21).",
            k
        );
        std::process::exit(1);
    }

    // 1) Count k-mers in the reference
    println!("Counting k-mers in reference FASTA...");
    let ref_kmer_counts =
        count_kmers_in_fasta(ref_path, k).expect("Failed to count k-mers in reference FASTA");

    // 2) Count k-mers in BAM
    println!("Counting k-mers in BAM...");
    let bam_kmer_counts = count_kmers_in_bam(bam_path, k).expect("Failed to count k-mers in BAM");

    // 3) Generate the KmerProportion vector
    println!("Calculating k-mer enrichment...");
    let kmer_props = generate_kmer_proportions(ref_kmer_counts, bam_kmer_counts, k);

    // Partition into two groups based on whether the `kmer` has 'N' or '-'
    let (canonical, with_n_or_dash): (Vec<KmerProportion>, Vec<KmerProportion>) =
        kmer_props.into_iter().partition(|kp| {
            // Return `true` if it's canonical => no 'N', no '-'
            // Return `false` if it should go to the "with_n_or_dash" list
            !(kp.kmer.contains('N') || kp.kmer.contains('-'))
        });

    // 4) Write JSON
    println!("Writing JSON output to {}", output_path);
    write_kmer_proportions_json(&canonical, "canonical.json")?;
    write_kmer_proportions_json(&with_n_or_dash, "with_n.json")?;

    println!("Done.");

    Ok(())
}
