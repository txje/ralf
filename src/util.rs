use std;
use bio::io::fasta;
use bio::alphabets;

pub fn bwt_alphabet() -> alphabets::Alphabet {
  alphabets::Alphabet::new(b"$ACGTN")
}
pub fn rank_alphabet() -> alphabets::Alphabet {
  alphabets::Alphabet::new(&[0u8, 1u8, 2u8, 3u8, 4u8, 5u8])
}


#[derive(Debug)]
pub struct Args {
  pub ref_fa: String,
  pub read_fa: String,
  pub k: usize,
  pub min_ordered_matches: usize,
  pub max_match_gap: u32,
  pub max_offset_variance: f32,
  pub max_abs_distance: u32,
  pub min_aln_len: usize,
  pub rep_limit: usize,
  pub bwt: bool
}

pub const USAGE: &'static str = "
Ralf: a rusty long read aligner

Usage:
  rust_aln <ref-fa> <read-fa> [options]

Options:
  -h --help                  Show this screen
  -k                         K-mer size to use for hash matching (default: 16)
  -m --min-ordered-matches   Minimum number of properly ordered hits required to report an alignment (default: 10)
  -g --max-match-gap         Maximum gap between adjacent matches (in query sequence) (default: 1000)
  -x --max-offset-variance   Maximum proportional difference between adjacent query and target match offsets (default: 0.2)
  -n --max-abs-distance      Maximum absolute difference between adjacent hits before proportional gap kicks in (default: 20)
  -l --min-aln-len           Minimum reported alignment length (default: 100)
  -r --rep-limit             Maximum occurrences of a k-mer before it's considered repetitive and ignored (default: 1000)
  --bwt                      Use BWT and FM-index for lookup (defaults to using a HashMap)
";

pub fn help_and_fail(msg:String) {
  if msg.len() > 0 {
    error!("{}", msg);
  }
  println!("{}", USAGE);
  std::process::exit(1);
}


pub fn cat_seq(ref_fa:&str, alphabet:&alphabets::Alphabet) -> (Vec<u8>, Vec<usize>, Vec<fasta::Record>) {
  let mut ref_seqs: Vec<fasta::Record> = Vec::new();

  // Iterate over a FASTA file, use the alphabet to validate read sequences
  let reader = fasta::Reader::from_file(ref_fa).unwrap();
  let mut full_seq: Vec<u8> = Vec::with_capacity(400000000); // approx. 2x total reference sequence length

  let mut multistring_idx: Vec<usize> = Vec::new(); // holds the ordered beginning indices of each reference sequence (and reverse) in the concatenated string

  for record in reader.records() {
    let record = record.unwrap();
    {
      let seq = record.seq();
      let slen = seq.len();

      //info!("Loading reference {}: {} ({}bp)", rid, record.id().unwrap(), seq.len());

      if alphabet.is_word(seq) {
        multistring_idx.push(full_seq.len());
        full_seq.extend_from_slice(seq); // fw
        full_seq.push(36); // ascii for '$'
        multistring_idx.push(full_seq.len());
        full_seq.extend_from_slice(&alphabets::dna::revcomp(seq)); // and rv
        full_seq.push(36);

      } else {
        info!("Reference sequence {} contains something that is not in (A,C,G,T,N). The *entire* sequence is being ignored.", record.id().unwrap());
      }
    }
    ref_seqs.push(record);
  }

  return (full_seq, multistring_idx, ref_seqs);
}
