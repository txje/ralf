use std;

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
