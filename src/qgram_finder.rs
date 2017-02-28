extern crate time;

use std::vec::Vec;
use std::collections::HashMap;
use std::ops::Deref;

use bio::io::fasta;
use bio::alphabets;
use bio::data_structures::qgram_index;

//use bwt::{bwt, less, Occ, BWT, Less};
use overlapper::{Overlapper, Position, KmerMatch};
use util::{cat_seq};

/* Full k-mer index using rust-bio's QGramIndex
 * - the sequence will be bit-packed using the given alphabet
 *    since our data may include Ns, we use the ACGTN alphabet,
 *    which must be encoded in 3 bits instead of 2
 */

pub struct QGramFinder {
  qgi: QGramIndex,
  sequences: Vec<fasta::Record>,
  k: usize,
  rep_limit: usize
}

impl QGramFinder {
  pub fn new(ref_fa:&str, k:usize, rep_limit:usize, alphabet:&alphabets::Alphabet) -> QGramFinder {

    let t_0:u64 = time::precise_time_ns();

    info!("Concatenating reference sequence from {}", ref_fa);
    let (cat_seq, seq_idx, ref_seqs) = cat_seq(&ref_fa[..], &alphabet);
    let t_1:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_1 - t_0) as f64 / 1000000000.0);

    info!("Building Q-gram index", ref_fa);
    let qgram_index = qgram_index::QGramIndex::with_max_count(k, &cat_seq, &alphabet, rep_limit);
    let t_2:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_2 - t_1) as f64 / 1000000000.0);

    QGramFinder{qgi:qgram_index, sequences:ref_seqs, k:k, rep_limit:rep_limit}
  }
}

impl Overlapper for QGramFinder {
  fn ovlSeq(&self, seq:&[u8]) -> HashMap<u32,Vec<KmerMatch>> {
    let mut hit_hash:HashMap<u32,Vec<KmerMatch>> = HashMap::new();

    for i in 0..(seq.len()-self.k+1) {
      let kmer = &seq[i..i+self.k]; // just using string slice for the bwt

      let matches = qgram_index.matches(kmer, 1);

      // convert to Vec<Position>
      let pos_vec:Vec<Position> = Vec::new();

      if pos_vec.len() < self.rep_limit {
        for p in 0..pos_vec.len() {
          if !hit_hash.contains_key(&pos_vec[p].rid) {
            let match_vec = Vec::new();
            hit_hash.insert(pos_vec[p].rid, match_vec);
          }
          let mut match_vec = hit_hash.get_mut(&pos_vec[p].rid).unwrap();
          match_vec.push(KmerMatch{qpos:i as u32, tpos: pos_vec[p].pos});
        }
      }
    }

    hit_hash
  }

  fn sequences(&self) -> &Vec<fasta::Record> {
    &self.sequences
  }
}
