extern crate time;

use std::vec::Vec;
use std::collections::HashMap;
use std::ops::Deref;
use std::rc::Rc;

use bio::io::fasta;
use bio::alphabets;
use bio::alphabets::dna::revcomp;
use bio::alphabets::RankTransform;
use bio::data_structures::fmindex::{FMIndex, FMDIndex};
use bio::data_structures::suffix_array::{RawSuffixArray, SuffixArray, suffix_array};
use bio::data_structures::bwt::{bwt, less, Occ, BWT, Less};

//use bwt::{bwt, less, Occ, BWT, Less};
use overlapper::{Overlapper, Position, KmerMatch};
use util::{cat_seq, bwt_alphabet, rank_alphabet};


/* -----------------------------------------------------------------
 * BWT/FM-index k-mer lookup machine
 * -----------------------------------------------------------------*/

pub struct FMIFinder {
  sa: RawSuffixArray,
  // Rc<T> is a single-thread reference-counted pointer
  // -- critically, it takes ownership of its contents
  // -- use an Arc for thread-safe atomic reference counting (around the FMIIndex if need be)
  fmd: FMDIndex<Rc<BWT>,Rc<Less>,Rc<Occ>>,
  sequences: Vec<fasta::Record>,
  ref_seq: Vec<u8>,
  ref_idx: Vec<usize>,
  k: usize,
  rep_limit: usize
}

impl FMIFinder {
  pub fn new(ref_fa:&str, k:usize, rep_limit:usize, alphabet:&alphabets::Alphabet) -> FMIFinder {

    let t_0:u64 = time::precise_time_ns();

    info!("Concatenating reference sequence from {}", ref_fa);
    let (cat_seq, seq_idx, ref_seqs) = cat_seq(&ref_fa[..], &alphabet);
    let t_1:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_1 - t_0) as f64 / 1000000000.0);

    info!("Computing rank transform");
    let rank_trans = RankTransform::new(&bwt_alphabet());
    let rank_alphabet = rank_alphabet();
    let seq = rank_trans.transform(&cat_seq);
    let t_2:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_2 - t_1) as f64 / 1000000000.0);

    // this all has to be built here for unavoidable (I think) lifetime issues
    // (the FMIndex *requires* references to bwt,less,occ, which necessarily go out of scope)
    info!("Building suffix array");
    let sa = suffix_array(&seq); // slice the &[u8] out of a Vec<u8>
    let t_3:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_3 - t_2) as f64 / 1000000000.0);

    info!("Building BWT");
    let bwt = bwt(&seq, &sa);
    let t_4:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_4 - t_3) as f64 / 1000000000.0);

    info!("Computing less");
    let less = less(&bwt, &rank_alphabet);
    let t_5:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_5 - t_4) as f64 / 1000000000.0);

    info!("Computing occurrences -- should take up O(n / 3[sampling] * 6[characters]) space");
    let occ = Occ::new(&bwt, 3, &rank_alphabet); // sampling rate of every 3rd position
    let t_6:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_6 - t_5) as f64 / 1000000000.0);

    info!("Building FM-index");
    let fm = FMIndex::new(Rc::new(bwt), Rc::new(less), Rc::new(occ), rank_trans); // Rc (single-thread reference counter) gains ownership of BWT components, derefs and clones
    let fmd = FMDIndex::from(fm);
    let t_7:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_7 - t_6) as f64 / 1000000000.0);

    FMIFinder{sa:sa, fmd:fmd, sequences:ref_seqs, ref_seq:cat_seq, ref_idx:seq_idx, k:k, rep_limit:rep_limit}
  }
}

impl Overlapper for FMIFinder {
  fn ovlSeq(&self, seq:&[u8]) -> HashMap<u32,Vec<KmerMatch>> {
    let mut hit_hash:HashMap<u32,Vec<KmerMatch>> = HashMap::new();

    for i in 0..(seq.len()-self.k+1) {
      let kmer = &seq[i..i+self.k]; // just using string slice for the bwt

      let intervals = self.fmd.smems_transform(kmer, 2);
      let forward_positions = intervals[0].forward().occ(&self.sa);
      let revcomp_positions = intervals[0].revcomp().occ(&self.sa);

      debug!("fw: {:?}", forward_positions);
      debug!("rv: {:?}", revcomp_positions);
      debug!("kmer {:?} matched fw[0] {:?} and rv[0] {:?}", kmer, &self.ref_seq[forward_positions[0]..forward_positions[0]+self.k], &revcomp(&self.ref_seq[revcomp_positions[0]..revcomp_positions[0]+self.k]));

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
