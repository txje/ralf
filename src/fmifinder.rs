extern crate time;

use std::vec::Vec;
use std::collections::HashMap;
use std::ops::Deref;
use std::rc::Rc;

use bio::io::fasta;
use bio::alphabets;
use bio::data_structures::fmindex::{FMIndex, FMDIndex};
use bio::data_structures::suffix_array::{RawSuffixArray, SuffixArray, suffix_array};
use bio::data_structures::bwt::{bwt, less, Occ, BWT, Less};

//use bwt::{bwt, less, Occ, BWT, Less};
use overlapper::{Overlapper, Position, KmerMatch};
use util::{cat_seq};


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

    // this all has to be built here for unavoidable (I think) lifetime issues
    // (the FMIndex *requires* references to bwt,less,occ, which necessarily go out of scope)
    info!("Building suffix array");
    let sa = suffix_array(&cat_seq); // slice the [u8] out of a Vec<u8>, I think
    let t_2:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_2 - t_1) as f64 / 1000000000.0);

    info!("Building BWT");
    let bwt = bwt(&cat_seq, &sa);
    let t_3:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_3 - t_2) as f64 / 1000000000.0);

    info!("Computing less");
    let less = less(&bwt, &alphabet);
    let t_4:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_4 - t_3) as f64 / 1000000000.0);

    info!("Computing occurrences");
    let occ = Occ::new(&bwt, 3, &alphabet); // sampling rate of every 3rd position
    let t_5:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_5 - t_4) as f64 / 1000000000.0);

    info!("Building FM-index");
    //let fm = FMIndex::new(DerefBox{value: bwt}, DerefBox{value: less}, DerefBox{value: occ});
    let fm = FMIndex::new(Rc::new(bwt), Rc::new(less), Rc::new(occ));
    let fmd = FMDIndex::from(fm);
    let t_6:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_6 - t_5) as f64 / 1000000000.0);

    FMIFinder{sa:sa, fmd:fmd, sequences:ref_seqs, k:k, rep_limit:rep_limit}
  }
}

impl Overlapper for FMIFinder {
  fn ovlSeq(&self, seq:&[u8]) -> HashMap<u32,Vec<KmerMatch>> {
    let mut hit_hash:HashMap<u32,Vec<KmerMatch>> = HashMap::new();

    for i in 0..(seq.len()-self.k+1) {
      let kmer = &seq[i..i+self.k]; // just using string slice for the bwt

      let intervals = self.fmd.smems(kmer, 2);
      let forward_positions = intervals[0].forward().occ(&self.sa);
      let revcomp_positions = intervals[0].revcomp().occ(&self.sa);

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
