extern crate time;

use std::vec::Vec;
use std::collections::HashMap;

use bio::io::fasta;
use bio::alphabets;
use bio::data_structures::fmindex::{FMIndex, FMDIndex};
use bio::data_structures::suffix_array::{RawSuffixArray, SuffixArray, suffix_array};
use bio::data_structures::bwt::{bwt, less, Occ, BWT, Less};

use overlapper::{Overlapper, Position, KmerMatch};
use util::{Args,help_and_fail};


/* -----------------------------------------------------------------
 * BWT/FM-index k-mer lookup machine
 * -----------------------------------------------------------------*/

pub struct FMIFinder<'a> {
  sa: RawSuffixArray,
  bwt: BWT,
  less: Less,
  occ: Occ,
  fmd: FMDIndex<&'a BWT,&'a Less,&'a Occ>,
  args:&'a Args
}

impl<'a> FMIFinder<'a> {
  pub fn new(ref_fa:&str, args:&'a Args, alphabet:&alphabets::Alphabet) -> FMIFinder<'a> {
    // try to intialize an "empty" FMIFinder, so that we can fill it as we make things
    let mut fmi:FMIFinder = FMIFinder{
      sa:RawSuffixArray::new(),
      bwt:BWT::new(),
      less:Less::new(),
      occ:Occ::new(&BWT::new(), 0, &alphabet),
      fmd:FMDIndex::from(FMIndex::new(&BWT::new(),&Less::new(),&Occ::new(&BWT::new(), 0, &alphabet))),
      args:args
    };
    fmi
  }

  pub fn make(&'a mut self, ref_fa:&str, args:&'a Args, alphabet:&alphabets::Alphabet) {
    let t_0:u64 = time::precise_time_ns();

    info!("Concatenating reference sequence from {}", args.ref_fa);
    let (cat_seq, seq_idx, ref_seqs) = load_seq(&args.ref_fa[..], &args, &alphabet);
    let t_1:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_1 - t_0) as f64 / 1000000000.0);

    // this all has to be built here for unavoidable (I think) lifetime issues
    // (the FMIndex *requires* references to bwt,less,occ, which necessarily go out of scope)
    info!("Building suffix array");
    self.sa = suffix_array(&cat_seq); // slice the [u8] out of a Vec<u8>, I think
    let t_2:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_2 - t_1) as f64 / 1000000000.0);

    info!("Building BWT");
    self.bwt = bwt(&cat_seq, &self.sa);
    let t_3:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_3 - t_2) as f64 / 1000000000.0);

    info!("Computing less");
    self.less = less(&self.bwt, &alphabet);
    let t_4:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_4 - t_3) as f64 / 1000000000.0);

    info!("Computing occurrences");
    self.occ = Occ::new(&self.bwt, 3, &alphabet); // sampling rate of every 3rd position
    let t_5:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_5 - t_4) as f64 / 1000000000.0);

    info!("Building FM-index");
    let fm = FMIndex::new(&self.bwt, &self.less, &self.occ);
    self.fmd = FMDIndex::from(fm);
    let t_6:u64 = time::precise_time_ns();
    info!("  {:.4} seconds", (t_6 - t_5) as f64 / 1000000000.0);
  }
}

fn load_seq(ref_fa:&str, args:&Args, alphabet:&alphabets::Alphabet) -> (Vec<u8>, Vec<usize>, Vec<fasta::Record>) {
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

impl<'a> Overlapper for FMIFinder<'a> {
  fn ovlSeq(&self, seq:&[u8]) -> HashMap<u32,Vec<KmerMatch>> {
    let mut hit_hash:HashMap<u32,Vec<KmerMatch>> = HashMap::new();

    for i in 0..(seq.len()-self.args.k+1) {
      let kmer = &seq[i..i+self.args.k]; // just using string slice for the bwt

      let intervals = self.fmd.smems(kmer, 2);
      let forward_positions = intervals[0].forward().occ(&self.sa);
      let revcomp_positions = intervals[0].revcomp().occ(&self.sa);

      debug!("{:?}", forward_positions);
      debug!("{:?}", revcomp_positions);

      // convert to Vec<Position>
      let pos_vec:Vec<Position> = Vec::new();

      if pos_vec.len() < self.args.rep_limit {
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
}
