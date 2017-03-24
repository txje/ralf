/* -----------------------------------------------------------------
 * Simple HashMap-based k-mer lookup machine
 * -----------------------------------------------------------------*/

use std::fs;
use std::vec::Vec;
use std::collections::HashMap;

use bio::io::fasta;
use bio::alphabets;

use overlapper::{Overlapper, Position, KmerMatch};
use util::{help_and_fail};

pub struct SeqHash {
  hash: HashMap<u64,Vec<Position>>,
  sequences: Vec<fasta::Record>,
  k:usize,
  rep_limit:usize
}

impl SeqHash {
  pub fn new(ref_fa:&str, k:usize, rep_limit:usize, alphabet:&alphabets::Alphabet) -> SeqHash {
    let mut rid:u32 = 0;

    // ideally the initial capacity will reflect the total # of unique k-mers
    // in practice, 2x the total reference genome size accounts for exclusively
    //   unique k-mers on both strands

    let mut hash_size = 0; // zero will never be used - if it can't infer size from the fasta file, the program will fail
    match fs::metadata(ref_fa) {
      Ok(md) => {hash_size = md.len() as usize*2;},
      Err(e) => help_and_fail(format!("Error reading file '{}'", ref_fa))
    };
    debug!("Creating SeqHash with capacity {}", hash_size);
    let mut sh:SeqHash = SeqHash{hash: HashMap::with_capacity(hash_size), sequences: Vec::new(), k:k, rep_limit:rep_limit};

    // Iterate over a FASTA file, use the alphabet to validate read sequences
    let reader = fasta::Reader::from_file(ref_fa).unwrap();

    for record in reader.records() {
      // it's critical to do this in two steps, otherwise the Record from step 1 will go out of scope immediately
      // and the borrowed seq will crash the world.
      let record = record.unwrap();
      { // block to deal with hash seq, allowing the borrow of record to end before it's added to the ref_hash
        let seq = record.seq();
        let slen = seq.len();

        //debug!("Loading reference {}: {} ({}bp)", rid, record.id().unwrap(), seq.len());

        if !alphabet.is_word(seq) {
          info!("Reference sequence {} contains something that is not in (A,C,G,T,N). The *entire* sequence is being ignored.", record.id().unwrap());
        } else if slen < k {
          info!("Reference sequence {} is <{} bp, it cannot be hashed and will be ignored.", record.id().unwrap(), k);
        } else {
          sh.kmer_pack_hash(&seq, rid, slen);
        }
      }
      sh.sequences.push(record);
      rid += 1;
    }
    sh
  }

  fn kmer_pack_hash(&mut self, seq:&[u8], rid:u32, slen:usize) {
    let mut kmer:u64 = 0;
    let mut rev_kmer:u64 = 0;
    let fw_mod:u64 = 1 << ((2*self.k) as u64);
    let rev_shift:u8 = 2*(self.k as u8-1);

    // set up sliding kmer (and reverse) as uint64
    for i in 0..self.k-1 {
      kmer = (kmer << 2) + ((seq[i] >> 1) & 3) as u64;
      rev_kmer = (rev_kmer >> 2) + ((((seq[i] as u64 / 2) & 3) ^ 2) << rev_shift);
    }

    for i in 0..(seq.len()-self.k+1) {
      kmer = (kmer << 2) + ((seq[i+self.k-1] >> 1) & 3) as u64;
      if self.k < 32 {
        kmer = kmer % fw_mod;
      }
      rev_kmer = (rev_kmer >> 2) + ((((seq[i+self.k-1] as u64 / 2) & 3) ^ 2) << rev_shift);

      if !self.hash.contains_key(&kmer) {
        let pos_vec:Vec<Position> = Vec::new();
        self.hash.insert(kmer, pos_vec);
      }
      {
        let mut pos_vec = self.hash.get_mut(&kmer).unwrap();
        pos_vec.push(Position{rid:rid<<1, pos:i as u32});
      }

      if !self.hash.contains_key(&rev_kmer) {
        let pos_vec:Vec<Position> = Vec::new();
        self.hash.insert(rev_kmer, pos_vec);
      }
      {
        let mut pos_vec = self.hash.get_mut(&rev_kmer).unwrap();
        // when aligned to the reverse target strand, positions are RELATIVE TO 3' (opposite) end strand
        pos_vec.push(Position{rid:(rid<<1)+1, pos:(slen-i-self.k) as u32});
      }
    }
  }
}

impl Overlapper for SeqHash {
  fn ovlSeq(&self, seq:&[u8]) -> HashMap<u32,Vec<KmerMatch>> {
    let fw_mod:u64 = 1 << ((2*self.k) as u64);
    let mut hit_hash:HashMap<u32,Vec<KmerMatch>> = HashMap::new();

    let mut kmer:u64 = 0;
    for i in 0..self.k-1 {
      kmer = (kmer << 2) + ((seq[i] >> 1) & 3) as u64;
    }

    for i in 0..(seq.len()-self.k+1) {
      kmer = (kmer << 2) + ((seq[i+self.k-1] >> 1) & 3) as u64;
      if self.k < 32 {
        kmer = kmer % fw_mod;
      }

      match self.hash.get(&kmer) {
        Some(pos_vec) => {
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
        },
        None => {}
      };
    }

    hit_hash
  }

  fn sequences(&self) -> &Vec<fasta::Record> {
    &self.sequences
  }
}
