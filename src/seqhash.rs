/* -----------------------------------------------------------------
 * Simple HashMap-based k-mer lookup machine
 * -----------------------------------------------------------------*/

use std::fs;
use std::vec::Vec;
use std::collections::HashMap;

use bio::io::fasta;
use bio::alphabets;

use overlapper::{Overlapper, Position, KmerMatch};
use util::{Args,help_and_fail};

pub struct SeqHash<'a> {
  hash: HashMap<u64,Vec<Position>>,
  sequences: Vec<fasta::Record>,
  fw_mod: u64,
  args:&'a Args
}

impl<'a> SeqHash<'a> {
  pub fn new(ref_fa:&str, args:&'a Args, fw_mod:u64, rev_shift:u8, alphabet:&alphabets::Alphabet) -> SeqHash<'a> {
    let mut rid:u32 = 0;

    // ideally the initial capacity will reflect the total # of unique k-mers
    // in practice, 2x the total reference genome size accounts for exclusively
    //   unique k-mers on both strands

    let mut hash_size = 0; // this will never be used - if it can't infer size from the fasta file, the program will fail
    match fs::metadata(ref_fa) {
      Ok(md) => {hash_size = md.len() as usize*2;},
      Err(e) => help_and_fail(format!("Error reading file '{}'", ref_fa))
    };
    debug!("Creating SeqHash with capacity {}", hash_size);
    let mut sh:SeqHash = SeqHash{hash: HashMap::with_capacity(hash_size), sequences: Vec::new(), fw_mod:fw_mod, args:args};

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

        if alphabet.is_word(seq) {
          sh.kmer_pack_hash(&seq, &args, fw_mod, rev_shift, rid, slen);
        } else {
          info!("Reference sequence {} contains something that is not in (A,C,G,T,N). The *entire* sequence is being ignored.", record.id().unwrap());
        }
      }
      sh.sequences.push(record);
      rid += 1;
    }
    sh
  }

  fn kmer_pack_hash(&mut self, seq:&[u8], args:&Args, fw_mod:u64, rev_shift:u8, rid:u32, slen:usize) {
    let mut kmer:u64 = 0;
    let mut rev_kmer:u64 = 0;

    // set up sliding kmer (and reverse) as uint64
    for i in 0..args.k-1 {
      kmer = (kmer << 2) + ((seq[i] >> 1) & 3) as u64;
      rev_kmer = (rev_kmer >> 2) + ((((seq[i] as u64 / 2) & 3) ^ 2) << rev_shift);
    }

    for i in 0..(seq.len()-args.k+1) {
      kmer = (kmer << 2) + ((seq[i+args.k-1] >> 1) & 3) as u64;
      if args.k < 32 {
        kmer = kmer % fw_mod;
      }
      rev_kmer = (rev_kmer >> 2) + ((((seq[i+args.k-1] as u64 / 2) & 3) ^ 2) << rev_shift);

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
        // when aligned to the reverse target strand, positions are RELATIVE TO 5' (opposite) end strand
        pos_vec.push(Position{rid:(rid<<1)+1, pos:(slen-i-args.k) as u32});
      }
    }
  }
}

impl<'a> Overlapper for SeqHash<'a> {
  fn ovlSeq(&self, seq:&[u8]) -> HashMap<u32,Vec<KmerMatch>> {
    let mut hit_hash:HashMap<u32,Vec<KmerMatch>> = HashMap::new();

    let mut kmer:u64 = 0;
    for i in 0..self.args.k-1 {
      kmer = (kmer << 2) + ((seq[i] >> 1) & 3) as u64;
    }

    for i in 0..(seq.len()-self.args.k+1) {
      kmer = (kmer << 2) + ((seq[i+self.args.k-1] >> 1) & 3) as u64;
      if self.args.k < 32 {
        kmer = kmer % self.fw_mod;
      }

      match self.hash.get(&kmer) {
        Some(pos_vec) => {
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
        },
        None => {}
      };
    }

    hit_hash
  }
}
