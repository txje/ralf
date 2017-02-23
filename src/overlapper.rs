extern crate bio;

use std::fmt;
use std::collections::HashMap;

use bio::io::fasta;

/*
 * Defines a trait returning a list of matching k-mers (Vec<Position>) from some structure
 * - ex. HashMap or BWT/FM-Index
 */

pub struct Position {
  pub rid: u32,
  pub pos: u32
}

pub struct KmerMatch {
  pub qpos: u32,
  pub tpos: u32
}

impl fmt::Display for KmerMatch {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{}:{}", self.qpos, self.tpos)
  }
}

impl fmt::Debug for KmerMatch {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{}:{}", self.qpos, self.tpos)
  }
}

pub trait Overlapper {
  fn ovlSeq(&self, &[u8]) -> HashMap<u32,Vec<KmerMatch>>;
  fn sequences(&self) -> &Vec<fasta::Record>;
}
