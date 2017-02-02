extern crate bio;
extern crate docopt;

// I guess this is required to enable info!, debug!, etc macros defined in log crate
#[macro_use] extern crate log;
mod easylog;

//use argparse::{ArgumentParser, StoreTrue, Store};
use docopt::Docopt;

use std::collections::HashMap;
use std::vec::Vec;
use std::fmt;
use std::str::FromStr;

use bio::alphabets;
use bio::io::fasta;


struct Position {
  rid: u32,
  pos: u32
}

struct KmerMatch {
  qpos: u32,
  tpos: u32
}

impl fmt::Display for KmerMatch {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{}:{}", self.qpos, self.tpos)
  }
}

struct Args {
  k: usize,
  min_ordered_matches: usize,
  max_match_gap: u32,
  max_offset_variance: f32,
  min_aln_len: usize,
}

// you have to be SUPER careful about this: see docopt.org
const USAGE: &'static str = "
Usage: rust_aln [options] <ref_fa> <read_fa>

ref_fa                     reference fasta file
read_fa                    reads fasta file
-k <k>                     K-mer size to use for hash matching (default: 16).
-m <min_ordered_matches>   Minimum number of properly ordered hits required to report an alignment
-g <max_match_gap>         Maximum gap between adjacent matches (in query sequence)
-x <max_offset_variance>   Maximum proportional difference between adjacent query and target match offsets
-l <min_aln_len>           Minimum reported alignment length
";


fn main() {
  easylog::init().unwrap();
  

  let mut args = Args {
    k: 16,
    min_ordered_matches: 10, // at least 10 kmer hits
    max_match_gap: 1000, // maximum gap between kmer hits on query, bp
    max_offset_variance: 0.2, // 20% offset variance allowed
    min_aln_len: 100, // minimum reported alignment path, bp
  };

  let arguments = Docopt::new(USAGE).and_then(|d| d.argv(std::env::args().into_iter()).parse()).unwrap_or_else(|e| e.exit());

  let ref_fa = arguments.get_str("<ref_fa>");
  let read_fa = arguments.get_str("<read_fa>");

  let argk = arguments.get_str("-k");
  args.k = match usize::from_str(argk) {
    Ok(argk) => argk,
    Err(_) => {
      warn!("Invalid value for -k, using default: 16");
      args.k
    }
  };

  let argm = arguments.get_str("-m");
  args.min_ordered_matches = match usize::from_str(argm) {
    Ok(argm) => argm,
    Err(_) => {
      warn!("Invalid value for -m, using default: 10");
      args.min_ordered_matches
    }
  };

  let argg = arguments.get_str("-g");
  args.max_match_gap = match u32::from_str(argg) {
    Ok(argg) => argg,
    Err(_) => {
      warn!("Invalid value for -g, using default: 1000");
      args.max_match_gap
    }
  };

  let mov = arguments.get_str("-x");
  args.max_offset_variance = match f32::from_str(mov) {
    Ok(mov) => mov,
    Err(_) => {
      warn!("Invalid value for -x, using default: 0.2");
      args.max_offset_variance
    }
  };

  let argl = arguments.get_str("-l");
  args.min_aln_len = match usize::from_str(argl) {
    Ok(argl) => argl,
    Err(_) => {
      warn!("Invalid value for -l, using default: 100");
      args.min_aln_len
    }
  };


  let alphabet = alphabets::dna::alphabet();

  let rev_shift:u8 = 2*(args.k as u8-1);
  let fw_mod:u64 = 1 << ((2*args.k) as u64);

  // ideally the initial capacity will reflect the total # of unique k-mers
  // in practice, 2x the total reference genome size accounts for exclusively
  //   unique k-mers on both strands
  let mut ref_hash:HashMap<u64,Position> = HashMap::with_capacity(400000000);

  info!("Building reference hash from {}", ref_fa);
  let ref_seqs = hash_ref(&ref_fa[..], &mut ref_hash, &args, fw_mod, rev_shift, &alphabet);

  info!("Overlapping reads from {}", read_fa);
  ovl_reads(&read_fa[..], ref_seqs, &mut ref_hash, &args, fw_mod, &alphabet);
}


fn ovl_reads(read_fa:&str, ref_seqs:Vec<fasta::Record>, ref_hash:&mut HashMap<u64,Position>, args:&Args, fw_mod:u64, alphabet:&alphabets::Alphabet) {

  let reader = fasta::Reader::from_file(read_fa).unwrap();

  for record in reader.records() {

    let record = record.unwrap();
    let seq = record.seq();

    if !alphabet.is_word(seq) || seq.len() < args.k {
      continue;
    }

    let mut hit_hash:HashMap<u32,Vec<KmerMatch>> = HashMap::new();

    let mut kmer:u64 = 0;
    for i in 0..args.k-1 {
      kmer = (kmer << 2) + ((seq[i] >> 1) & 3) as u64;
    }

    for i in 0..(seq.len()-args.k+1) {
      kmer = (kmer << 2) + ((seq[i+args.k-1] >> 1) & 3) as u64;
      if args.k < 32 {
        kmer = kmer % fw_mod;
      }

      match ref_hash.get(&kmer) {
        Some(pos) => {
          if !hit_hash.contains_key(&pos.rid) {
            let match_vec = Vec::new();
            hit_hash.insert(pos.rid, match_vec);
          }
          let mut match_vec = hit_hash.get_mut(&pos.rid).unwrap();
          match_vec.push(KmerMatch{qpos:i as u32, tpos: pos.pos});
        },
        None => {}
      };
    }

    for (rid, matches) in hit_hash.iter_mut() {

      let rev = rid & 1 == 1;

      // reverse match sets on rv strand
      if rev {
        matches.reverse();
      }

      // don't bother if there aren't a minimum of matches to begin with
      if matches.len() < args.min_ordered_matches {
        continue;
      }

      debug!("read {} ({} bp) hit ref {} ({}) {} times q{}:t{} -> q{}:t{}", record.id().unwrap(), seq.len(), rid/2, rid&1, matches.len(), matches[0].qpos, matches[0].tpos, matches[matches.len()-1].qpos, matches[matches.len()-1].tpos);

      let mut m_start = 0;

      let mut best_st = 0;
      let mut best_en = 0;
      for i in 1..matches.len() {
        let qdiff = if matches[i].qpos > matches[i-1].qpos {matches[i].qpos - matches[i-1].qpos} else {matches[i-1].qpos - matches[i].qpos};
        let tdiff = if matches[i].tpos > matches[i-1].tpos {matches[i].tpos - matches[i-1].tpos} else {matches[i-1].tpos - matches[i].tpos};
        if tdiff as f32 > qdiff as f32 * (1.0 + args.max_offset_variance)
          || qdiff as f32 > tdiff as f32 * (1.0 + args.max_offset_variance) {
          if i-m_start > best_en-best_st+1 {
            best_st = m_start;
            best_en = i-1;
          }
          m_start += 1;
        }
      }

      debug!("  best is {} hits from q{}-{}, t{}-{}", best_en-best_st+1, matches[best_st].qpos, matches[best_en].qpos, matches[best_st].tpos, matches[best_en].tpos);
/*
      for i in 1..matches.len() {

        // enforce maximum offset variance between adjacent matches
        let qdiff = if !rev {matches[i].qpos - matches[i-1].qpos} else {matches[i-1].qpos - matches[i].qpos};

        let tdiff = matches[i].tpos - matches[i-1].tpos;
        let (maxdiff,mindiff) = if qdiff > tdiff {(qdiff,tdiff)} else {(tdiff,qdiff)};
        let offset = maxdiff - mindiff;

        // if this gap is too big
        if maxdiff > args.max_match_gap || (offset as f32 > (mindiff as f32) * args.max_offset_variance) {
          // if there are enough preceding matches, output them
          if (if !rev {matches[i-1].qpos - matches[m_start].qpos + args.k as u32} else {matches[m_start].qpos - matches[i-1].qpos + args.k as u32}) >= args.min_aln_len as u32
            && i-m_start >= args.min_ordered_matches {
            report(&record, &ref_seqs[*rid as usize/2], &matches[m_start], &matches[i-1], i-m_start, rev, args);
          }
          m_start = i;
        }
      }

      // report last match set
      if (if !rev {matches[matches.len()-1].qpos - matches[m_start].qpos + args.k as u32} else {matches[m_start].qpos - matches[matches.len()-1].qpos + args.k as u32}) >= args.min_aln_len as u32
        && matches.len()-m_start >= args.min_ordered_matches {
        report(&record, &ref_seqs[*rid as usize/2], &matches[m_start], &matches[matches.len()-1], matches.len()-m_start, rev, args);
      }
*/
    }
  }
}

/*
 format:
 0   qName qSeqLength qStart qEnd qStrand
 5   tName tSeqLength tStart tEnd tStrand
 10  numKmerMatches
*/
fn report (query:&fasta::Record, target:&fasta::Record, first_match:&KmerMatch, last_match:&KmerMatch, nmatches:usize, rev:bool, args:&Args) {

  // it's probably NOT slow to have to get seq() only to find it's len
  // it should be returning a reference to a slice of the underlying string, so takes O(1) time
  let qlen = query.seq().len();

  print!("{} ", query.id().unwrap());
  print!("{} ", qlen);
  // reported position will be relative to the appropriate strand
  print!("{} ", if !rev {first_match.qpos} else {qlen as u32 - (first_match.qpos + args.k as u32)});
  print!("{} ", (if !rev {last_match.qpos + args.k as u32 - 1} else {qlen as u32 - 1 - last_match.qpos})); // right now, end position is INCLUSIVE
  print!("{} ", if !rev {0} else {1});

  print!("{} ", target.id().unwrap());
  print!("{} ", target.seq().len());
  print!("{} ", first_match.tpos);
  print!("{} ", last_match.tpos + args.k as u32 - 1);
  print!("0 "); // target will always show fw strand

  println!("{}", nmatches);
}


// longest increasing subsequence
// a la https://en.wikipedia.org/wiki/Longest_increasing_subsequence
fn lis(arr:&Vec<KmerMatch>) -> Vec<KmerMatch> {
  // matches are already ordered by qpos
  //let mut lis: Vec<KmerMatch> = Vec::new();
  let n = arr.len();

  // these need to act as fixed, variable-size arrays
  let mut p:Vec<usize> = Vec::with_capacity(n);
  let mut m:Vec<usize> = Vec::with_capacity(n + 1);
  // so we have to initialize them
  for _ in 0..n {
    p.push(0);
    m.push(0);
  }
  m.push(0);

  let mut l = 0;
  for i in 0..n {
    // Binary search for the largest positive j ≤ L
    // such that X[M[j]] < X[i]
    let mut lo = 1;
    let mut hi = l;
    while lo <= hi {
       let mid = if (lo+hi) & 1 == 1 {(lo+hi)/2 + 1} else {(lo+hi)/2};
       if arr[m[mid]].tpos < arr[i].tpos {
         lo = mid+1;
       } else {
         hi = mid-1;
       }
    }

    // After searching, lo is 1 greater than the
    // length of the longest prefix of X[i]

    // The predecessor of X[i] is the last index of 
    // the subsequence of length lo-1
    p[i] = m[lo-1];
    m[lo] = i;

    if lo > l {
      // If we found a subsequence longer than any we've
      // found yet, update L
      l = lo;
    }
  }

  // Reconstruct the longest increasing subsequence
  let mut s:Vec<KmerMatch> = Vec::new();
  let mut k = m[l];
  if l > 1 {
    for _ in 0..l {
      s.push(KmerMatch{qpos:arr[k].qpos, tpos:arr[k].tpos}); // poor man's clone
      k = p[k];
    }
  }

  s.reverse();
  s
}


fn hash_ref(ref_fa:&str, ref_hash:&mut HashMap<u64,Position>, args:&Args, fw_mod:u64, rev_shift:u8, alphabet:&alphabets::Alphabet) -> Vec<fasta::Record> {
  let mut rid:u32 = 0;
  let mut ref_seqs: Vec<fasta::Record> = Vec::new();

  // Iterate over a FASTA file, use the alphabet to validate read sequences
  let reader = fasta::Reader::from_file(ref_fa).unwrap();

  for record in reader.records() {
    // it's critical to do this in two steps, otherwise the Record from step 1 will go out of scope immediately
    // and the borrowed seq will crash the world.
    let record = record.unwrap();
    { // block to deal with hash seq, allowing the borrow of record to end before it's added to the ref_hash
      let seq = record.seq();

      //info!("Loading reference {}: {} ({}bp)", rid, record.id().unwrap(), seq.len());

      if alphabet.is_word(seq) {

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

          ref_hash.insert(kmer, Position{rid:rid<<1, pos:i as u32});
          ref_hash.insert(rev_kmer, Position{rid:(rid<<1)+1, pos:i as u32});
        }
      }
    }
    ref_seqs.push(record);
    rid += 1;
  }

  ref_seqs
}
