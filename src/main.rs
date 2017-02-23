extern crate bio;
extern crate time;

// I guess this is required to enable info!, debug!, etc macros defined in log crate
#[macro_use] extern crate log;
mod easylog;

mod overlapper;
use overlapper::{Overlapper, Position, KmerMatch};
/*
mod fmifinder;
use fmifinder::{FMIFinder};
*/
mod seqhash;
use seqhash::SeqHash;
mod util;
use util::{Args,help_and_fail};

// for HashMap method
use std::collections::HashMap;

use std::vec::Vec;
use std::str::FromStr;
use std::env;

use bio::alphabets;
use bio::io::fasta;

fn main() {
  easylog::init().unwrap();
  
  // --------------------------------------------
  // Argument parser
  // --------------------------------------------
  let mut args:Args = Args{
    ref_fa:String::new(),
    read_fa:String::new(),
    k: 16,
    min_ordered_matches: 10, // at least 10 kmer hits
    max_match_gap: 1000, // maximum gap between kmer hits on query, bp
    max_offset_variance: 0.2, // 20% offset variance allowed
    max_abs_distance: 20, // 20 nt offset allowed before offset variance kicks in
    min_aln_len: 100, // minimum reported alignment path, bp
    rep_limit: 1000, // maximum occurrences of a single k-mer before it's considered repetitive and ignored
    bwt:false};
  let mut arg_iter = env::args();
  arg_iter.next(); // get rid of the command
  let n_args = arg_iter.len(); // only because it's an ExactSizeIterator
  let mut i:usize = 0;
  let mut positional:usize = 0;
  while i < n_args {
    let arg = arg_iter.next().unwrap_or(String::new());
    i += 1;
    if arg.len() == 0 {
      continue;
    }
    if arg.starts_with("--") {
      match arg.trim_left_matches("--").as_ref() {
        "help" => { help_and_fail(String::new()); },
        "min-ordered-matches" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<usize>() {
          Ok(n) => args.min_ordered_matches = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "max-match-gap" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<u32>() {
          Ok(n) => args.max_match_gap = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "max-offset-variance" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<f32>() {
          Ok(n) => args.max_offset_variance = n,
          Err(e) => help_and_fail(format!("'{}' requires a decimal value, got {}", arg, val)),
        }; },
        "max-abs-distance" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<u32>() {
          Ok(n) => args.max_abs_distance = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "min-aln-len" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<usize>() {
          Ok(n) => args.min_aln_len = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "rep-limit" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<usize>() {
          Ok(n) => args.rep_limit = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "bwt" => { args.bwt = true; },
        _ => { help_and_fail(format!("Unknown argument '{}'", arg)); }
      }
    } else if arg.starts_with("-") {
      match arg.trim_left_matches("-").as_ref() {
        "h" => { help_and_fail(String::new()); },
        "k" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<usize>() {
          Ok(n) => args.k = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "m" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<usize>() {
          Ok(n) => args.min_ordered_matches = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "g" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<u32>() {
          Ok(n) => args.max_match_gap = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "x" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<f32>() {
          Ok(n) => args.max_offset_variance = n,
          Err(e) => help_and_fail(format!("'{}' requires a decimal value, got {}", arg, val)),
        }; },
        "n" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<u32>() {
          Ok(n) => args.max_abs_distance = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "l" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<usize>() {
          Ok(n) => args.min_aln_len = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        "r" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<usize>() {
          Ok(n) => args.rep_limit = n,
          Err(e) => help_and_fail(format!("'{}' requires a positive integer value, got {}", arg, val)),
        }; },
        _ => { help_and_fail(format!("Unknown argument '{}'", arg)); }
      }
    } else {
      if positional == 0 {
        args.ref_fa = arg;
      } else if positional == 1 {
        args.read_fa = arg;
      } else {
        help_and_fail(format!("Too many positional arguments '{}'", arg));
      }
      positional += 1;
    }
  }
  if positional < 2 {
    help_and_fail("Not enough positional arguments".to_string());
  }

  //debug!("{:?}", args);


  let t_0:u64 = time::precise_time_ns();

  let alphabet = alphabets::dna::n_alphabet();

  let rev_shift:u8 = 2*(args.k as u8-1);
  let fw_mod:u64 = 1 << ((2*args.k) as u64);


  info!("Building SeqHash from {}", args.ref_fa);
  let finder = SeqHash::new(&args.ref_fa[..], &args, fw_mod, rev_shift, &alphabet);

  /*
  info!("Building BWT+FMD-Index from {}", args.ref_fa);
  let finder = FMIFinder::new(&args.ref_fa[..], &args, &alphabet);
  */


  let t_1:u64 = time::precise_time_ns();
  info!("  {:.4} seconds", (t_1 - t_0) as f64 / 1000000000.0);

  info!("Overlapping reads from {}", args.read_fa);
  ovl_reads(&args.read_fa[..], &finder, &args, fw_mod, &alphabet);

  let t_2:u64 = time::precise_time_ns();
  info!("  {:.4} seconds", (t_2 - t_1) as f64 / 1000000000.0);

  info!("Total time: {:.4} seconds", (t_2 - t_0) as f64 / 1000000000.0);
}


fn ovl_reads(read_fa:&str, finder:&Overlapper, args:&Args, fw_mod:u64, alphabet:&alphabets::Alphabet) {

  let reader = fasta::Reader::from_file(read_fa).unwrap();

  for record in reader.records() {

    let record = record.unwrap();
    let seq = record.seq();

    debug!("Read {} ({} bp)", record.id().unwrap(), seq.len());

    if !alphabet.is_word(seq) || seq.len() < args.k {
      continue;
    }

    let mut hit_hash = finder.ovlSeq(seq);

    for (rid, matches) in hit_hash.iter_mut() {

      //let rev = rid & 1 == 1;

      //matches.sort_by(|a, b| a.tpos.cmp(&b.tpos)); // sort by target position
      let mut matches = lis(&matches);
      matches = nonovl(&matches, &args);

      // don't bother if there aren't a minimum of matches to begin with
      if matches.len() < args.min_ordered_matches {
        continue;
      }

      let mut m_start = 0;

      let mut best_st = 0;
      let mut best_en = 0;
      for i in 1..matches.len() {
        // both must necessarily be sorted since query positions were already sorted and we got the LIS of target positions
        //let qdiff = if matches[i].qpos > matches[i-1].qpos {matches[i].qpos - matches[i-1].qpos} else {matches[i-1].qpos - matches[i].qpos};
        let qdiff = matches[i].qpos - matches[i-1].qpos;
        //let tdiff = if matches[i].tpos > matches[i-1].tpos {matches[i].tpos - matches[i-1].tpos} else {matches[i-1].tpos - matches[i].tpos};
        let tdiff = matches[i].tpos - matches[i-1].tpos;
        if (if qdiff > tdiff {qdiff - tdiff} else {tdiff - qdiff}) > args.max_abs_distance
          && (tdiff as f32 > qdiff as f32 * (1.0 + args.max_offset_variance) || qdiff as f32 > tdiff as f32 * (1.0 + args.max_offset_variance)) {
          if i-m_start > best_en-best_st+1 {
            best_st = m_start;
            best_en = i-1;
          }
          m_start = i;
        }
      }
      if matches.len()-m_start > best_en-best_st+1 {
        best_st = m_start;
        best_en = matches.len()-1;
      }


      if best_en-best_st+1 >= args.min_ordered_matches {
        debug!("  hit ref {} ({}) {} times at q{}-{}, t{}-{}", rid/2, rid&1, best_en-best_st+1, matches[best_st].qpos, matches[best_en].qpos, matches[best_st].tpos, matches[best_en].tpos);
      }
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


// greedily keep only nonoverlapping hits to the query sequence
fn nonovl(arr:&Vec<KmerMatch>, args:&Args) -> Vec<KmerMatch> {
  // matches are already ordered by qpos
  let mut s:Vec<KmerMatch> = Vec::new();
  let mut last:usize = 0;
  for i in 0..arr.len() {
     if i == 0 || arr[i].qpos >= arr[last].qpos + args.k as u32 {
       s.push(KmerMatch{qpos:arr[i].qpos, tpos:arr[i].tpos});
       last = i;
     }
  }
  s
}
