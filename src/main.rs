extern crate bio;

extern crate time;

// I guess this is required to enable info!, debug!, etc macros defined in log crate
#[macro_use] extern crate log;
mod easylog;

mod overlapper;
use overlapper::{Overlapper, Position, KmerMatch};
mod fmifinder;
use fmifinder::{FMIFinder};
mod qgramfinder;
use qgramfinder::{QGramFinder};
mod seqhash;
use seqhash::SeqHash;
mod util;
use util::{Args, help_and_fail, bwt_alphabet};
mod dotplot;
use dotplot::{DotPlot, draw_dp};

use std::vec::Vec;
use std::env;
use std::panic;

use bio::alphabets;
use bio::alphabets::dna::revcomp;
use bio::io::fasta;

use bio::alignment::pairwise::Aligner;
use bio::alignment::{AlignmentOperation, Alignment};
use bio::alignment::AlignmentOperation::*; // enum: Match, Subst, Del, Ins

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
    dp_block_size: 1000, // bin size for dot plot
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
        "dp-block" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<usize>() {
          Ok(n) => args.dp_block_size = n,
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
        "d" => { let val = arg_iter.next().unwrap_or(String::new()); match val.parse::<usize>() {
          Ok(n) => args.dp_block_size = n,
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


  let t_0:u64 = time::precise_time_ns();

  let alphabet = bwt_alphabet();


  info!("Building SeqHash from {}", args.ref_fa);
  let finder = SeqHash::new(&args.ref_fa[..], args.k, args.rep_limit, &alphabet);

  /*
  // takes way too much memory unless you use my modified version that uses a reduced ACGTN$ alphabet and ranktransform
  // in which case it only uses a bunch of memory
  info!("Building BWT+FMD-Index from {}", args.ref_fa);
  let finder = FMIFinder::new(&args.ref_fa[..], args.k, args.rep_limit, &alphabet);
  */

  /*
  // requires some 10-20x (alphabet size)^k -- this is a lot and can only be practically used for k<~10
  info!("Building QGram index from {}", args.ref_fa);
  let finder = QGramFinder::new(&args.ref_fa[..], args.k, args.rep_limit, &alphabet);
  */


  let t_1:u64 = time::precise_time_ns();
  info!("  {:.4} seconds", (t_1 - t_0) as f64 / 1000000000.0);

  info!("Overlapping reads from {}", args.read_fa);
  ovl_reads(&args.read_fa[..], &finder, &args, &alphabet);

  let t_2:u64 = time::precise_time_ns();
  info!("  {:.4} seconds", (t_2 - t_1) as f64 / 1000000000.0);

  info!("Total time: {:.4} seconds", (t_2 - t_0) as f64 / 1000000000.0);
}


fn ovl_reads(read_fa:&str, finder:&Overlapper, args:&Args, alphabet:&alphabets::Alphabet) {

  let reader = fasta::Reader::from_file(read_fa).unwrap();

  // a 2d matrix for dots for every combination of query and target, hence 4D
  let mut dp = DotPlot{matrix: Vec::with_capacity(finder.sequences().len())};
  dp.matrix.resize(finder.sequences().len(), Vec::new()); // this actually fills it with 0s
  let dp_block_size = args.dp_block_size;

  let mut r = 0; // read counter
  for record in reader.records() {

    let record = record.unwrap();
    let seq = record.seq();

    // create each dp matrix
    for i in 0..dp.matrix.len() {
      dp.matrix[i].push(Vec::with_capacity(finder.sequences()[i].seq().len()/dp_block_size + 1));
      dp.matrix[i][r].resize(finder.sequences()[i].seq().len()/dp_block_size + 1, Vec::new()); // this actually fills it with 0s
      for j in 0..dp.matrix[i][r].len() {
        dp.matrix[i][r][j].resize(seq.len()/dp_block_size + 1, 0); // this actually fills it with 0s
      }
    }

    debug!("Read {} ({} bp)", record.id().unwrap(), seq.len());

    if !alphabet.is_word(seq) || seq.len() < args.k {
      continue;
    }

    let mut hit_hash = finder.ovlSeq(seq);

    for (rid, matches) in hit_hash.iter_mut() {

      // fill dotplot
      for i in 0..matches.len() {
        dp.matrix[(rid/2) as usize][r as usize][(if rid&1==0 {matches[i].tpos} else {finder.sequences()[(rid/2) as usize].seq().len() as u32 - (matches[i].tpos + args.k as u32)}) as usize / dp_block_size][matches[i].qpos as usize/dp_block_size] += 1;
      }

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
        let qdiff = matches[i].qpos - matches[i-1].qpos;
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
        let aln = align(&matches[best_st..best_en+1], &seq, rid, &finder.sequences()[(rid/2) as usize].seq(), args.k);
        report(&record, &finder.sequences()[(rid/2) as usize], &matches[best_st], &matches[best_en], best_en-best_st+1, rid&1==1, args.k, &aln);
      }
    }
    r += 1; // read counter
  }
  draw_dp(dp);
}

fn align(matches:&[KmerMatch], query_seq:&[u8], rid:&u32, target_seq:&[u8], k:usize) -> Alignment {
  let rev = rid&1 == 1;

  // this is a vector of the AlignmentOperation enum: Match, Subst, Ins, Del
  // initialize to a reasonably expected capacity for long read alignment:
  //   length of query string + 20%
  // this doesn't take into account extra alignment 5' and 3' of the terminal k-mer matches
  let max_query_len = matches[matches.len()-1].qpos - matches[0].qpos;
  let mut path:Vec<AlignmentOperation> = Vec::with_capacity((max_query_len as f32 * 1.2) as usize);
  let mut score:i32 = 0;

  // first k-mer match
  /*
  for _ in 0..k {
    path.push(Match);
  }
  score += k as i32;
  */

  let mut m_start:usize = 0;

  // compute DP alignment between matches
  for i in 1..matches.len() {
    let qdiff = matches[i].qpos - matches[i-1].qpos;
    let tdiff = matches[i].tpos - matches[i-1].tpos;
    // qdiff must be at least k apart - there is no such requirement for tdiff, so we'll align including the proximal k-mer match

    // lambda scoring function: 1 if match, -1 if mismatch
    let score_func = |a: u8, b: u8| if a == b {1i32} else {-1i32};
    let gap_open:i32 = -1;
    let gap_extend:i32 = -1;
    let mut aligner = Aligner::with_capacity(qdiff as usize, tdiff as usize, gap_open, gap_extend, &score_func);
    let alignment = if !rev {
      aligner.global(&query_seq[matches[i-1].qpos as usize..matches[i].qpos as usize], &target_seq[matches[i-1].tpos as usize..matches[i].tpos as usize])
    } else {
      aligner.global(&query_seq[matches[i-1].qpos as usize..matches[i].qpos as usize], &revcomp(target_seq)[matches[i-1].tpos as usize..matches[i].tpos as usize])
    };
    path.extend(alignment.operations);
    score += alignment.score;

    // k-mer match
    /*
    for _ in 0..k {
      path.push(Match);
    }
    score += k as i32;
    */
  }
  // last k-mer match
  for _ in 0..k {
    path.push(Match);
  }
  score += k as i32;

  let qalnlen = matches[matches.len()-1].qpos as usize - matches[0].qpos as usize + k;
  let talnlen = matches[matches.len()-1].tpos as usize - matches[0].tpos as usize + k;
  // ends are 1 PAST the last index
  Alignment{score:score, ystart:0, xstart:0, yend:talnlen, xend:qalnlen, xlen:qalnlen, operations:path}
}

/*
 format:
 0   qName qSeqLength qStart qEnd qStrand
 5   tName tSeqLength tStart tEnd tStrand
 10  numKmerMatches numMatches numSubst numDel numIns
 15  queryAln alnString targetAln
*/
fn report (query:&fasta::Record, target:&fasta::Record, first_match:&KmerMatch, last_match:&KmerMatch, nmatches:usize, rev:bool, k:usize, aln:&Alignment) {

  // it's probably NOT slow to have to get seq() only to find it's len
  // it should be returning a reference to a slice of the underlying string, so takes O(1) time
  let tlen = target.seq().len();

  print!("{} ", query.id().unwrap());
  print!("{} ", query.seq().len());
  print!("{} ", first_match.qpos);
  print!("{} ", last_match.qpos + k as u32 - 1);
  print!("0 "); // query is always the fw strand

  print!("{} ", target.id().unwrap());
  print!("{} ", tlen);
  // reported position will be relative to the appropriate strand
  print!("{} ", if !rev {first_match.tpos} else {tlen as u32 - (first_match.tpos + k as u32)});
  print!("{} ", (if !rev {last_match.tpos + k as u32 - 1} else {tlen as u32 - 1 - last_match.tpos})); // right now, end position is INCLUSIVE
  print!("{} ", if !rev {0} else {1});

  print!("{} ", nmatches);

  // compute simple stats
  let mut n_match:usize = 0;
  let mut n_mis:usize = 0;
  let mut n_del:usize = 0;
  let mut n_ins:usize = 0;
  for op in &aln.operations {
    match op {
      &Match => {n_match+=1;},
      &Subst => {n_mis+=1;},
      &Del => {n_del+=1;},
      &Ins => {n_ins+=1;}
    }
  }
  let acc = n_match as f32 / aln.operations.len() as f32;
  print!("{} {} {} {} {:.4} ", n_match, n_mis, n_del, n_ins, acc);

  let qseq = &query.seq()[first_match.qpos as usize .. last_match.qpos as usize + k];
  let tseq = &target.seq()[first_match.tpos as usize .. last_match.tpos as usize + k];
  let rv_tseq = &revcomp(tseq);
  let pretty_fmt = aln.pretty(qseq, if rev {rv_tseq} else {tseq}); // pretty() reports a string with 3 rows: query, alignment, target
  let mut parts = pretty_fmt.split("\n");
  print!("{} ", parts.next().unwrap());
  let aln_string = parts.next().unwrap().replace(" ", "*");
  print!("{} ", aln_string);
  print!("{} ", parts.next().unwrap());
  println!("");
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
