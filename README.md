Ralf is an experimental long-read sequence aligner in Rust.

Ralf depends heavily on the excellent work that is [rust-bio](https://github.com/rust-bio/rust-bio).

I maintain that Ralf is *not* an acronym, but you may remember it as Rust Aligner (Fast, Fun, FASTA)... if that helps.

Ralf is my very first Rust program, so I welcome any and all input.


Instal and Build
================

Rust must be installed to build from the source, see: https://www.rust-lang.org/

It has been tested on v1.7 Stable

    git clone https://github.com/txje/ralf.git
    cd ralf
    cargo build


Usage
=====

After build completes, you can run it as follow (on Unix systems):

    ./target/debug/rust_aln [options] <ref_fa> <read_fa>
  
Options:

    <ref_fa>                   reference fasta file (required)
    <read_fa>                  reads fasta file (required)
    -k <k>                     K-mer size to use for hash matching (default: 16).
    -m <min_ordered_matches>   Minimum number of properly ordered hits required to report an alignment (default: 10)
    -g <max_match_gap>         Maximum gap between adjacent matches (in query sequence) (default: 1000)
    -x <max_offset_variance>   Maximum difference between adjacent query and target match offsets, as (default: 0.2)
                               a fraction of total offset
    -o <min_allowable_offset>  Minimum difference to allow between adjacent query and target match offsets (default: 50)
    -l <min_aln_len>           Minimum reported alignment path (default: 100)

