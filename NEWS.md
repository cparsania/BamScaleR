# BamScaleR 0.99.6

## Build

* Only C++ code now uses OpenMP. The parallel block compression of the
  vendored libBigWig writer moved from C (`bwWrite.c`) to C++
  (`src/bw_compress.cpp`), so `src/Makevars` and `src/Makevars.win` no longer
  put `$(SHLIB_OPENMP_CFLAGS)` in `PKG_CFLAGS`. This removes the
  `R CMD check --as-cran` NOTE "SHLIB_OPENMP_CFLAGS is included in PKG_CFLAGS
  but not in PKG_LIBS", and follows Writing R Extensions (section 1.2.1.1),
  which notes that using OpenMP from more than one of C, C++ and Fortran in a
  package is not portable. bigWig output is unchanged (byte-identical).
* `configure` (macOS) no longer probes C OpenMP flags, and
  `BAMSCALER_OPENMP_CFLAGS` is no longer read. OpenMP flags that a user's own
  configuration applies to C code (for example `CFLAGS += -Xclang -fopenmp`)
  no longer stop the installation.

# BamScaleR 0.99.5

* No code changes. Version bump so the Bioconductor submission bot
  reports the build of this version (its report labelled 0.99.4 showed
  the 0.99.3 build).

# BamScaleR 0.99.4

## Bug fixes

* A test of `bam_coverage_bigwig()` failed on Windows: it built an output
  path ending in a separator with `file.path(x, "")`, which drops the
  trailing separator on Windows. The package code is unchanged.

# BamScaleR 0.99.3

Changes in response to the Bioconductor review (BiocContributions#226).

## Bug fixes

* macOS: BamScaleR could fail to load with `symbol not found in flat namespace
  '___kmpc_critical'`. CRAN's R for macOS does not enable OpenMP itself; when
  it is enabled locally with `-Xclang -fopenmp` (for example in
  `~/.R/Makevars`), OpenMP code is generated but the OpenMP runtime (`libomp`)
  is not linked. A new `configure` script now tests candidate compile/link
  flags by building and `dyn.load()`ing a small library that uses the same
  OpenMP constructs as the package (including dynamic loop schedules), and
  falls back to a build without OpenMP (one thread per file, otherwise fully
  functional) when none works. With CRAN's R the runtime linked is always the
  `libomp` that ships with R; an OpenMP install from mac.r-project.org or
  Homebrew supplies the header. If R's runtime is too old for the compiler
  (Apple clang 17 needs LLVM 19 or later), the installation log says so. If
  the user's own flags switch OpenMP on for every package, so that a build
  without OpenMP is impossible, `configure` stops with a message saying so.
  Linux and Windows builds are unchanged. Flags can be set by hand with both
  `BAMSCALER_OPENMP_CXXFLAGS` and `BAMSCALER_OPENMP_LIBS`.
* Multi-file calls no longer modify the supplied `BPPARAM`. Previously the
  worker count could be reduced in place with `bpworkers<-`; because
  `BiocParallelParam` objects have reference semantics this changed the
  caller's object for later use and, for an already-started cluster, left its
  reported worker count out of step with the running workers. The files are
  passed to `bplapply()` directly, as before, except that a `SnowParam` that
  is not running is started with only the workers needed (no more than the
  files, or than an `auto_threads` plan allows) and stopped again afterwards;
  that under an `auto_threads` cap a `MulticoreParam` that is not running is
  replaced by a copy with fewer workers; and that under a cap the files of a
  running cluster, or of another kind of `BPPARAM`, are grouped into that many
  tasks (an error then names the file, and `bptry()` returns per-group
  results).
* Linux: once OpenMP threads had been used in the R session (by BamScaleR or
  another package, such as data.table), a multi-file call with
  `MulticoreParam` and `threads > 1` could hang forever, because GCC's OpenMP
  runtime cannot start threads in a forked process once its parent has used
  them. BamScaleR now releases the runtime's idle threads before each `fork()`
  (with libgomp 9 or later, looked up at run time), so forked workers use the
  requested threads. This applies to every `fork()` in the session and may
  discard other libraries' OpenMP thread-private data. With an older libgomp,
  a worker forked after BamScaleR has used threads reads on one thread, and
  one forked after another package used OpenMP threads can still hang; use
  `SnowParam(type = "SOCK")` there.
* `param$which` kept only records whose start position lay in a range, so
  alignments that start before a range and extend into it were missing. It now
  keeps every alignment that overlaps a range, as `Rsamtools` does (each
  record once, under the first range it overlaps). As in `Rsamtools`, a
  zero-width range now selects nothing and a range on a sequence that is not
  in the BAM header is an error (it silently matched nothing). An empty
  `which` now means no restriction instead of an error; ranges with missing
  values, or with an end less than start - 1, are an error; repeated range
  labels are made unique (so a range given twice gets an empty second batch).
* `bam_coverage_bigwig()`:
  - compressed blocks on every available core whatever `threads` was; it now
    uses `threads`;
  - a failure while finishing the file (for example a full disk) left a
    truncated bigWig and returned normally; it is now an error and the partial
    file is removed;
  - opened the output in text mode on Windows (the vendored libBigWig ignored
    the binary mode it was given), which would corrupt the file;
  - with `parallel = FALSE`, a write failure while writing the data (for
    example a full disk) overran the vendored libBigWig's write buffer, which
    could crash R; it is now an error and the partial file is removed;
  - crashed R for an `outfile` in a directory that does not exist or cannot be
    written. That, an `outfile` that is a directory (or ends in a path
    separator) or cannot be overwritten,
    repeated `outfile` paths (also through symbolic links) and an `outfile`
    that is one of the input BAMs are now errors before any work. Output paths
    are made absolute (and `~` expanded) before they reach the workers.
* `bam_read(as = "GAlignmentPairs")` paired the i-th first mate with the i-th
  last mate of each read name, which mismatched mates when a read had
  secondary or supplementary alignments. Mates are now matched by
  `GenomicAlignments::makeGAlignmentPairs()` (read name, mate position,
  strand, first/last, proper-pair and secondary flags); without `which` this gives the
  same pairs as `readGAlignmentPairs()`, ordered by first-mate position. When
  nothing can be paired, an empty `GAlignmentPairs` is returned instead of an
  error.
* `mapq_dist()` no longer counts unmapped reads (at MAPQ 0), and
  `fragment_sizes()` with `drop_mate_unmapped = TRUE` now also drops unmapped
  reads, as `scanBam()` reports `NA` for both.
* Settings of `ScanBamParam` that are not supported (`tagFilter`,
  `simpleCigar`, `reverseComplement`, and `which` in the summary functions),
  and unknown entries of a list-form `param`, now give a warning instead of
  being ignored silently. A list-form `param` now applies `flag`
  (a `scanBamFlag()` value, which was silently dropped) and matches entry
  names exactly.
* `use.names = TRUE` now names `GAlignments` output even when `what` does
  not include `"qname"`, and `with.which_label = TRUE` adds a column only
  when `which` is set. Input paths starting with `~` are now expanded, and
  relative input paths (also of `BamFile` and `BamFileList` inputs) are made
  absolute, so `SnowParam` workers started in another directory read the
  intended files.

## Improvements

* The `GAlignments` fast path no longer calls the non-exported
  `S4Vectors:::new2()`. The object is assembled by slot assignment on a valid
  prototype, and the C++ reader counts the records whose validity it cannot
  guarantee by construction (a CIGAR with an op code above 8 or a zero-length
  op, and an alignment outside its sequence); if there are any,
  `validObject()` is run, so an invalid object is never returned and
  out-of-bound alignments warn as with `readGAlignments()`.
* Requesting `threads` > 1 on a build without OpenMP now warns once per
  session instead of silently running single-threaded, and `auto_threads`
  then keeps all `BPPARAM` workers.

## Documentation

* The main vignette has a new section with evaluated examples of
  `fragment_sizes()`, `mapq_dist()`, `bam_coverage()` and
  `bam_coverage_bigwig()`, each compared with the standard Bioconductor
  route; an installation section covering OpenMP on macOS; and a section on
  what BamScaleR adds to ompBAM.
* The benchmark vignette now computes every table and number from CSV files
  shipped in `inst/extdata/` (canonical run of 2026-10-08), with their
  provenance documented in `inst/scripts/make_benchmark_extdata.R` and a table of
  the ENCODE accessions used. Two thread-scaling statements that had been
  carried over from an earlier, smaller run were corrected.
* Expanded the DESCRIPTION `Description` field.

## Dependencies

* `parallel` is now declared in `Imports`; `testthat` (>= 3.1.7) and `withr`
  (used by the test setup) in `Suggests`.
* Removed unused `Suggests`: `chipseqDBData`, `dplyr`, `ExperimentHub`,
  `readr`, `tibble`, `tidyr`. Packages needed only by the benchmark scripts in
  `inst/benchmarks/` (including `chipseqDBData` and `ATACseqQC`) are listed in
  `inst/benchmarks/README.md` and are intentionally not declared in
  DESCRIPTION.

# BamScaleR 0.99.2

## Documentation

* README and the benchmark vignette now report the canonical benchmark run of
  2026-10-08, measured on the full 25-file public ENCODE ATAC-seq corpus
  (single-file target ENCFF191LNL, 151,057,932 alignments). The previous tables
  came from a much smaller early run and overstated some results: the
  `GAlignments` read speedup was shown as 4.03x against a measured 2.50x, and
  coverage to an in-memory `RleList` as 3.29x against 2.49x. Read-pattern
  speedups are 1.94-2.50x, ATAC fragment-size QC 4.39x end-to-end, and coverage
  to bigWig 1.52x.
* Added the measured in-reader aggregation results to the README
  (5.47-7.06x over the standard Bioconductor routes; 3.45x over megadepth at
  matched core counts) and the memory characteristic of `bam_coverage()`.
* Documented the rename from BamScale and the distinction from the unrelated
  NCBI command-line tool BAMscale.
* Installation instructions no longer advertise a Bioconductor install that
  does not yet work; BamScaleR is still under review.
* Recorded two further limitations: sublinear intra-file thread scaling, and
  that object-free command-line decoders reach higher raw throughput.

# BamScaleR 0.99.1

## Bug fixes

* `bam_coverage()` and `bam_coverage_bigwig()` no longer allocate memory
  proportional to the genome declared in the BAM header. The coverage core
  previously built a dense per-base difference array for every header sequence
  before reading a single alignment, so a BAM carrying few alignments still cost
  several GB whenever its header declared a whole genome -- 12.6 GB for a 1.4 MB
  file with a GRCh38 header -- and the collapse pass then walked every base of
  every contig. Coverage is now accumulated as per-thread interval endpoints and
  swept per contig, so memory and time scale with the alignments actually
  present rather than with the reference. Contigs without alignments cost
  nothing and still emit the correct full-length zero run.

  Measured on a 2.5 GB ATAC-seq BAM (195 contigs): peak RSS 14.2 GB -> 2.5 GB
  and 71.0 s -> 14.9 s at 24 threads. On a 1.4 MB BAM with a whole-genome
  header: 12.6 GB -> 0.55 GB (R's baseline) and 17.6 s -> 7.8 s. Output is
  unchanged and remains `identical()` to `coverage(readGAlignments())` at every
  thread count.

  Removing the shared difference array also removes the atomic updates from the
  read loop, which is where most of the speedup comes from.

# BamScaleR 0.99.0

## Package rename

* The package has been renamed from **BamScale** to **BamScaleR**. The former
  name collided with the existing NCBI command-line tool BAMscale
  (<https://github.com/ncbi/BAMscale>; Pongor et al., *Epigenetics & Chromatin*
  2020, doi:10.1186/s13072-020-00343-x), which performs a different task:
  normalised ("scaled") coverage-track generation and peak quantification. The
  rename removes that ambiguity.
* No user-facing functionality changed. All exported functions keep their names,
  arguments and behaviour.
* The option `BamScale.ga_fastpath` is now `BamScaleR.ga_fastpath`.
* BamScale was accepted into Bioconductor devel on 4 May 2026 and was never part
  of a numbered Bioconductor release, so no deprecation cycle applies. BamScaleR
  is submitted as its renamed continuation and version numbering restarts at
  0.99.0, as required for a new submission. The entries below document this
  package's history under its former name.

## Earlier history (as BamScale)

The entries below document this package's development under its former name,
**BamScale**, before the rename described above. Those version numbers belong to
the former name's series: BamScale was accepted into Bioconductor devel on
4 May 2026 and never shipped in a numbered release, so version numbering
restarts at 0.99.0 under the new name. The code described below is carried
forward into BamScaleR unchanged.

### BamScale 0.99.14

#### New features

* Four new exported functions compute common BAM summaries entirely inside the
  multithreaded C++ reader (per-thread accumulators folded in the OpenMP region;
  only the compact result crosses into R -- no per-read R objects are ever
  materialised), each verified byte-identical to its standard Bioconductor
  equivalent:
  - `fragment_sizes()`: paired-end fragment-size (insert-size) distribution;
    identical to `table(abs(scanBam(isize)))`. `drop_mate_unmapped = TRUE`
    (default) replicates `scanBam`'s `isize = NA` rule for mate-unmapped reads.
  - `mapq_dist()`: mapping-quality distribution; identical to
    `table(scanBam(mapq))` (MAPQ 255 reported as 255, where scanBam uses `NA`).
  - `bam_coverage()`: per-base coverage as a `SimpleRleList`; `identical()` to
    `GenomicAlignments::coverage(readGAlignments(...))`, including the CIGAR
    conventions (`M/=/X/D` cover, `N` gaps, `drop.D.ranges = FALSE`) and the
    readGAlignments record filter (only unmapped reads dropped).
  - `bam_coverage_bigwig()`: single-pass BAM -> coverage bigWig written natively
    via the bundled libBigWig -- no R-side coverage object, no
    `rtracklayer::export.bw` round-trip. Data and zoom blocks are compressed
    across the OpenMP threads (`parallel = TRUE`, default); the parallel output
    is bit-for-bit identical to the serial writer. A `compress_level` argument
    exposes the zlib deflate level, and `verbose = TRUE` prints a per-phase
    timing breakdown.
* `bam_read(as = "GAlignments")` now uses a dedicated fast path that assembles
  the `GAlignments` slots in C++: seqnames/strand accumulated as run-length
  pairs inside the parallel region, the cigar column materialised through a
  CHARSXP cache that exploits CIGAR-string redundancy, and qname decoded only
  when actually requested. Output is `identical()` to
  `GenomicAlignments::readGAlignments()` (verified at 226M-read scale) and to
  the previous construction path. `options(BamScale.ga_fastpath = FALSE)`
  restores the old path.
* Vendored libBigWig 0.4.8 (MIT, Devon Ryan) under `src/libBigWig/` (built with
  `-DNOCURL`, local files only). Local modifications, each commented in
  `bwWrite.c`: a tunable zlib deflate level, deferred parallel block
  compression, and zoom-level construction from in-memory coverage runs instead
  of re-reading the just-written data section.

#### Performance

* Fixed a quadratic reallocation pattern in the reader's batch merge: per-batch
  exact-capacity `reserve()` calls copied all accumulated data on every batch.
  Buffers now grow geometrically. This affects every read path -- e.g. on a
  229M-read BAM at 48 threads, `as = "GAlignments"` dropped from ~300 s to
  ~56 s, and plain data-frame reads improved 2.8-4.4x.
* CIGAR strings are formatted with an allocation-free digit writer instead of
  per-op `std::to_string`.
* Representative large-file results (229M-read ATAC BAM, 48 threads, warm
  cache): `bam_coverage()` 7.2x vs `coverage(readGAlignments())`;
  `bam_coverage_bigwig()` 2.6x vs megadepth at equal compression;
  `fragment_sizes()` 5.4x vs `scanBam` + `table()`; `as = "GAlignments"` 4.8x
  vs `readGAlignments()` (and faster even single-threaded).

#### Bug fixes

* `bam_coverage_bigwig()` now works on Windows: the output stream is opened in
  binary mode (`"w+b"`; text mode corrupted the bigWig), and zoom levels are
  always built from the in-memory coverage runs instead of re-reading the
  just-written stream (Windows' C runtime is strict about write-to-read
  transitions on update streams; the in-memory path is byte-identical to the
  re-read path, verified at 226M-read scale). Vendored libBigWig diagnostics now
  go through `REprintf()`, clearing the compiled-code check NOTE.
* Removed a long-standing `^src/Makevars$` entry from `.Rbuildignore` (dating to
  the pkgdown setup): built tarballs -- including the Bioconductor builds --
  shipped without `src/Makevars` and therefore compiled WITHOUT OpenMP (and
  without the explicit `-lz`), silently running single-threaded when installed
  from a tarball. In-place `R CMD INSTALL .` was unaffected, which is why local
  installs always threaded. Tarball builds now carry the full build
  configuration (OpenMP, zlib, `-DNOCURL`, and the bundled libBigWig objects).

#### Benchmarks

* `run_workflow_benchmark.R` gains first-class arms for the four new functions
  (`--include-fastcov/-bigwig/-fragsize/-mapq`), each behind an `identical()`
  correctness gate, with an optional megadepth context arm (`--megadepth-bin`).
* `run_server_benchmark.R`: BiocParallel PSOCK clusters are now started and
  their namespaces pre-loaded outside the timed region in all multi-file arms
  (previously cluster spin-up was charged to every iteration, diluting
  comparator ratios at high worker counts); added `--include-single` so single
  and multi scenarios can run in separate stages.
* New reproducible manuscript-benchmark tooling: `select_encode_atac.R`
  (ENCODE portal selection -> accession/md5 manifest), manifest-driven
  `download_atac_data.R` with md5 verification, `samtools_reference.sh`,
  `merge_runs.R`, and the staged `run_manuscript_final.sh` orchestrator.

#### Tests

* New testthat gates: byte-identity of `fragment_sizes()` / `mapq_dist()` /
  `bam_coverage()` against their Rsamtools/GenomicAlignments equivalents;
  `bam_coverage_bigwig()` value-identity after re-import plus bit-for-bit
  parallel-vs-serial output; six `GAlignments` fast-path identity gates
  (threads 1/4, mcols order, `use.names`, flag/mapq filters, fast == slow,
  `gctorture`).

### BamScale 0.99.13

#### New features

* On the compatible read path, `seq` and `qual` are now built directly in C as
  `DNAStringSet` / `PhredQuality` objects from the packed BAM byte buffer (via the
  IRanges/XVector C-callable constructors), avoiding a character round-trip. Output
  remains byte-identical to `Rsamtools::scanBam`.
* Added an end-to-end workflow benchmark harness under `inst/benchmarks/`
  (coverage → bigWig and ATAC fragment-size QC) and a self-contained
  `benchmark-results` vignette.

#### Bug fixes

* `GAlignments` output now carries the BAM-header sequence lengths in its `Seqinfo`.
  Previously the seqlengths were `NA`, so `coverage()` extended each seqlevel only to
  its maximum end and `export.bw()` wrote truncated chromosome sizes, diverging from
  `GenomicAlignments::readGAlignments`.
* Force-load the C-callable providers (`S4Vectors`, `IRanges`, `XVector`, `Biostrings`)
  and `GenomicAlignments` in `.onLoad`, so the compatible seq/qual and `GAlignments`
  paths resolve their constructors and classes on fresh `SnowParam` SOCK workers.

#### Documentation

* Refreshed the README with the latest benchmark results (read-pattern and end-to-end
  workflow speedups) and Bioconductor installation instructions.
