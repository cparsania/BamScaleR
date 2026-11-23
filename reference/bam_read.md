# Fast BAM reading with Bioconductor-compatible arguments

`bam_read()` is a multithreaded sequential BAM reader built on top of
`ompBAM`. The interface is designed to be familiar to users of
[`Rsamtools::scanBam()`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html),
[`GenomicAlignments::readGAlignments()`](https://rdrr.io/pkg/GenomicAlignments/man/readGAlignments.html),
and
[`GenomicAlignments::readGAlignmentPairs()`](https://rdrr.io/pkg/GenomicAlignments/man/readGAlignments.html).

## Usage

``` r
bam_read(
  file,
  param = NULL,
  what = NULL,
  tag = NULL,
  as = c("DataFrame", "data.frame", "GAlignments", "GAlignmentPairs", "scanBam"),
  seqqual_mode = c("compatible", "compact"),
  threads = 1L,
  BPPARAM = BiocParallel::bpparam(),
  auto_threads = FALSE,
  use.names = FALSE,
  with.which_label = FALSE,
  include_unmapped = TRUE
)
```

## Arguments

- file:

  A BAM input. Supported values are:

  - a single BAM path (`character(1)`) or multiple BAM paths,

  - a
    [`Rsamtools::BamFile`](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html),

  - a
    [`Rsamtools::BamFileList`](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html).

- param:

  Optional
  [`Rsamtools::ScanBamParam`](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html)
  (or a compatible list for lightweight use). The following fields are
  honored: `mapqFilter`, `flag`, `which`, `what`, and `tag`.
  `tagFilter`, `simpleCigar` and `reverseComplement` are not supported;
  setting them gives a warning and they are ignored.

- what:

  Character vector of fields to return, similar to `scanBam(what=...)`.
  Supported fields are `qname`, `flag`, `rname`, `strand`, `pos`,
  `qwidth`, `mapq`, `cigar`, `mrnm`, `mpos`, `isize`, `seq`, `qual`.

- tag:

  Character vector of 2-letter tag names to extract.

- as:

  Output format:

  - `"DataFrame"`: returns
    [`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
    (default),

  - `"data.frame"`: returns base `data.frame`,

  - `"GAlignments"`: returns
    [`GenomicAlignments::GAlignments`](https://rdrr.io/pkg/GenomicAlignments/man/GAlignments-class.html),

  - `"GAlignmentPairs"`: returns
    [`GenomicAlignments::GAlignmentPairs`](https://rdrr.io/pkg/GenomicAlignments/man/GAlignmentPairs-class.html),

  - `"scanBam"`: returns a `scanBam()`-shaped list-of-lists.

- seqqual_mode:

  Controls representation of `seq`/`qual` when those fields are
  requested:

  - `"compatible"` (default): return `seq` and `qual` as
    [`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
    and
    [`Biostrings::PhredQuality`](https://rdrr.io/pkg/Biostrings/man/XStringQuality-class.html)
    objects, identical to those from
    [`Rsamtools::scanBam()`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html).
    These are best held in a `DataFrame` (`as = "DataFrame"`); a base
    `data.frame` can carry them but cannot print them,

  - `"compact"`: return lower-level raw list-columns for
    faster/lower-overhead extraction. In compact mode, `seq` is returned
    as one `raw` vector per read containing BAM-native packed sequence
    bytes (two bases per byte), and `qual` is returned as one `raw`
    vector per read containing per-base Phred bytes. These are not plain
    character strings. `qwidth` is needed to decode compact `seq` back
    to base letters, and `qual` values of `255` correspond to missing
    quality values. This mode is currently supported for
    `as = "data.frame"` or `as = "DataFrame"`.

- threads:

  Requested number of OpenMP threads used for reading/decompression. May
  be capped when `auto_threads = TRUE`. Values above 1 require a build
  with OpenMP (see the Installation section of
  [`vignette("BamScaleR-intro")`](https://cparsania.github.io/BamScaleR/articles/BamScaleR-intro.md));
  without it a warning is issued once per session and each file is read
  on one thread.

- BPPARAM:

  `BiocParallel` parameter used when `file` contains more than one BAM.
  Defaults to
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).
  Set to `NULL` to force serial file processing. The object is used as
  supplied and is never modified, so it can be reused afterwards; an
  already-started cluster stays running with its configured number of
  workers.

- auto_threads:

  Logical; when `TRUE` and `BPPARAM` has multiple workers, BamScaleR
  adaptively avoids oversubscription by preserving higher per-file
  OpenMP thread counts when possible and reducing the number of
  concurrently active file workers before shrinking per-file threads.

- use.names:

  Passed to alignment object conversion. When `TRUE`, read names
  (`qname`) are used as object names (they are read for this even when
  `what` does not include `"qname"`).

- with.which_label:

  Logical; if `TRUE` and `param` includes `which`, an extra
  `which_label` column is returned. It is a character column (or mcol),
  where `readGAlignments()` returns a factor `Rle`.

- include_unmapped:

  Logical; whether unmapped records are retained (subject to
  `param$flag` constraints).

## Value

If `file` is length 1: one object in the format specified by `as`. If
`file` has length \> 1 (or is a `BamFileList`): a named list of outputs,
one per BAM file.

## Details

`bam_read()` is intentionally column-compatible with common BAM fields
used by Bioconductor workflows and can be used in place of
[`Rsamtools::scanBam()`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html)
for the supported fields and settings; the differences are listed below.

Parallelism model:

- `BPPARAM` parallelizes across files (one file per BiocParallel
  worker).

- `threads` parallelizes within each file via OpenMP.

- Effective total concurrency is approximately
  `min(length(file), BiocParallel::bpnworkers(BPPARAM)) * threads`.

- If `auto_threads = TRUE` and `BPPARAM` has multiple workers, BamScaleR
  first limits the number of concurrently active workers to preserve the
  requested per-file thread count within the detected core budget, then
  caps per-file OpenMP threads only if a single file would still
  oversubscribe the machine. `BPPARAM` itself is never changed. A
  `MulticoreParam` that is not running is replaced by a copy with fewer
  workers; a `SnowParam` that is not running is started with only the
  workers needed (also when there are fewer files than workers) and
  stopped again afterwards. For a running cluster, or another kind of
  `BPPARAM`, the files are grouped into that many `bplapply()` tasks; an
  error then stops its whole group, and
  [`BiocParallel::bptry()`](https://rdrr.io/pkg/BiocParallel/man/bptry.html)
  returns one result per group rather than per file.

- GCC's OpenMP runtime (the usual one on Linux) cannot start threads in
  a process created by `fork()` once its parent has used OpenMP threads,
  and such a forked worker
  ([`BiocParallel::MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html))
  would wait forever. BamScaleR therefore releases the runtime's idle
  threads just before each `fork()` in the session, so forked workers
  can use `threads` as requested, whichever package used OpenMP first.
  (This may also discard other libraries' OpenMP thread-private data.)
  With an OpenMP runtime too old to support this (GCC's libgomp before
  version 9), a worker forked after BamScaleR has used several threads
  reads on one thread, and a worker forked after another package used
  OpenMP threads can still hang; use
  [`BiocParallel::SnowParam()`](https://rdrr.io/pkg/BiocParallel/man/SnowParam-class.html)
  (type `"SOCK"`) there.

Compatibility notes:

- `rname`, `mrnm`, and `strand` are returned as `factor` columns,
  matching
  [`Rsamtools::scanBam()`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html).
  `rname`/`mrnm` levels are the BAM header reference names; unmapped or
  unset references are `NA` (not `"*"`). `strand` has levels
  `c("+", "-", "*")`.

- Region filtering via `param$which` is a sequential filter (the whole
  file is read; the BAM index is not used). As in `Rsamtools`, a record
  is kept when its alignment overlaps a range. Unlike `Rsamtools`, a
  record that overlaps several ranges is returned once, labelled with
  the first of them; records come in file order rather than grouped by
  range; and unmapped records are not returned when `which` is set. As
  in `Rsamtools`, a zero-width range selects nothing, and a range on a
  sequence that is not in the BAM header is an error. Repeated range
  labels are made unique with
  [`make.unique()`](https://rdrr.io/r/base/make.unique.html), one batch
  per range; a range given twice therefore gets an empty second batch,
  where `scanBam()` repeats the records.

- Flag filtering uses `ScanBamFlag` semantics by converting logical flag
  requirements into required-set and required-unset bit masks.

- Tag values are returned with their native BAM-derived type, matching
  `scanBam()`: integer-valued tags (`c`, `C`, `s`, `S`, `i`, `I`) become
  integer columns (or double when a value exceeds the R integer range),
  floating-point (`f`) tags become double columns, and `A`/`Z`/`H`/`B`
  tags become character columns. `B` tags are comma-separated character
  values. Absent tag values are `NA`.

- `seqqual_mode = "compact"` is optimized for throughput-oriented
  benchmarking and returns raw list-columns for `seq`/`qual`, not
  ordinary sequence or quality strings. In this representation, `seq`
  contains BAM-packed nucleotide bytes and `qual` contains raw Phred
  bytes. Compact output is intended for users who want to defer or avoid
  full string-materialization costs; use
  [`decode_compact_seq()`](https://cparsania.github.io/BamScaleR/reference/decode_compact_seq.md),
  [`decode_compact_qual()`](https://cparsania.github.io/BamScaleR/reference/decode_compact_qual.md),
  or
  [`decode_seqqual_compact()`](https://cparsania.github.io/BamScaleR/reference/decode_seqqual_compact.md)
  to decode compact output back to standard string form when needed.

- `"GAlignments"` and `"GAlignmentPairs"` output exclude unmapped
  records.

- For plain `as = "GAlignments"` requests (no `which`, tags, `seq` or
  `qual`), the object is assembled from slots built in C++. Every
  invariant checked by the `GAlignments` validity method holds by
  construction, except for two kinds of record that the reader counts: a
  CIGAR with an op code above 8 or a zero-length op (which the validity
  method may reject), and an alignment that extends outside its
  sequence. If either occurs,
  [`methods::validObject()`](https://rdrr.io/r/methods/validObject.html)
  is run on the result, so an invalid object is never returned and
  out-of-bound alignments give the same warning as
  [`GenomicAlignments::readGAlignments()`](https://rdrr.io/pkg/GenomicAlignments/man/readGAlignments.html).
  Set `options(BamScaleR.ga_fastpath = FALSE)` to always use the
  validating `GAlignments` constructor instead.

- For `as = "GAlignmentPairs"`, mates are matched by
  [`GenomicAlignments::makeGAlignmentPairs()`](https://rdrr.io/pkg/GenomicAlignments/man/findMateAlignment.html),
  on the read name, each record's position against its mate's
  `mrnm`/`mpos`, the strand, first/last, proper-pair and secondary
  flags. Without `which`, the result contains the same pairs as
  [`GenomicAlignments::readGAlignmentPairs()`](https://rdrr.io/pkg/GenomicAlignments/man/readGAlignments.html),
  ordered by the position of the first mate (pairs whose first mates
  start at the same position may come in a different order). With
  `which`, records are filtered before pairing, so only pairs with both
  mates in the ranges are returned, each once, whereas
  `readGAlignmentPairs()` returns every pair with at least one mate in a
  range. Ambiguous mates are dropped, and `BamFile` read-name options
  (`qnamePrefixEnd`, `qnameSuffixStart`) are not used. The mates' mcols
  hold the requested `what` fields and tags that are not alignment
  slots, as for `as = "GAlignments"`.

- `as = "scanBam"` returns a strict scan-like list-of-lists: without
  `param$which`, it returns one unnamed batch; with `param$which`, it
  returns one batch per range label (including empty ranges), with
  requested `what` fields and `tag` values under `$tag`. In this output
  mode, `seq` and `qual` are returned as
  [`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  and
  [`Biostrings::PhredQuality`](https://rdrr.io/pkg/Biostrings/man/XStringQuality-class.html)
  for closer `scanBam()` compatibility.

## Examples

``` r
bam <- ompBAM::example_BAM("Unsorted")

# Familiar scanBam-like field selection
x <- bam_read(bam, what = c("qname", "flag", "rname", "pos", "cigar"))

# Include sequence + quality
y <- bam_read(bam, what = c("qname", "seq", "qual"), threads = 2)

# scanBam-shaped output
z <- bam_read(bam, what = c("qname", "flag"), tag = "NM", as = "scanBam")
```
