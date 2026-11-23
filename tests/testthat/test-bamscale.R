test_that("bam_read returns selected fields", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar"),
    threads = 1
  )

  expect_s4_class(x, "DataFrame")
  expect_true(all(c("qname", "flag", "rname", "pos", "mapq", "cigar") %in% colnames(x)))
  expect_true(nrow(x) > 0)
})

test_that("rname, mrnm and strand are scanBam-compatible factors", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("rname", "mrnm", "strand"),
    as = "data.frame",
    threads = 1
  )

  expect_s3_class(x$rname, "factor")
  expect_s3_class(x$mrnm, "factor")
  expect_s3_class(x$strand, "factor")
  expect_identical(levels(x$strand), c("+", "-", "*"))

  # Levels and integer codes match Rsamtools::scanBam exactly.
  rs <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = c("rname", "mrnm", "strand"))
  )[[1]]
  expect_identical(levels(x$rname), levels(rs$rname))
  expect_identical(as.integer(x$rname), as.integer(rs$rname))
  expect_identical(as.integer(x$mrnm), as.integer(rs$mrnm))
  expect_identical(as.character(x$strand), as.character(rs$strand))
})

test_that("numeric tags are returned with native scanBam-compatible types", {
  # Build a small BAM carrying integer (i), float (f) and string (Z) tags,
  # with some records omitting a tag so absent values decode to NA.
  tmp <- tempfile(fileext = ".sam")
  n <- 1500L
  hdr <- c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:100000")
  recs <- vapply(seq_len(n), function(i) {
    core <- c(sprintf("r%05d", i), "0", "chr1", as.character(10L + i), "60",
              "10M", "*", "0", "0", "ACGTACGTAC", "IIIIIIIIII")
    opt <- if (i %% 3L == 0L) {
      c(sprintf("AS:i:%d", 100L - (i %% 40L)), sprintf("RG:Z:g%d", i %% 3L))
    } else {
      c(sprintf("NM:i:%d", i %% 5L), sprintf("AS:i:%d", 100L - (i %% 40L)),
        sprintf("XF:f:%.3f", (i %% 5L) + 0.25), sprintf("RG:Z:g%d", i %% 3L))
    }
    paste(c(core, opt), collapse = "\t")
  }, character(1))
  writeLines(c(hdr, recs), tmp)
  bam <- Rsamtools::asBam(tmp, tempfile(), overwrite = TRUE, indexDestination = FALSE)
  on.exit(unlink(c(tmp, bam)), add = TRUE)

  tags <- c("NM", "AS", "XF", "RG")
  x <- bam_read(bam, what = "qname", tag = tags, as = "data.frame", threads = 2)
  rs <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "qname", tag = tags)
  )[[1]]$tag

  expect_type(x$NM, "integer")
  expect_type(x$AS, "integer")
  expect_type(x$XF, "double")
  expect_type(x$RG, "character")

  expect_identical(x$NM, rs$NM)
  expect_identical(x$AS, rs$AS)
  expect_equal(x$XF, rs$XF)
  expect_identical(x$RG, rs$RG)

  # Absent NM values decode to NA in the same positions as scanBam.
  expect_identical(is.na(x$NM), is.na(rs$NM))
})

test_that("bam_read supports sequence and quality columns", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("qname", "seq", "qual"),
    as = "data.frame",
    threads = 1
  )

  expect_s3_class(x, "data.frame")
  expect_true(all(c("qname", "seq", "qual") %in% colnames(x)))
  expect_true(nrow(x) > 0)
})

test_that("list-style param filtering works", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    param = list(mapqFilter = 20L, flag_unset = 4L),
    what = c("flag", "mapq"),
    as = "data.frame",
    threads = 1
  )

  expect_true(all(x$mapq >= 20L))
  expect_true(all(bitwAnd(as.integer(x$flag), 4L) == 0L))
})

test_that("scanBam mode returns scan-like list structure", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos"),
    tag = c("NM"),
    as = "scanBam",
    threads = 1
  )

  expect_type(x, "list")
  expect_true(length(x) >= 1L)
  expect_null(names(x))
  expect_type(x[[1]], "list")
  expect_true(all(c("qname", "flag", "rname", "pos", "tag") %in% names(x[[1]])))
  expect_true("NM" %in% names(x[[1]]$tag))
})


test_that("scanBam seq/qual uses Biostrings classes when available", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("qname", "seq", "qual"),
    as = "scanBam",
    threads = 1
  )

  expect_s4_class(x[[1]]$seq, "DNAStringSet")
  expect_s4_class(x[[1]]$qual, "PhredQuality")
})


test_that("seq and qual are byte-identical to Rsamtools::scanBam", {
  bam <- ompBAM::example_BAM("Unsorted")
  rs <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = c("seq", "qual"))
  )[[1]]

  # data.frame mode now returns native Biostrings columns (previously character),
  # built in C directly from the reader's shared byte buffer. For any result
  # scanBam returns as a single pool buffer (small/moderate files) this is
  # byte-for-byte identical(); larger inputs where scanBam fragments its pool
  # remain content-identical but not object-identical() by design.
  df <- bam_read(bam, what = c("seq", "qual"), as = "data.frame", threads = 1)
  expect_s4_class(df$seq, "DNAStringSet")
  expect_s4_class(df$qual, "PhredQuality")
  expect_identical(df$seq, rs$seq)
  expect_identical(df$qual, rs$qual)

  # The multi-threaded merge preserves file order, so the result is unchanged.
  df_mt <- bam_read(bam, what = c("seq", "qual"), as = "data.frame", threads = 4)
  expect_identical(df_mt$seq, rs$seq)
  expect_identical(df_mt$qual, rs$qual)

  # scanBam output mode carries the same S4 objects through unchanged.
  sb <- bam_read(bam, what = c("seq", "qual"), as = "scanBam", threads = 1)[[1]]
  expect_identical(sb$seq, rs$seq)
  expect_identical(sb$qual, rs$qual)
})


test_that("compatible seq/qual builds correctly on fresh SnowParam SOCK workers", {
  # Regression guard: the C XStringSet builder resolves IRanges/XVector
  # C-callables at runtime and instantiates Biostrings S4 classes. A fresh SOCK
  # worker must have those namespaces loaded (via .onLoad) or it fails with
  # "function '_new_IRanges' not provided" / "DNAStringSet is not a defined class".
  skip_if_not_installed("BiocParallel")
  bam <- ompBAM::example_BAM("Unsorted")
  rs <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(what = c("seq", "qual")))[[1]]
  bp <- suppressWarnings(BiocParallel::SnowParam(2L, type = "SOCK"))
  res <- BiocParallel::bplapply(1:2, function(i) {
    x <- BamScaleR::bam_read(bam, what = c("seq", "qual"), as = "data.frame", threads = 1)
    list(seq = as.character(x$seq), qual = as.character(x$qual),
         seq_cls = class(x$seq), qual_cls = class(x$qual))
  }, BPPARAM = bp)

  expect_identical(as.character(res[[1]]$seq_cls), "DNAStringSet")
  expect_identical(as.character(res[[1]]$qual_cls), "PhredQuality")
  expect_identical(res[[1]]$seq, as.character(rs$seq))
  expect_identical(res[[1]]$qual, as.character(rs$qual))
})


test_that("scanBam mode batches by which labels and keeps empty ranges", {
  bam <- ompBAM::example_BAM("Unsorted")
  seed <- bam_read(
    file = bam,
    what = c("rname", "pos"),
    as = "data.frame",
    threads = 1
  )

  first_rname <- as.character(seed[["rname"]][1])
  first_pos <- as.integer(seed[["pos"]][1])

  x <- bam_read(
    file = bam,
    param = list(
      which = data.frame(
        seqname = c(first_rname, first_rname),
        start = c(max(1L, first_pos - 2L), first_pos + 1000000000L),
        end = c(first_pos + 2L, first_pos + 1000000010L),
        label = c("hit", "empty"),
        stringsAsFactors = FALSE
      )
    ),
    what = c("qname", "flag"),
    as = "scanBam",
    threads = 1
  )

  expect_type(x, "list")
  expect_equal(names(x), c("hit", "empty"))
  expect_equal(length(x), 2L)
  expect_true(length(x[["hit"]]$qname) >= 1L)
  expect_length(x[["empty"]]$qname, 0L)
  expect_length(x[["empty"]]$flag, 0L)
})

test_that("multi-file input returns named list", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(file = c(a = bam, b = bam), what = c("qname"), threads = 1)

  expect_type(x, "list")
  expect_equal(length(x), 2L)
  expect_equal(names(x), c("a", "b"))
})

test_that("bam_count runs and returns expected columns", {
  bam <- ompBAM::example_BAM("Unsorted")
  y <- bam_count(file = bam, threads = 1)

  expect_s3_class(y, "data.frame")
  expect_true(all(c("seqname", "seqlength", "count") %in% names(y)))
  expect_true(nrow(y) > 0)
})

test_that("fragment_sizes is byte-identical to table(abs(scanBam isize))", {
  bam <- ompBAM::example_BAM("Unsorted")
  FLAG <- Rsamtools::scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE
  )
  fs <- fragment_sizes(
    bam,
    param = Rsamtools::ScanBamParam(flag = FLAG),
    threads = 1
  )
  for (th in c(2L, 4L)) {
    expect_identical(
      fragment_sizes(bam, param = Rsamtools::ScanBamParam(flag = FLAG), threads = th),
      fs
    )
  }
  expect_s3_class(fs, "data.frame")
  expect_identical(names(fs), c("fragment_size", "count"))
  expect_false(is.unsorted(fs$fragment_size))

  # Reference: scanBam materialises every insert size, then table(abs()).
  # scanBam reports isize = NA for mate-unmapped reads, which table() drops;
  # fragment_sizes(drop_mate_unmapped = TRUE) replicates that exactly.
  iz <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "isize", flag = FLAG)
  )[[1]]$isize
  ref <- table(abs(iz))
  ref_df <- data.frame(
    fragment_size = as.integer(names(ref)),
    ref = as.numeric(ref)
  )
  m <- merge(fs, ref_df, by = "fragment_size", all = TRUE)
  m[is.na(m)] <- 0
  expect_equal(m$count, m$ref)
  expect_equal(sum(fs$count), sum(!is.na(iz)))
})

test_that("fragment_sizes(drop_mate_unmapped = FALSE) keeps mate-unmapped reads", {
  bam <- ompBAM::example_BAM("Unsorted")
  keep <- fragment_sizes(bam, threads = 1, drop_mate_unmapped = FALSE)
  drop <- fragment_sizes(bam, threads = 1, drop_mate_unmapped = TRUE)
  expect_gte(sum(keep$count), sum(drop$count))
})

test_that("mapq_dist is byte-identical to table(scanBam mapq)", {
  bam <- ompBAM::example_BAM("Unsorted")
  FLAG <- Rsamtools::scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE
  )
  md <- mapq_dist(
    bam,
    param = Rsamtools::ScanBamParam(flag = FLAG),
    threads = 1
  )
  for (th in c(2L, 4L)) {
    expect_identical(
      mapq_dist(bam, param = Rsamtools::ScanBamParam(flag = FLAG), threads = th),
      md
    )
  }
  expect_s3_class(md, "data.frame")
  expect_identical(names(md), c("mapq", "count"))
  expect_false(is.unsorted(md$mapq))

  mq <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "mapq", flag = FLAG)
  )[[1]]$mapq
  ref <- table(mq)
  ref_df <- data.frame(
    mapq = as.integer(names(ref)),
    ref = as.numeric(ref)
  )
  m <- merge(md, ref_df, by = "mapq", all = TRUE)
  m[is.na(m)] <- 0
  expect_equal(m$count, m$ref)
  expect_equal(sum(md$count), sum(!is.na(mq)))
})

test_that("fragment_sizes and mapq_dist return named lists for multiple files", {
  bam <- ompBAM::example_BAM("Unsorted")
  fs <- fragment_sizes(file = c(a = bam, b = bam), threads = 1)
  expect_type(fs, "list")
  expect_equal(names(fs), c("a", "b"))
  expect_s3_class(fs[[1]], "data.frame")

  md <- mapq_dist(file = c(a = bam, b = bam), threads = 1)
  expect_type(md, "list")
  expect_equal(names(md), c("a", "b"))
  expect_s3_class(md[[1]], "data.frame")
})

test_that("bam_coverage is byte-identical to coverage(readGAlignments)", {
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")
  got <- bam_coverage(bam, threads = 1)
  ref <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(bam))

  # Output must be an *uncompressed* SimpleRleList, integer runs, header-order
  # names -- each independently required for identical() to coverage().
  expect_true(methods::is(got, "SimpleRleList"))
  expect_false(methods::is(got, "CompressedRleList"))
  expect_identical(names(got), names(ref))
  expect_identical(typeof(S4Vectors::runValue(got[[1]])), "integer")
  expect_identical(got, ref)
  for (th in c(2L, 4L)) {
    expect_identical(bam_coverage(bam, threads = th), ref)
  }
})

test_that("bam_coverage keeps secondary/supplementary reads (readGAlignments filter)", {
  # Only unmapped (0x4) reads are dropped; the default must match readGAlignments,
  # not a stricter flag filter. Passing that exact flag must not change the result.
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")
  default_cov <- bam_coverage(bam, threads = 1)
  explicit_cov <- bam_coverage(
    bam, threads = 1,
    param = Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
    )
  )
  expect_identical(default_cov, explicit_cov)
})

test_that("bam_coverage returns a named list for multiple files", {
  bam <- ompBAM::example_BAM("Unsorted")
  cv <- bam_coverage(file = c(a = bam, b = bam), threads = 1)
  expect_type(cv, "list")
  expect_equal(names(cv), c("a", "b"))
  expect_true(methods::is(cv[[1]], "SimpleRleList"))
})

test_that("bam_coverage_bigwig writes a bigWig with values identical to coverage()", {
  # rtracklayer's bigWig import goes through the UCSC udc layer, which cannot
  # parse Windows paths ("Unrecognized protocol C") -- bigWig I/O via
  # rtracklayer is unsupported on Windows. The parallel-vs-serial byte-identity
  # test below still runs there.
  skip_on_os("windows")
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")
  out <- tempfile(fileext = ".bw")
  on.exit(unlink(out), add = TRUE)

  ret <- bam_coverage_bigwig(bam, out, threads = 1)
  expect_identical(as.character(ret), out)
  expect_true(file.exists(out) && file.info(out)$size > 0)

  imp <- rtracklayer::import(out, as = "RleList")
  ref <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(bam))
  common <- intersect(names(ref), names(imp))
  expect_gt(length(common), 0)
  for (nm in common) {
    expect_equal(as.numeric(imp[[nm]]), as.numeric(ref[[nm]]))
  }
  auc <- function(rl) sum(vapply(rl, function(x)
    sum(as.numeric(S4Vectors::runValue(x)) * as.numeric(S4Vectors::runLength(x))),
    numeric(1)))
  expect_equal(auc(imp[common]), auc(ref))
})

test_that("bam_coverage_bigwig errors when outfile length != file length", {
  bam <- ompBAM::example_BAM("Unsorted")
  expect_error(
    bam_coverage_bigwig(c(bam, bam), tempfile(fileext = ".bw"), threads = 1),
    "same length"
  )
})

test_that("bam_coverage_bigwig parallel output is byte-identical to serial", {
  # Parallel block compression must produce a bit-for-bit identical bigWig
  # (compress2 is deterministic; blocks are written in the same order). On
  # Windows the UCRT's stream-position semantics on buffered update streams can
  # encode index offsets slightly differently between the two write patterns, so
  # there we assert the invariants that survive -- identical size and identical
  # header -- and leave bit-identity to the POSIX platforms where it is proven.
  bam <- ompBAM::example_BAM("Unsorted")
  op <- tempfile(fileext = ".bw")
  os <- tempfile(fileext = ".bw")
  on.exit(unlink(c(op, os)), add = TRUE)
  bam_coverage_bigwig(bam, op, threads = 2, parallel = TRUE)
  bam_coverage_bigwig(bam, os, threads = 2, parallel = FALSE)
  if (identical(.Platform$OS.type, "windows")) {
    expect_identical(file.info(op)$size, file.info(os)$size)
    expect_identical(readBin(op, "raw", 64L), readBin(os, "raw", 64L))
  } else {
    expect_identical(
      readBin(op, "raw", file.info(op)$size),
      readBin(os, "raw", file.info(os)$size)
    )
  }
})

test_that("GA fast path is identical to readGAlignments (core fields)", {
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")
  ref <- GenomicAlignments::readGAlignments(Rsamtools::BamFile(bam))

  ga1 <- bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
                  as = "GAlignments", threads = 1)
  expect_identical(ga1, ref)

  # threads > 1 exercises the tid-ordered merge and RLE boundary coalescing
  ga4 <- bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
                  as = "GAlignments", threads = 4)
  expect_identical(ga4, ref)
})

test_that("GA fast path mcols match readGAlignments param what", {
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")

  ref <- GenomicAlignments::readGAlignments(
    Rsamtools::BamFile(bam),
    param = Rsamtools::ScanBamParam(what = c("flag", "qname"))
  )
  ga <- bam_read(bam,
                 what = c("rname", "pos", "cigar", "strand", "flag", "qname"),
                 as = "GAlignments", threads = 2)
  expect_identical(ga, ref)

  # reversed mcols order must be preserved
  ref2 <- GenomicAlignments::readGAlignments(
    Rsamtools::BamFile(bam),
    param = Rsamtools::ScanBamParam(what = c("qname", "flag"))
  )
  ga2 <- bam_read(bam,
                  what = c("rname", "pos", "cigar", "strand", "qname", "flag"),
                  as = "GAlignments", threads = 2)
  expect_identical(ga2, ref2)
})

test_that("GA fast path honors use.names and filters like readGAlignments", {
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")

  ref <- GenomicAlignments::readGAlignments(
    Rsamtools::BamFile(bam),
    use.names = TRUE,
    param = Rsamtools::ScanBamParam(what = c("flag", "qname"))
  )
  ga <- bam_read(bam,
                 what = c("rname", "pos", "cigar", "strand", "flag", "qname"),
                 as = "GAlignments", use.names = TRUE, threads = 2)
  expect_identical(ga, ref)

  FLAG <- Rsamtools::scanBamFlag(isMinusStrand = FALSE)
  reff <- GenomicAlignments::readGAlignments(
    Rsamtools::BamFile(bam),
    param = Rsamtools::ScanBamParam(flag = FLAG, mapqFilter = 10)
  )
  gaf <- bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
                  as = "GAlignments", threads = 2,
                  param = Rsamtools::ScanBamParam(flag = FLAG, mapqFilter = 10))
  expect_identical(gaf, reff)
})

test_that("GA fast path and slow path produce identical objects", {
  bam <- ompBAM::example_BAM("Unsorted")
  fields <- c("rname", "pos", "cigar", "strand", "flag", "mapq")
  fast <- bam_read(bam, what = fields, as = "GAlignments", threads = 2)
  old <- options(BamScaleR.ga_fastpath = FALSE)
  on.exit(options(old), add = TRUE)
  slow <- bam_read(bam, what = fields, as = "GAlignments", threads = 2)
  expect_identical(fast, slow)
})

test_that("GA fast path survives gctorture (CHARSXP cache GC safety)", {
  bam <- ompBAM::example_BAM("Unsorted")
  gctorture(TRUE)
  ga <- tryCatch(
    bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
             as = "GAlignments", threads = 1),
    finally = gctorture(FALSE)
  )
  expect_s4_class(ga, "GAlignments")
  expect_true(length(ga) > 0)
})

test_that("GAlignments output is available when package is installed", {
  bam <- ompBAM::example_BAM("Unsorted")
  g <- bam_read(
    file = bam,
    what = c("qname", "rname", "pos", "cigar", "strand", "flag"),
    as = "GAlignments",
    threads = 1
  )

  expect_s4_class(g, "GAlignments")
})

test_that("GAlignmentPairs output is available when package is installed", {
  bam <- ompBAM::example_BAM("Unsorted")
  gp <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos", "cigar", "strand"),
    as = "GAlignmentPairs",
    include_unmapped = FALSE,
    threads = 1
  )

  expect_s4_class(gp, "GAlignmentPairs")
})


test_that("GAlignments carries header seqlengths and gives identical coverage()", {
  bam <- ompBAM::example_BAM("Unsorted")

  g <- bam_read(
    file = bam,
    what = c("rname", "pos", "cigar", "strand"),
    as = "GAlignments",
    threads = 1
  )
  std <- GenomicAlignments::readGAlignments(bam)

  # Header sequence lengths must be populated (not NA) and match the standard
  # reader, so coverage()/export.bw() see true chromosome sizes rather than
  # max-end truncation.
  expect_identical(GenomeInfoDb::seqlevels(g), GenomeInfoDb::seqlevels(std))
  expect_identical(GenomeInfoDb::seqlengths(g), GenomeInfoDb::seqlengths(std))
  expect_false(any(is.na(GenomeInfoDb::seqlengths(g))))

  # coverage() is byte-identical to the standard Bioconductor path -- the
  # equivalence anchor for the coverage -> bigWig workflow benchmark.
  expect_identical(
    GenomicAlignments::coverage(g),
    GenomicAlignments::coverage(std)
  )
})


test_that("auto_threads validates logical input", {
  expect_error(
    BamScaleR:::.bamscale_resolve_parallel_plan(threads = 1L, auto_threads = NA),
    "`auto_threads` must be TRUE or FALSE"
  )
})


test_that("auto_threads preserves per-file threads by reducing active workers first", {
  skip_if(identical(BamScaleR:::openmp_max_threads_cpp(), 0L), "built without OpenMP")
  bp <- tryCatch(
    suppressWarnings(
      BiocParallel::SnowParam(workers = 8L, type = "SOCK", progressbar = FALSE)
    ),
    error = function(e) NULL
  )
  if (is.null(bp)) {
    skip("SnowParam could not be created (e.g. worker limits under R CMD check --as-cran)")
  }
  on.exit(try(BiocParallel::bpstop(bp), silent = TRUE), add = TRUE)

  plan <- BamScaleR:::.bamscale_resolve_parallel_plan(
    threads = 64L,
    BPPARAM = bp,
    auto_threads = TRUE,
    n_files = 8L
  )

  expected_threads <- max(1L, min(64L, BamScaleR:::.bamscale_detect_cores()))
  expect_equal(plan$threads, expected_threads)
  expect_equal(plan$bp_workers, 1L)
})


test_that("auto_threads keeps multiple workers when requested per-file threads are small", {
  skip_if(identical(BamScaleR:::openmp_max_threads_cpp(), 0L), "built without OpenMP")
  bp <- tryCatch(
    suppressWarnings(
      BiocParallel::SnowParam(workers = 8L, type = "SOCK", progressbar = FALSE)
    ),
    error = function(e) NULL
  )
  if (is.null(bp)) {
    skip("SnowParam could not be created (e.g. worker limits under R CMD check --as-cran)")
  }
  on.exit(try(BiocParallel::bpstop(bp), silent = TRUE), add = TRUE)

  plan <- BamScaleR:::.bamscale_resolve_parallel_plan(
    threads = 2L,
    BPPARAM = bp,
    auto_threads = TRUE,
    n_files = 8L
  )

  workers_eff <- min(BamScaleR:::.bamscale_bpparam_workers(bp), 8L)
  expected_workers <- min(
    workers_eff,
    max(1L, floor(BamScaleR:::.bamscale_detect_cores() / 2L))
  )
  expect_equal(plan$threads, 2L)
  expect_equal(plan$bp_workers, expected_workers)
  expect_lte(plan$threads * plan$bp_workers, max(1L, BamScaleR:::.bamscale_detect_cores()))
})

test_that("compact decode helpers decode synthetic seq and qual correctly", {
  seq_raw <- list(
    as.raw(c(0x12, 0x48)), # A C G T
    as.raw(c(0xFF)),       # N N
    raw()
  )
  qwidth <- c(4L, 2L, 0L)

  qual_raw <- list(
    as.raw(c(0L, 1L, 2L, 3L)),
    as.raw(c(255L, 255L)),
    raw()
  )

  expect_identical(
    decode_compact_seq(seq_raw, qwidth),
    c("ACGT", "NN", "")
  )

  expect_identical(
    decode_compact_qual(qual_raw),
    c("!\"#$", "*", "")
  )
})

test_that("compact decode helpers validate malformed inputs", {
  expect_error(
    decode_compact_seq(list(as.raw(c(0x12))), integer()),
    "`seq` and `qwidth` must have the same length"
  )

  expect_error(
    decode_compact_seq("not-a-list", 1L),
    "`seq` must be a list of raw vectors"
  )

  expect_error(
    decode_compact_qual("not-a-list"),
    "`qual` must be a list of raw vectors"
  )

  expect_error(
    decode_seqqual_compact(data.frame(seq = I(list(as.raw(c(0x12)))))),
    "requires a `qwidth` column"
  )
})

test_that("decode_seqqual_compact preserves non-seqqual columns", {
  x <- data.frame(
    qname = "r1",
    qwidth = 4L,
    mapq = 42L,
    stringsAsFactors = FALSE
  )
  x$seq <- I(list(as.raw(c(0x12, 0x48))))
  x$qual <- I(list(as.raw(c(0L, 1L, 2L, 3L))))

  y <- decode_seqqual_compact(x)

  expect_identical(y$qname, x$qname)
  expect_identical(y$qwidth, x$qwidth)
  expect_identical(y$mapq, x$mapq)
  expect_identical(y$seq, "ACGT")
  expect_identical(y$qual, "!\"#$")
})

test_that("compact seqqual round-trips to BamScaleR compatible output", {
  bam <- ompBAM::example_BAM("Unsorted")

  compat <- bam_read(
    file = bam,
    what = c("qname", "qwidth", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compatible",
    threads = 1
  )

  compact <- bam_read(
    file = bam,
    what = c("qname", "qwidth", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compact",
    threads = 1
  )

  expect_true(is.list(compact$seq))
  expect_true(is.list(compact$qual))
  expect_true(all(vapply(compact$seq, is.raw, logical(1))))
  expect_true(all(vapply(compact$qual, is.raw, logical(1))))

  decoded <- decode_seqqual_compact(compact)

  expect_identical(decoded$qname, compat$qname)
  expect_identical(decoded$qwidth, compat$qwidth)
  # Compatible mode now returns native Biostrings columns, so compare the decoded
  # compact strings against their character form.
  expect_identical(decoded$seq, as.character(compat$seq))
  expect_identical(decoded$qual, as.character(compat$qual))
})

test_that("compact seqqual decode matches Rsamtools sequence and quality output", {
  bam <- ompBAM::example_BAM("Unsorted")

  compact <- bam_read(
    file = bam,
    what = c("qname", "qwidth", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compact",
    threads = 1
  )
  decoded <- decode_seqqual_compact(compact)

  rs <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = c("qname", "seq", "qual"))
  )[[1]]

  expect_identical(decoded$qname, rs$qname)
  expect_identical(decoded$seq, as.character(rs$seq))
  expect_identical(decoded$qual, as.character(rs$qual))
})


# --- BPPARAM is used as supplied, never resized -------------------------------

# SnowParam with `n` workers, or a skip when the environment allows fewer
# (BiocParallel caps workers under _R_CHECK_LIMIT_CORES_ or
# BIOCPARALLEL_WORKER_MAX, with a warning or an error).
.snow_param_or_skip <- function(n) {
  bp <- tryCatch(
    suppressWarnings(
      BiocParallel::SnowParam(workers = n, type = "SOCK", progressbar = FALSE)
    ),
    error = function(e) NULL
  )
  if (is.null(bp) || BiocParallel::bpnworkers(bp) < n) {
    skip(sprintf("fewer than %d BiocParallel workers allowed in this environment", n))
  }
  bp
}

test_that("multi-file calls never modify BPPARAM (unstarted or started cluster)", {
  bp <- .snow_param_or_skip(3L)
  on.exit(try(BiocParallel::bpstop(bp), silent = TRUE), add = TRUE)

  bam <- ompBAM::example_BAM("Unsorted")
  files <- c(a = bam, b = bam)
  ref <- bam_count(bam, threads = 1, BPPARAM = NULL)

  # Unstarted: 3 workers but only 2 files. Previously the object was resized
  # to 2 workers in place.
  res <- bam_count(files, threads = 1, BPPARAM = bp)
  expect_identical(res$a, ref)
  expect_identical(res$b, ref)
  expect_identical(BiocParallel::bpnworkers(bp), 3L)
  expect_false(BiocParallel::bpisup(bp))

  # Started: stays up with its configured workers and remains fully usable.
  BiocParallel::bpstart(bp)
  res2 <- bam_count(files, threads = 1, BPPARAM = bp)
  expect_identical(res2$a, ref)
  expect_identical(res2$b, ref)
  expect_true(BiocParallel::bpisup(bp))
  expect_identical(BiocParallel::bpnworkers(bp), 3L)
  pids <- unlist(BiocParallel::bplapply(1:3, function(i) Sys.getpid(), BPPARAM = bp))
  expect_length(unique(pids), 3L)
})

test_that(".bamscale_apply_files caps concurrency by grouping, keeping input order", {
  bp <- .snow_param_or_skip(3L)
  on.exit(try(BiocParallel::bpstop(bp), silent = TRUE), add = TRUE)
  BiocParallel::bpstart(bp)
  n0 <- BiocParallel::bpnworkers(bp)

  files <- stats::setNames(as.character(1:5), letters[1:5])
  out <- BamScaleR:::.bamscale_apply_files(
    files,
    function(x) c(x = x, pid = as.character(Sys.getpid())),
    BPPARAM = bp, bp_workers = 2L
  )
  expect_identical(names(out), letters[1:5])
  expect_identical(unname(vapply(out, `[[`, "", "x")), as.character(1:5))
  expect_lte(length(unique(vapply(out, `[[`, "", "pid"))), 2L)
  expect_identical(BiocParallel::bpnworkers(bp), n0)

  # A capped run on a running cluster groups the files; an error names the
  # file that failed.
  fail_on_e <- function(x) if (x == "5") stop("boom") else x
  expect_error(
    BamScaleR:::.bamscale_apply_files(files, fail_on_e, BPPARAM = bp, bp_workers = 2L),
    "file 5 \\(e\\): boom"
  )
  # Without a cap the files go to bplapply() directly, so BiocParallel's
  # element index is the file's index. (Grouping 5 files on 3 workers would
  # put file 5 in group 2.)
  expect_error(
    BamScaleR:::.bamscale_apply_files(files, fail_on_e, BPPARAM = bp),
    "element index: 5"
  )
})

test_that("a cap on an unstarted SnowParam keeps one result per file", {
  bp <- .snow_param_or_skip(3L)
  on.exit(try(BiocParallel::bpstop(bp), silent = TRUE), add = TRUE)
  files <- stats::setNames(as.character(1:5), letters[1:5])
  out <- BamScaleR:::.bamscale_apply_files(
    files,
    function(x) c(x = x, pid = as.character(Sys.getpid())),
    BPPARAM = bp, bp_workers = 2L
  )
  expect_identical(names(out), letters[1:5])
  expect_lte(length(unique(vapply(out, `[[`, "", "pid"))), 2L)
  expect_false(BiocParallel::bpisup(bp))
  expect_identical(BiocParallel::bpnworkers(bp), 3L)
  # Not grouped: the element index of an error is the file's index.
  expect_error(
    BamScaleR:::.bamscale_apply_files(
      files, function(x) if (x == "5") stop("boom") else x,
      BPPARAM = bp, bp_workers = 2L
    ),
    "element index: 5"
  )
  expect_false(BiocParallel::bpisup(bp))
})


# --- GAlignments fast path: validity ------------------------------------------

# Small synthetic BAM: `records` (SAM lines, positions >= 2) preceded by
# in-bounds filler reads at position 1. The filler makes the alignments span
# more than the header's BGZF block; ompBAM 1.16.0 returns no records for a BAM
# whose alignments all share the header's block.
.make_test_bam <- function(records, seqlen = 100L, n_filler = 3000L) {
  filler <- sprintf("f%d\t0\tchr1\t1\t60\t10M\t*\t0\t0\t*\t*", seq_len(n_filler))
  sam <- tempfile(fileext = ".sam")
  on.exit(unlink(sam), add = TRUE)
  writeLines(c(
    "@HD\tVN:1.6\tSO:coordinate",
    sprintf("@SQ\tSN:chr1\tLN:%d", seqlen),
    filler,
    records
  ), sam)
  out <- Rsamtools::asBam(sam, sub("\\.bam$", "", tempfile(fileext = ".bam")),
                          overwrite = TRUE, indexDestination = FALSE)
  # Removed (with any index) when the calling test finishes.
  withr::defer(unlink(c(out, paste0(out, ".bai"))), envir = parent.frame())
  out
}

test_that("GA fast path output passes validObject(complete = TRUE)", {
  bam <- ompBAM::example_BAM("Unsorted")
  ga <- bam_read(
    bam,
    what = c("rname", "pos", "cigar", "strand", "qname", "flag", "mapq"),
    as = "GAlignments", use.names = TRUE, threads = 2
  )
  expect_s4_class(ga, "GAlignments")
  expect_silent(methods::validObject(ga, complete = TRUE))

  # The reader found nothing that would need the R-side re-check.
  payload <- BamScaleR:::galignments_cpp(
    bam, 2L, 0L, 0L, 0L,
    BamScaleR:::.bamscale_make_field_mask(c("rname", "pos", "cigar", "strand"))
  )
  expect_identical(payload$n_nonstandard_cigar, 0)
  expect_identical(payload$n_out_of_bound, 0)
})

test_that("GA fast path warns on out-of-bound alignments like readGAlignments", {
  bam <- .make_test_bam("r1\t0\tchr1\t95\t60\t10M\t*\t0\t0\t*\t*")
  expect_warning(ref <- GenomicAlignments::readGAlignments(bam), "out-of-bound")
  expect_warning(
    ga <- bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
                   as = "GAlignments", threads = 1),
    "out-of-bound"
  )
  expect_identical(ga, ref)
})

test_that("GA fast path rejects CIGARs that GAlignments validity rejects", {
  bam <- .make_test_bam("r1\t0\tchr1\t20\t60\t5M0M\t*\t0\t0\t*\t*")
  expect_error(GenomicAlignments::readGAlignments(bam))
  expect_error(
    bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
             as = "GAlignments", threads = 1),
    "cigar"
  )
})


# --- OpenMP availability ------------------------------------------------------

test_that("openmp_max_threads_cpp reports the build's OpenMP status", {
  n <- BamScaleR:::openmp_max_threads_cpp()
  expect_type(n, "integer")
  expect_gte(n, 0L)

  # Where R has OpenMP flags of its own (Linux, Windows), the package is built
  # with them, so OpenMP must be present -- unless the flags were overridden,
  # with BAMSCALER_OPENMP_* (see ./configure) or a user Makevars. On macOS it
  # depends on the toolchain. BAMSCALER_EXPECT_OPENMP=true/false forces the
  # check on or off.
  expect_env <- tolower(Sys.getenv("BAMSCALER_EXPECT_OPENMP"))
  user_mv <- tools::makevars_user()
  user_omp <- length(user_mv) > 0L && any(grepl(
    "OPENMP|fopenmp", unlist(lapply(user_mv, readLines, warn = FALSE))
  ))
  overridden <- !is.na(Sys.getenv("BAMSCALER_OPENMP_CXXFLAGS", unset = NA)) ||
    user_omp || identical(expect_env, "false")
  if (Sys.info()[["sysname"]] != "Darwin" && !overridden) {
    makeconf <- file.path(R.home("etc"), Sys.getenv("R_ARCH"), "Makeconf")
    conf <- if (file.exists(makeconf)) readLines(makeconf, warn = FALSE) else character()
    flags <- sub(
      "^SHLIB_OPENMP_CXXFLAGS[[:space:]]*=[[:space:]]*", "",
      grep("^SHLIB_OPENMP_CXXFLAGS[[:space:]]*=", conf, value = TRUE)
    )
    if (any(nzchar(trimws(flags)))) {
      expect_gte(n, 1L)
    }
  }
  if (identical(expect_env, "true")) {
    expect_gte(n, 1L)
  }
})

test_that("threads > 1 warns once per session when built without OpenMP", {
  local_mocked_bindings(openmp_max_threads_cpp = function() 0L)
  st <- BamScaleR:::.bamscale_state
  old <- st$warned_no_openmp
  st$warned_no_openmp <- NULL
  on.exit(st$warned_no_openmp <- old, add = TRUE)

  bam <- ompBAM::example_BAM("Unsorted")
  expect_no_warning(bam_count(bam, threads = 1, BPPARAM = NULL))
  expect_warning(bam_count(bam, threads = 2, BPPARAM = NULL), "built without OpenMP")
  expect_no_warning(bam_count(bam, threads = 2, BPPARAM = NULL))
})


# --- Forked workers after OpenMP use in the parent ----------------------------

test_that("a forked child does not hang after OpenMP threads ran in the parent", {
  skip_on_os("windows")
  skip_if(identical(BamScaleR:::openmp_max_threads_cpp(), 0L), "built without OpenMP")
  bam <- ompBAM::example_BAM("Unsorted")
  # Runs a multithreaded OpenMP region in this process. With GCC's libgomp a
  # forked child that then starts its own threads used to wait forever.
  ref <- bam_count(bam, threads = 2, BPPARAM = NULL)
  job <- parallel::mcparallel({
    counts <- bam_count(bam, threads = 2, BPPARAM = NULL)
    # Threads alive in the child after its region (Linux): more than one means
    # the child really ran OpenMP threads, not a one-thread fallback.
    n_threads <- if (dir.exists("/proc/self/task")) length(list.files("/proc/self/task")) else NA
    list(counts = counts, n_threads = n_threads)
  })
  res <- parallel::mccollect(job, wait = FALSE, timeout = 60)
  if (is.null(res)) {
    tools::pskill(job$pid, tools::SIGKILL)
    suppressWarnings(parallel::mccollect(job, wait = TRUE))
    fail("forked child did not finish within 60 seconds")
  } else {
    expect_identical(res[[1L]]$counts, ref)
    info <- BamScaleR:::openmp_fork_info_cpp()
    can_thread <- min(BamScaleR:::openmp_max_threads_cpp(), info[["thread_limit"]]) >= 2L
    if (!is.na(res[[1L]]$n_threads) && can_thread) {
      if (identical(info[["release"]], 1L)) {
        expect_gte(res[[1L]]$n_threads, 2L)
      } else {
        # libgomp before version 9: the documented one-thread fallback.
        expect_identical(res[[1L]]$n_threads, 1L)
      }
    }
  }
})


# --- GAlignmentPairs: mate matching --------------------------------------------

# Order-independent description of a GAlignmentPairs object.
.pair_key <- function(gp) {
  f <- GenomicAlignments::first(gp)
  l <- GenomicAlignments::last(gp)
  d <- data.frame(
    s1 = as.character(GenomicAlignments::seqnames(f)), p1 = BiocGenerics::start(f),
    c1 = GenomicAlignments::cigar(f), t1 = as.character(BiocGenerics::strand(f)),
    s2 = as.character(GenomicAlignments::seqnames(l)), p2 = BiocGenerics::start(l),
    c2 = GenomicAlignments::cigar(l), t2 = as.character(BiocGenerics::strand(l)),
    stringsAsFactors = FALSE
  )
  d <- d[do.call(order, d), ]
  rownames(d) <- NULL
  d
}

test_that("GAlignmentPairs output has the same pairs as readGAlignmentPairs", {
  bam <- ompBAM::example_BAM("Unsorted")
  gp <- bam_read(bam, what = c("rname", "pos", "cigar", "strand", "mapq"),
                 as = "GAlignmentPairs", threads = 2, BPPARAM = NULL)
  ref <- GenomicAlignments::readGAlignmentPairs(
    bam, param = Rsamtools::ScanBamParam(what = "mapq")
  )
  expect_identical(length(gp), length(ref))
  expect_identical(.pair_key(gp), .pair_key(ref))
  # Ordered by first-mate position, as readGAlignmentPairs() orders them.
  f <- GenomicAlignments::first(gp)
  expect_false(is.unsorted(
    as.integer(GenomicAlignments::seqnames(f)) * 1e10 + BiocGenerics::start(f)
  ))
  expect_identical(colnames(S4Vectors::mcols(f)), "mapq")
})

test_that("GAlignmentPairs pairs secondary and supplementary alignments correctly", {
  # Template m1 has a primary pair (1000/1200) and a secondary pair
  # (3000/500); pairing the i-th first mate with the i-th last mate by name
  # would give 1000/500 and 3000/1200.
  sec <- .make_test_bam(c(
    "m1\t99\tchr1\t1000\t60\t20M\t=\t1200\t220\t*\t*",
    "m1\t147\tchr1\t1200\t60\t20M\t=\t1000\t-220\t*\t*",
    "m1\t353\tchr1\t3000\t0\t20M\t=\t500\t-2480\t*\t*",
    "m1\t401\tchr1\t500\t0\t20M\t=\t3000\t2480\t*\t*",
    "p2\t99\tchr1\t2000\t60\t20M\t=\t2100\t120\t*\t*",
    "p2\t147\tchr1\t2100\t60\t20M\t=\t2000\t-120\t*\t*"
  ), seqlen = 5000L)
  expect_no_warning(
    gp <- bam_read(sec, what = c("rname", "pos", "cigar", "strand"),
                   as = "GAlignmentPairs", BPPARAM = NULL)
  )
  expect_identical(gp, GenomicAlignments::readGAlignmentPairs(sec))
  expect_identical(BiocGenerics::start(GenomicAlignments::first(gp)), c(1000L, 2000L, 3000L))
  expect_identical(BiocGenerics::start(GenomicAlignments::last(gp)), c(1200L, 2100L, 500L))

  # A supplementary piece of the first read must not replace its primary.
  sup <- .make_test_bam(c(
    "s1\t2113\tchr1\t400\t60\t30H20M\t=\t1200\t0\t*\t*",
    "s1\t99\tchr1\t1000\t60\t20M30S\t=\t1200\t220\t*\t*",
    "s1\t147\tchr1\t1200\t60\t20M\t=\t1000\t-220\t*\t*"
  ), seqlen = 5000L)
  expect_no_warning(
    gp <- bam_read(sup, what = c("rname", "pos", "cigar", "strand"),
                   as = "GAlignmentPairs", BPPARAM = NULL)
  )
  expect_identical(gp, GenomicAlignments::readGAlignmentPairs(sup))
  expect_identical(BiocGenerics::start(GenomicAlignments::first(gp)), 1000L)
})


# --- Summaries with unmapped reads, and unsupported param settings ------------

test_that("mapq_dist and fragment_sizes skip unmapped reads, like scanBam", {
  um <- .make_test_bam(c(
    sprintf("u%d\t133\tchr1\t50\t0\t*\t=\t50\t0\t*\t*", 1:11),
    sprintf("u%d\t73\tchr1\t50\t60\t20M\t=\t50\t0\t*\t*", 1:11),
    sprintf("z%d\t0\tchr1\t60\t0\t10M\t*\t0\t0\t*\t*", 1:7)
  ))
  sb <- Rsamtools::scanBam(um, param = Rsamtools::ScanBamParam(what = c("mapq", "isize")))[[1]]

  mq <- mapq_dist(um, BPPARAM = NULL)
  ref <- table(sb$mapq)
  expect_identical(mq$mapq, as.integer(names(ref)))
  expect_identical(mq$count, as.numeric(ref))

  fs <- fragment_sizes(um, BPPARAM = NULL)
  ref <- table(abs(sb$isize))
  expect_identical(fs$fragment_size, as.integer(names(ref)))
  expect_identical(fs$count, as.numeric(ref))
})

test_that("unsupported ScanBamParam settings give a warning", {
  bam <- ompBAM::example_BAM("Unsorted")
  region <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges("1:1-1000000"))
  expect_warning(mapq_dist(bam, param = region, BPPARAM = NULL), "which")
  expect_warning(fragment_sizes(bam, param = region, BPPARAM = NULL), "which")
  expect_warning(bam_coverage(bam, param = region, BPPARAM = NULL), "which")
  expect_no_warning(bam_count(bam, param = region, BPPARAM = NULL))
  expect_warning(
    bam_read(bam, what = "qname", BPPARAM = NULL,
             param = Rsamtools::ScanBamParam(tagFilter = list(NM = 0L))),
    "tagFilter"
  )
})

test_that("bam_coverage_bigwig rejects repeated output paths", {
  bam <- ompBAM::example_BAM("Unsorted")
  out <- tempfile(fileext = ".bw")
  dup <- file.path(dirname(out), ".", basename(out))
  expect_error(
    bam_coverage_bigwig(c(bam, bam), c(out, dup), BPPARAM = NULL),
    "must be unique"
  )
})


# --- which: overlap semantics --------------------------------------------------

test_that("which keeps alignments that overlap a range, as Rsamtools does", {
  src <- ompBAM::example_BAM("Unsorted")
  s <- Rsamtools::sortBam(src, tempfile("sorted"))
  Rsamtools::indexBam(s)
  on.exit(unlink(c(s, paste0(s, ".bai"))), add = TRUE)

  # Two reads start before this range and extend into it.
  one <- GenomicRanges::GRanges("1", IRanges::IRanges(156221565, 156221585))
  ref <- Rsamtools::scanBam(s, param = Rsamtools::ScanBamParam(which = one, what = "qname"))[[1]]$qname
  got <- bam_read(s, param = Rsamtools::ScanBamParam(which = one), what = "qname",
                  as = "data.frame", BPPARAM = NULL)$qname
  expect_length(ref, 2L)
  expect_identical(sort(got), sort(ref))

  # Disjoint windows: the same alignments as readGAlignments() (which repeats
  # an alignment once per window it overlaps; BamScaleR returns it once).
  ga_all <- GenomicAlignments::readGAlignments(s)
  st <- sort(BiocGenerics::start(ga_all[as.character(GenomicAlignments::seqnames(ga_all)) == "1"]))
  st <- st[seq(1L, length(st), by = 25L)] - 300L
  st <- st[c(TRUE, diff(st) > 1000L)]
  gr <- GenomicRanges::GRanges("1", IRanges::IRanges(st, width = 1000L))
  ref <- GenomicAlignments::readGAlignments(s, param = Rsamtools::ScanBamParam(which = gr, what = "qname"))
  got <- bam_read(s, param = Rsamtools::ScanBamParam(which = gr),
                  what = c("rname", "pos", "cigar", "strand", "qname"),
                  as = "GAlignments", BPPARAM = NULL)
  key <- function(g) {
    sort(paste(BiocGenerics::start(g), GenomicAlignments::cigar(g), S4Vectors::mcols(g)[["qname"]]))
  }
  expect_identical(key(got), unique(key(ref)))
  expect_identical(
    sum(bam_count(s, param = Rsamtools::ScanBamParam(which = gr), BPPARAM = NULL)$count),
    as.numeric(length(got))
  )
})

test_that("an empty which means no restriction, as in Rsamtools", {
  bam <- ompBAM::example_BAM("Unsorted")
  all <- bam_count(bam, BPPARAM = NULL)
  expect_identical(
    bam_count(bam, param = Rsamtools::ScanBamParam(which = IRanges::IRangesList(`1` = IRanges::IRanges())),
              BPPARAM = NULL),
    all
  )
  expect_identical(
    bam_count(bam, param = list(which = data.frame(seqname = character(), start = integer(),
                                                   end = integer())),
              BPPARAM = NULL),
    all
  )
})


# --- Empty GAlignmentPairs, bigWig output paths, list param --------------------

test_that("GAlignmentPairs output is an empty object when nothing can be paired", {
  bam <- ompBAM::example_BAM("Unsorted")
  empty_region <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges("2", IRanges::IRanges(90000, 90100)))
  gp <- bam_read(bam, param = empty_region, what = c("rname", "pos", "cigar", "strand", "mapq"),
                 as = "GAlignmentPairs", BPPARAM = NULL)
  expect_s4_class(gp, "GAlignmentPairs")
  expect_length(gp, 0L)
  expect_length(GenomeInfoDb::seqlevels(gp), 25L)

  single_end <- .make_test_bam(character())  # filler reads only, all single-end
  gp <- bam_read(single_end, as = "GAlignmentPairs", BPPARAM = NULL)
  expect_s4_class(gp, "GAlignmentPairs")
  expect_length(gp, 0L)
})

test_that("bam_coverage_bigwig checks output paths before writing", {
  bam <- ompBAM::example_BAM("Unsorted")
  expect_error(
    bam_coverage_bigwig(bam, file.path(tempfile(), "no_such_dir", "out.bw"), BPPARAM = NULL),
    "cannot write"
  )
  victim <- tempfile(fileext = ".bam")
  on.exit(unlink(victim), add = TRUE)
  file.copy(bam, victim)
  expect_error(bam_coverage_bigwig(victim, victim, BPPARAM = NULL), "input BAM")
  expect_identical(file.size(victim), file.size(bam))
})

test_that("list-form param warns about settings it does not know", {
  bam <- ompBAM::example_BAM("Unsorted")
  expect_warning(bam_count(bam, param = list(mapqfilter = 30), BPPARAM = NULL), "mapqfilter")
  expect_no_warning(bam_count(bam, param = list(mapqFilter = 30), BPPARAM = NULL))
})



# --- Round-4 regressions and edge cases -----------------------------------------

test_that("a capped MulticoreParam works inside forked workers and is not modified", {
  skip_on_os("windows")
  bp <- tryCatch(suppressWarnings(BiocParallel::MulticoreParam(3L)), error = function(e) NULL)
  if (is.null(bp) || BiocParallel::bpnworkers(bp) < 3L) {
    skip("fewer than 3 BiocParallel workers allowed in this environment")
  }
  files <- stats::setNames(as.character(1:3), letters[1:3])
  out <- BamScaleR:::.bamscale_apply_files(
    files, function(x) c(x = x, pid = as.character(Sys.getpid())),
    BPPARAM = bp, bp_workers = 2L
  )
  expect_identical(names(out), letters[1:3])
  expect_lte(length(unique(vapply(out, `[[`, "", "pid"))), 2L)
  expect_identical(BiocParallel::bpnworkers(bp), 3L)
  expect_false(BiocParallel::bpisup(bp))
  # A capped call inside forked workers (it failed when the cap started a
  # socket-based fork cluster).
  nested <- parallel::mclapply(1:2, function(j) {
    unname(unlist(BamScaleR:::.bamscale_apply_files(files, identity, BPPARAM = bp, bp_workers = 2L)))
  }, mc.cores = 2L)
  expect_identical(nested, list(as.character(1:3), as.character(1:3)))
})

test_that("which matches Rsamtools at range boundaries, deletions and splices", {
  # Records in coordinate order, so that the BAM can be indexed.
  bam <- .make_test_bam(c(
    "before\t0\tchr1\t80\t60\t10M\t*\t0\t0\t*\t*",
    "end_at_start\t0\tchr1\t91\t60\t10M\t*\t0\t0\t*\t*",
    "start_at_end\t0\tchr1\t120\t60\t10M\t*\t0\t0\t*\t*",
    "after\t0\tchr1\t121\t60\t10M\t*\t0\t0\t*\t*",
    "deletion\t0\tchr1\t285\t60\t5M20D5M\t*\t0\t0\t*\t*",
    "splice\t0\tchr1\t500\t60\t5M100N5M\t*\t0\t0\t*\t*"
  ), seqlen = 1000L)
  Rsamtools::indexBam(bam)
  on.exit(unlink(c(bam, paste0(bam, ".bai"))), add = TRUE)
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 300, 550), c(120, 302, 560)))
  ref <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(which = gr, what = "qname"))
  got <- bam_read(bam, param = Rsamtools::ScanBamParam(which = gr), what = "qname",
                  as = "scanBam", BPPARAM = NULL)
  expect_identical(
    lapply(got, function(b) sort(b$qname)),
    stats::setNames(lapply(ref, function(b) sort(b$qname)), names(got))
  )
  expect_identical(sort(got[[1]]$qname), c("end_at_start", "start_at_end"))
})

test_that("which: zero-width ranges, missing sequences, NA ranges and repeated labels", {
  bam <- ompBAM::example_BAM("Unsorted")
  zero <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges("1", IRanges::IRanges(630300, width = 0)))
  expect_identical(sum(bam_count(bam, param = zero, BPPARAM = NULL)$count), 0)
  expect_identical(
    nrow(bam_read(bam, param = zero, what = "qname", as = "data.frame", BPPARAM = NULL)), 0L
  )
  sb <- bam_read(bam, param = zero, what = "qname", as = "scanBam", BPPARAM = NULL)
  expect_identical(names(sb), "1:630300-630299")
  expect_length(sb[[1]]$qname, 0L)
  expect_error(
    bam_count(bam, param = list(which = data.frame(seqname = "1", start = 300L, end = 100L)),
              BPPARAM = NULL),
    "end is less than start - 1"
  )
  expect_error(
    bam_count(bam, param = Rsamtools::ScanBamParam(which = GenomicRanges::GRanges("chr1:1-100")),
              BPPARAM = NULL),
    "not in BAM header"
  )
  expect_error(
    bam_count(bam, param = list(which = data.frame(seqname = "1", start = NA, end = 1e6)),
              BPPARAM = NULL),
    "missing"
  )
  two <- data.frame(seqname = c("1", "10"), start = c(1L, 1L), end = c(1000000L, 3000000L),
                    label = c("x", "x"))
  got <- bam_read(bam, param = list(which = two), what = "qname", as = "scanBam", BPPARAM = NULL)
  expect_identical(names(got), c("x", "x.1"))
  one <- function(i) {
    bam_read(bam, param = list(which = two[i, c("seqname", "start", "end")]), what = "qname",
             as = "scanBam", BPPARAM = NULL)[[1]]$qname
  }
  expect_identical(got[[1]]$qname, one(1))
  expect_identical(got[[2]]$qname, one(2))
  # with.which_label adds a column only when which is set
  expect_false("which_label" %in% names(
    bam_read(bam, what = "qname", with.which_label = TRUE, as = "data.frame", BPPARAM = NULL)
  ))
})

test_that("list-form param applies flag and matches names exactly", {
  bam <- ompBAM::example_BAM("Unsorted")
  fl <- Rsamtools::scanBamFlag(isMinusStrand = TRUE)
  expect_identical(
    bam_count(bam, param = list(flag = fl), BPPARAM = NULL),
    bam_count(bam, param = Rsamtools::ScanBamParam(flag = fl), BPPARAM = NULL)
  )
  # `flag` is combined with flag_set / flag_unset given alongside it
  both <- list(flag = fl, flag_set = 16L, flag_unset = 16L)
  expect_identical(sum(bam_count(bam, param = both, BPPARAM = NULL)$count), 0)
  expect_identical(
    bam_count(bam, param = list(flag = fl, flag_set = 64L), BPPARAM = NULL),
    bam_count(bam, param = Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE, isFirstMateRead = TRUE)
    ), BPPARAM = NULL)
  )
  expect_warning(
    x <- bam_read(bam, what = "qname", param = list(tagFilter = list(NM = 0L)),
                  as = "data.frame", BPPARAM = NULL),
    "tagFilter"
  )
  expect_identical(names(x), "qname")
})

test_that("use.names = TRUE names GAlignments even without qname in what", {
  bam <- ompBAM::example_BAM("Unsorted")
  core <- c("rname", "pos", "cigar", "strand")
  ref <- GenomicAlignments::readGAlignments(bam, use.names = TRUE)
  ga <- bam_read(bam, what = core, as = "GAlignments", use.names = TRUE, BPPARAM = NULL)
  expect_identical(names(ga), names(ref))
  expect_identical(ncol(S4Vectors::mcols(ga)), 0L)
  old <- options(BamScaleR.ga_fastpath = FALSE)
  on.exit(options(old), add = TRUE)
  ga2 <- bam_read(bam, what = core, as = "GAlignments", use.names = TRUE, BPPARAM = NULL)
  expect_identical(names(ga2), names(ref))
})

test_that("relative input paths are made absolute for every input type", {
  skip_on_os("windows")
  d <- tempfile("inputs")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  file.copy(ompBAM::example_BAM("Unsorted"), file.path(d, "real.bam"))
  file.symlink("real.bam", file.path(d, "link.bam"))
  old <- setwd(d)
  on.exit(setwd(old), add = TRUE, after = FALSE)
  want <- file.path(getwd(), "link.bam")
  expect_identical(unname(BamScaleR:::.bamscale_normalize_files("link.bam")), want)
  expect_identical(unname(BamScaleR:::.bamscale_normalize_files("./link.bam")), want)
  expect_identical(
    unname(BamScaleR:::.bamscale_normalize_files(Rsamtools::BamFile("link.bam"))), want
  )
  expect_identical(
    unname(BamScaleR:::.bamscale_normalize_files(Rsamtools::BamFileList("link.bam"))), want
  )
  expect_identical(names(BamScaleR:::.bamscale_normalize_files("link.bam")), "link.bam")
})

test_that("bam_coverage_bigwig rejects a directory and a link to another outfile", {
  bam <- ompBAM::example_BAM("Unsorted")
  expect_error(bam_coverage_bigwig(bam, tempdir(), BPPARAM = NULL), "directory")
  not_yet <- file.path(tempfile("tracks"), "")  # ends in "/", does not exist
  expect_error(bam_coverage_bigwig(bam, not_yet, BPPARAM = NULL), "directory")
  expect_false(file.exists(sub("/$", "", not_yet)))
  skip_on_os("windows")
  real <- tempfile(fileext = ".bw")
  link <- tempfile(fileext = ".bw")
  file.symlink(real, link)  # dangling until real is written
  on.exit(unlink(c(real, link)), add = TRUE)
  expect_error(bam_coverage_bigwig(c(bam, bam), c(real, link), BPPARAM = NULL), "unique")
})
