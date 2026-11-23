#!/usr/bin/env Rscript
# ============================================================================
# make_benchmark_extdata.R
#
# Regenerates the four benchmark summary files in inst/extdata/ that the
# "benchmark-results" vignette reads:
#
#   benchmark_summary.csv          per-case timings used by every table, figure
#                                  and inline number in the vignette
#   benchmark_encode_manifest.csv  the 25 public ENCODE ATAC-seq BAMs used
#   benchmark_run_info.csv         host, software versions and git commit
#   benchmark_correctness.csv      output-equivalence checks behind the timings
#
# Usage (from the package root):
#   Rscript inst/scripts/make_benchmark_extdata.R <RUN> [OUTDIR]
#     <RUN>    a finished benchmark run directory (see below); the published
#              files come from benchmark_results/final_20261008
#     OUTDIR   where to write the CSVs (default: inst/extdata)
# Re-running it on the same <RUN> reproduces the shipped files byte for byte.
# Only base R is used.
#
# ---------------------------------------------------------------------------
# Provenance of the published run (final_20261008)
# ---------------------------------------------------------------------------
# Code     BamScaleR 0.99.1, git commit 30f6194729b5d26b7de95e1d8b451480212e6183
#          (no tracked-file changes at run start); run started 2026-10-08.
# Host     one node, two Intel Xeon Gold 6252 sockets (2.10 GHz; 48 physical /
#          96 logical cores), 723.6 GiB RAM, Ubuntu 24.04.5 (Linux 6.8.0-139),
#          R 4.6.1, Rsamtools 2.28.0, GenomicAlignments 1.48.0, rtracklayer
#          1.72.0, ompBAM 1.16.0, BiocParallel 1.46.0, ATACseqQC 1.36.1.
#          OMP_PROC_BIND unset (threads not pinned). Before every stage each
#          input BAM was read once (cat > /dev/null) so timings are warm-cache.
#
# Step 1. Corpus selection -- inst/benchmarks/select_encode_atac.R
#   Rscript inst/benchmarks/select_encode_atac.R --out=<DATA>
#   Queries the ENCODE portal for released GRCh38 paired-end ATAC-seq alignment
#   BAMs and applies the script's default criteria: one "single" file of
#   5-9 GB whose ENCODE QC total_reads is closest above 150 million, plus 24
#   "multi" files of 2.5-4.5 GB, one per experiment, largest first (excluding
#   ENCFF646NWY and ENCFF415FEC). It writes <DATA>/encode_atac_manifest.csv
#   (accession, experiment, biosample, size, md5, QC read count, URL). The
#   portal changes over time, so re-running the query may select differently;
#   the manifest, shipped here as benchmark_encode_manifest.csv, pins the
#   corpus by accession and md5.
#
# Step 2. Download -- inst/benchmarks/download_atac_data.R
#   Rscript inst/benchmarks/download_atac_data.R \
#     --manifest=<DATA>/encode_atac_manifest.csv --role=single \
#     --dest=<DATA>/single
#   Rscript inst/benchmarks/download_atac_data.R \
#     --manifest=<DATA>/encode_atac_manifest.csv --role=multi \
#     --dest=<DATA>/multi24
#   Each BAM is md5-verified against the manifest, indexed, and checked to be
#   paired-end.
#
# Step 3. Benchmarks. Common arguments:
#   T="--threads=1,2,4,8,12,24,36,48,72 --budget-threads=72 \
#      --max-threads=72 --bp-backend=snow"
#
#   3a. Read patterns, single file -- inst/benchmarks/run_server_benchmark.R
#     Rscript inst/benchmarks/run_server_benchmark.R --profile=full \
#       --bam-dir=<DATA>/single --recursive=false \
#       --single-workloads=step1,galignments,seqqual --include-multi=false \
#       --include-galignments=true --include-rsamtools=true \
#       --include-seqqual=true --seqqual-compact=true --n-files=1 \
#       --allow-repeat-files=true $T --iterations=5 --outdir=<RUN>/read_single
#   3b. Read patterns, 24 files (merged into merged_read, but its rows are
#       not used here)
#     Rscript inst/benchmarks/run_server_benchmark.R --profile=full \
#       --bam-dir=<DATA>/multi24 --recursive=false --include-single=false \
#       --multi-workloads=step1 --include-galignments=false \
#       --include-rsamtools=true --include-seqqual=false \
#       --include-auto-threads=false --n-files=24 --workers=12,24 $T \
#       --iterations=2 --outdir=<RUN>/read_multi
#   3c. Workflows, single file -- inst/benchmarks/run_workflow_benchmark.R
#     Rscript inst/benchmarks/run_workflow_benchmark.R \
#       --include-coverage=true --include-atacqc=true \
#       --bam-dir=<DATA>/single --atac-bam-dir=<DATA>/single \
#       --include-multi=false --write-bigwig=true --atac-anchor=true \
#       $T --iterations=5 --outdir=<RUN>/wf_single
#   3d. Workflows, 24 files (same script)
#     Rscript inst/benchmarks/run_workflow_benchmark.R \
#       --include-coverage=true --include-atacqc=true \
#       --bam-dir=<DATA>/multi24 --atac-bam-dir=<DATA>/multi24 \
#       --include-single=false --n-files=24 --write-bigwig=true \
#       --multi-write-bigwig=false --workers=4,12,24 $T --iterations=2 \
#       --outdir=<RUN>/wf_multi
#   The run also holds in-reader aggregation stages (api_single, api_multi)
#   and a samtools reference (samtools_reference.sh); none of their output
#   is used here.
#
# Step 4. Merge -- inst/benchmarks/merge_runs.R
#   Each stage writes a time-stamped subdirectory (run_* or wf_run_*):
#   Rscript inst/benchmarks/merge_runs.R --out=<RUN>/merged_read \
#     --dirs=<read_single/run_*>,<read_multi/run_*>
#   Rscript inst/benchmarks/merge_runs.R --out=<RUN>/merged_cov \
#     --dirs=<wf_single/wf_run_*>,<wf_multi/wf_run_*> --workloads=coverage
#   Rscript inst/benchmarks/merge_runs.R --out=<RUN>/merged_atac \
#     --dirs=<wf_single/wf_run_*>,<wf_multi/wf_run_*> --workloads=atacqc
#   cp <DATA>/encode_atac_manifest.csv <RUN>/
#   <RUN>/RUN_NOTES.md records the git commit, package version and host at
#   run start.
#
# Step 5. This script reads, from <RUN>:
#   RUN_NOTES.md, encode_atac_manifest.csv,
#   merged_read/{summary,correctness_preflight}.csv,
#   merged_cov/{summary,correctness_workflow,host_info,pkg_versions}.csv,
#   merged_atac/{summary,correctness_workflow,atacseqqc_anchor}.csv,
#   merged_read/sessionInfo.txt
#
# ---------------------------------------------------------------------------
# Columns of benchmark_summary.csv (one row per benchmarked configuration)
# ---------------------------------------------------------------------------
#   benchmark         "read" (read-pattern benchmark) or "workflow"
#   workload          step1 (core fields qname/flag/rname/pos/mapq/cigar),
#                     galignments (GAlignments object), seqqual (sequence +
#                     base quality), coverage (reads -> coverage() -> bigWig),
#                     atacqc (ATAC fragment-size table)
#   scenario          "single" (one BAM) or "multi" (24 BAMs)
#   comparison_track  read rows: "fair" (BamScaleR returns the same object as
#                     the comparator); workflow rows: "single-file" or
#                     "matched-cores" (both arms use the same number of cores)
#   method_family     BamScaleR, Rsamtools or GenomicAlignments
#   method            the harness's label for the arm
#   threads           OpenMP threads per file (1 for the standard readers)
#   workers           BiocParallel workers (files processed concurrently)
#   cores_used        threads x workers
#   n_files, n_records  BAM files and alignments processed
#   iterations        timed iterations (after one untimed warm-up)
#   median_s          median wall-clock seconds of the whole timed case: the
#                     read call for read rows, the end-to-end workflow for
#                     workflow rows
#   read_s, compute_s, write_s
#                     single-file workflow rows only: median seconds of each
#                     phase (taken per phase, so they need not sum to median_s)
# Kept rows: read benchmark, status ok, single file, "fair" track (the
# compact-seqqual "optimized" track and the 24-file read rows are dropped);
# workflow benchmark, status ok, "single-file" and "matched-cores" tracks (the
# 24-file "full-budget" sweep is dropped). Times are rounded to 1e-4 s, which
# only removes floating-point noise (the timer resolution is 1 ms).
#
# ---------------------------------------------------------------------------
# Columns of the other files
# ---------------------------------------------------------------------------
# benchmark_encode_manifest.csv (one row per BAM, in manifest order)
#   role              "single" (target of every single-file benchmark) or
#                     "multi" (the 24-file cohort)
#   accession, experiment   ENCODE file (ENCFF...) and experiment (ENCSR...)
#   biosample         ENCODE biosample term name
#   file_size_bytes, md5sum   as reported by the ENCODE portal
#   total_reads       ENCODE alignment QC total_reads
#   url               download URL
# benchmark_run_info.csv (field, value): run id and start date, BamScaleR
#   version and git commit (from RUN_NOTES.md; git_diff_stat "none" means
#   `git diff --stat` was empty at run start), CPU model, logical cores, RAM
#   (GiB), OS, R version, single-file thread grid, multi-file core budget, and
#   the versions of the comparator packages.
# benchmark_correctness.csv (one row per equivalence check)
#   check             read_preflight (BamScaleR vs the standard reader on a
#                     1 kb window, per read pattern), workflow_equivalence
#                     (workflow output identical() to the standard pipeline's),
#                     atacseqqc_anchor (fragment-size table vs the installed
#                     ATACseqQC::fragSizeDist over the whole single file)
#   workload, comparator, passed, detail   as recorded by the harness
#   n_compared        records compared, where the harness reports it
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
    stop("usage: Rscript make_benchmark_extdata.R <RUN> [OUTDIR]",
         call. = FALSE)
}
run <- normalizePath(args[[1L]], mustWork = TRUE)
outdir <- if (length(args) == 2L) args[[2L]] else file.path("inst", "extdata")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rd <- function(...) {
    utils::read.csv(file.path(run, ...), stringsAsFactors = FALSE)
}
wr <- function(x, name) {
    rownames(x) <- NULL
    utils::write.csv(x, file.path(outdir, name), row.names = FALSE,
                     eol = "\n", fileEncoding = "UTF-8")
    message(sprintf("wrote %-30s %3d rows", name, nrow(x)))
}
sort_rows <- function(x, by) {
    x[do.call(order, c(unname(as.list(x[by])), method = "radix")), ,
      drop = FALSE]
}
tm <- function(x) round(x, 4L)

families <- c("BamScaleR", "Rsamtools", "GenomicAlignments")

# ---- benchmark_summary.csv -------------------------------------------------
rs <- rd("merged_read", "summary.csv")
rs <- rs[rs$status == "ok" & rs$scenario == "single" &
             rs$comparison_track == "fair", ]
read_rows <- data.frame(
    benchmark = "read", workload = rs$workload, scenario = rs$scenario,
    comparison_track = rs$comparison_track,
    method_family = rs$method_family, method = rs$method,
    threads = rs$threads_effective, workers = rs$bp_workers_effective,
    cores_used = rs$total_threads_effective,
    n_files = rs$n_files, n_records = rs$n_records,
    iterations = rs$iterations_completed, median_s = tm(rs$median_s),
    read_s = NA_real_, compute_s = NA_real_, write_s = NA_real_,
    stringsAsFactors = FALSE)

ws <- rbind(rd("merged_cov", "summary.csv"), rd("merged_atac", "summary.csv"))
ws <- ws[ws$status == "ok" &
             ws$comparison_track %in% c("single-file", "matched-cores"), ]
wf_rows <- data.frame(
    benchmark = "workflow", workload = ws$workload, scenario = ws$scenario,
    comparison_track = ws$comparison_track,
    method_family = ws$method_family, method = ws$method,
    threads = ws$threads_effective, workers = ws$bp_workers_effective,
    cores_used = ws$cores_used,
    n_files = ws$n_files, n_records = ws$n_records,
    iterations = ws$iterations_completed, median_s = tm(ws$total_s),
    read_s = tm(ws$read_s), compute_s = tm(ws$compute_s),
    write_s = tm(ws$write_s),
    stringsAsFactors = FALSE)

summary_tbl <- rbind(read_rows, wf_rows)
stopifnot(
    all(summary_tbl$method_family %in% families),
    !anyNA(summary_tbl$median_s),
    all(summary_tbl$cores_used == summary_tbl$threads * summary_tbl$workers),
    # every single-file case timed the same BAM
    length(unique(
        summary_tbl$n_records[summary_tbl$scenario == "single"])) == 1L,
    # single-file workflow rows carry the full phase decomposition
    !anyNA(summary_tbl[summary_tbl$benchmark == "workflow" &
                           summary_tbl$scenario == "single",
                       c("read_s", "compute_s", "write_s")])
)
summary_tbl <- sort_rows(summary_tbl,
    c("benchmark", "workload", "comparison_track", "method_family",
      "threads", "workers"))
wr(summary_tbl, "benchmark_summary.csv")

# ---- benchmark_encode_manifest.csv -----------------------------------------
man <- rd("encode_atac_manifest.csv")
manifest_tbl <- data.frame(
    role = man$role, accession = man$accession, experiment = man$experiment,
    biosample = man$biosample, file_size_bytes = man$file_size,
    md5sum = man$md5sum, total_reads = man$total_reads, url = man$url,
    stringsAsFactors = FALSE)
stopifnot(
    sum(manifest_tbl$role == "single") == 1L,
    !anyDuplicated(manifest_tbl$accession),
    !anyDuplicated(manifest_tbl$experiment),
    all(grepl("^[0-9a-f]{32}$", manifest_tbl$md5sum)),
    # the single-file benchmark target is the manifest's "single" BAM
    manifest_tbl$total_reads[manifest_tbl$role == "single"] ==
        unique(summary_tbl$n_records[summary_tbl$scenario == "single"])
)
wr(manifest_tbl, "benchmark_encode_manifest.csv")

# ---- benchmark_run_info.csv ------------------------------------------------
notes <- readLines(file.path(run, "RUN_NOTES.md"), encoding = "UTF-8",
                   warn = FALSE)
note <- function(key) {
    hit <- grep(paste0("^- ", key, ":"), notes, value = TRUE)
    if (length(hit) != 1L) stop("RUN_NOTES.md: no unique '", key, "' line")
    trimws(sub("^[^:]*:", "", hit))
}
run_date <- regmatches(notes[1L], regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}",
                                          notes[1L]))
git_diff <- note("working tree")
host <- rd("merged_cov", "host_info.csv")
pkgs <- rd("merged_cov", "pkg_versions.csv")
pkg_version <- function(p) pkgs$version[pkgs$package == p]
sess <- readLines(file.path(run, "merged_read", "sessionInfo.txt"),
                  warn = FALSE)
os <- trimws(sub("^Running under:", "", grep("^Running under:", sess,
                                             value = TRUE)))
stopifnot(nrow(host) == 1L, length(run_date) == 1L,
          identical(note("BamScaleR"), pkg_version("BamScaleR")))

info <- c(
    run_id            = basename(run),
    run_date          = run_date,
    bamscaler_version = note("BamScaleR"),
    git_head          = note("git HEAD"),
    git_diff_stat     = if (nzchar(git_diff)) git_diff else "none",
    cpu_model         = host$cpu_model,
    logical_cores     = as.character(host$logical_cores),
    total_mem_gib     = as.character(host$total_mem_gb),
    os                = paste0(os, " (Linux ", host$release, ")"),
    r_version         = host$r_version,
    thread_grid       = sub(" .*$", "", sub("^threads=", "", note("grids"))),
    core_budget       = sub("^.*budget=", "", note("grids")),
    Rsamtools_version         = pkg_version("Rsamtools"),
    GenomicAlignments_version = pkg_version("GenomicAlignments"),
    rtracklayer_version       = pkg_version("rtracklayer"),
    ompBAM_version            = pkg_version("ompBAM"),
    BiocParallel_version      = pkg_version("BiocParallel"),
    ATACseqQC_version         = pkg_version("ATACseqQC"))
stopifnot(all(lengths(info) == 1L), !anyNA(info), all(nzchar(info)))
wr(data.frame(field = names(info), value = unname(info),
              stringsAsFactors = FALSE), "benchmark_run_info.csv")

# ---- benchmark_correctness.csv ---------------------------------------------
num_after <- function(x, key) {
    m <- regmatches(x, regexpr(paste0(key, "=[0-9]+"), x))
    out <- rep(NA_real_, length(x))
    out[regexpr(paste0(key, "=[0-9]+"), x) > 0] <- as.numeric(sub(".*=", "", m))
    out
}
pf <- rd("merged_read", "correctness_preflight.csv")
stopifnot(all(pf$bam_n == pf$ref_n))
pre_rows <- data.frame(
    check = "read_preflight", workload = pf$workload,
    comparator = pf$comparator, passed = pf$passed, n_compared = pf$bam_n,
    detail = sprintf("%s:%d-%d window; columns %s; identical=%s",
                     pf$seqname, pf$start, pf$end, pf$columns_compared,
                     pf$identical_data),
    stringsAsFactors = FALSE)

cw <- rbind(rd("merged_cov", "correctness_workflow.csv"),
            rd("merged_atac", "correctness_workflow.csv"))
wf_chk <- data.frame(
    check = "workflow_equivalence", workload = cw$workload,
    comparator = cw$comparator, passed = cw$passed,
    n_compared = num_after(cw$detail, "std n"), detail = cw$detail,
    stringsAsFactors = FALSE)

an <- rd("merged_atac", "atacseqqc_anchor.csv")
stopifnot(identical(num_after(an$detail, "ours"), num_after(an$detail, "pkg")))
anchor_rows <- data.frame(
    check = "atacseqqc_anchor", workload = "atacqc", comparator = an$check,
    passed = an$passed, n_compared = num_after(an$detail, "pkg"),
    detail = an$detail, stringsAsFactors = FALSE)

correctness_tbl <- rbind(pre_rows, wf_chk, anchor_rows)
stopifnot(is.logical(correctness_tbl$passed), !anyNA(correctness_tbl$passed))
wr(correctness_tbl, "benchmark_correctness.csv")
