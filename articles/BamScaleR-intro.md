# Getting started with BamScaleR

## Introduction

`BamScaleR` provides multithreaded sequential BAM traversal for
R/Bioconductor workflows. The package is designed for users who already
rely on `Rsamtools`, `GenomicAlignments`, and `BiocParallel`, but need
faster BAM parsing without moving to a separate non-Bioconductor
workflow model.

The package targets four common classes of BAM access:

- metadata-oriented scans for filtering, fragment summaries, and quality
  control
- generation of alignment objects for downstream `GenomicAlignments`
  workflows
- sequence and quality extraction, including an optimized compact mode
- summaries computed inside the multithreaded reader (fragment-size and
  mapping-quality distributions, per-base coverage, bigWig tracks), so
  that no per-read R objects are created at all

The motivation for inclusion in Bioconductor is therefore
straightforward: `BamScaleR` accelerates a core data-access bottleneck
while preserving familiar Bioconductor inputs and output types.

`BamScaleR` is built on the `ompBAM` C++ library, which supplies
parallel BAM decompression. The section [BamScaleR and
ompBAM](#bamscaler-and-ompbam) below explains what each package
provides.

## Installation

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("BamScaleR")
```

### Multithreading and OpenMP

`BamScaleR` reads each BAM with several threads (the `threads` argument)
using OpenMP. OpenMP is optional. Without it, every function works and
returns the same results, but each BAM is read on a single thread.
Processing several files in parallel through `BPPARAM` (see below) is
unaffected. If `threads` greater than 1 is requested on a build without
OpenMP, `BamScaleR` warns once per session.

- **Linux and Windows.** The standard R toolchains (GCC, and Rtools on
  Windows) support OpenMP. Nothing extra is needed.
- **macOS.** Apple’s `clang` can compile OpenMP code, but Apple does not
  ship the OpenMP header (`omp.h`) or runtime library (`libomp`). CRAN’s
  R for macOS ships a `libomp` in `R.home("lib")`, but no header. When
  `BamScaleR` is installed from source, its `configure` script tries
  several compiler and linker settings. For each one it builds a small
  OpenMP test library and loads it into R, and it uses the first setting
  that loads. With CRAN’s R, the runtime that gets linked is always R’s
  own `libomp`, so only one OpenMP runtime is loaded into the R session.
  The chosen flags are printed in the installation log. If no setting
  works, the package is built without OpenMP and the log says why. (If
  your own compiler flags switch OpenMP on for every package, for
  example `CPPFLAGS += -Xclang -fopenmp` in `~/.R/Makevars`, a build
  without OpenMP is not possible; `configure` then stops and says so.)

To enable multithreading on macOS, install the OpenMP files that match
your Xcode version from <https://mac.r-project.org/openmp/>. Check the
Apple clang version with `clang --version` in a terminal. For example,
for Xcode 16.0 to 16.2 (Apple clang 1600):

``` sh
curl -O https://mac.r-project.org/openmp/openmp-17.0.6-darwin20-Release.tar.gz
sudo tar fvxz openmp-17.0.6-darwin20-Release.tar.gz -C /
```

Use `openmp-16.0.4-darwin20-Release.tar.gz` for Xcode 15 (Apple clang
1500), or `openmp-19.1.5-darwin20-Release.tar.gz` for Xcode 16.3 or
later (Apple clang 1700). The tarball installs `omp.h` under
`/usr/local/include` and `libomp.dylib` under `/usr/local/lib`. A
Homebrew installation (`brew install libomp`) is also detected.

Then reinstall `BamScaleR` from source. A plain
[`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
would skip a package that is already installed, or install the same
binary again, so ask for a source build explicitly (this needs the Xcode
Command Line Tools):

``` r

BiocManager::install("BamScaleR", type = "source", force = TRUE)
```

In a source checkout, rebuild with `R CMD INSTALL --preclean .` or
`devtools::load_all(recompile = TRUE)`, so that objects compiled without
OpenMP are not reused.

Apple clang 1700 (Xcode 16.3 and later) needs an OpenMP runtime from
LLVM 19 or later. If the `libomp` bundled with your R is older,
`configure` rejects it, the package is built without OpenMP, and the
installation log explains why.

To bypass the detection, set both flags yourself before installing, for
example:

``` r

Sys.setenv(
    BAMSCALER_OPENMP_CXXFLAGS = "-Xclang -fopenmp",
    BAMSCALER_OPENMP_LIBS = "-lomp"
)
BiocManager::install("BamScaleR", type = "source", force = TRUE)
```

Setting both variables to `""` builds without OpenMP.

## Basic usage

``` r

library(BamScaleR)

bam <- ompBAM::example_BAM("Unsorted")
```

The example file shipped with `ompBAM` is a small paired-end RNA-seq BAM
(10,000 alignments to the human genome, aligned with STAR).

### Metadata-oriented BAM access

The most common BAM-access pattern is extraction of alignment metadata
such as read name, flag, reference name, position, mapping quality, and
CIGAR string.

``` r

x <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar"),
    as = "data.frame",
    threads = 2
)

head(x)
#>                                       qname flag rname      pos mapq
#> 1  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  163    19   572614  255
#> 2  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309   83    19   579499  255
#> 3 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309  163     6 44252112  255
#> 4 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309   83     6 44252124  255
#> 5 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309   99    12 46185884  255
#> 6 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309  147    12 46185968  255
#>           cigar
#> 1 1S88M6798N61M
#> 2      1S148M1S
#> 3      60S86M4S
#> 4      1S146M3S
#> 5        1S149M
#> 6      7S142M1S
```

### Filtering with `ScanBamParam`

`BamScaleR` accepts
[`Rsamtools::ScanBamParam`](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html),
allowing existing filtering logic to be reused directly.

``` r

param <- Rsamtools::ScanBamParam(
    what = c("qname", "flag", "mapq"),
    mapqFilter = 20L,
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
)

filtered <- bam_read(
    file = bam,
    param = param,
    as = "data.frame",
    threads = 2
)

head(filtered)
#>                                       qname flag mapq
#> 1  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  163  255
#> 2  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309   83  255
#> 3 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309  163  255
#> 4 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309   83  255
#> 5 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309   99  255
#> 6 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309  147  255
```

### Alignment-object output

When downstream workflows expect `GenomicAlignments` objects,
[`bam_read()`](https://cparsania.github.io/BamScaleR/reference/bam_read.md)
can return those directly.

``` r

ga <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos", "cigar", "strand"),
    as = "GAlignments",
    threads = 2
)

ga
#> GAlignments object with 10000 alignments and 2 metadata columns:
#>           seqnames strand              cigar    qwidth     start       end
#>              <Rle>  <Rle>        <character> <integer> <integer> <integer>
#>       [1]       19      +      1S88M6798N61M       150    572614    579560
#>       [2]       19      -           1S148M1S       150    579499    579646
#>       [3]        6      +           60S86M4S       150  44252112  44252197
#>       [4]        6      -           1S146M3S       150  44252124  44252269
#>       [5]       12      +             1S149M       150  46185884  46186032
#>       ...      ...    ...                ...       ...       ...       ...
#>    [9996]        2      -         16S122M12S       150 186680457 186680578
#>    [9997]       15      +           4S143M3S       150  90501177  90501319
#>    [9998]       15      -        26M2I56M66S       150  90501301  90501382
#>    [9999]        4      + 3M336N95M752N51M1S       150 121801488 121802724
#>   [10000]        4      -      1S117M121N32M       150 121802677 121802946
#>               width     njunc |                  qname      flag
#>           <integer> <integer> |            <character> <integer>
#>       [1]      6947         1 | ST-E00600:137:H77Y3C..       163
#>       [2]       148         0 | ST-E00600:137:H77Y3C..        83
#>       [3]        86         0 | ST-E00600:137:H77Y3C..       163
#>       [4]       146         0 | ST-E00600:137:H77Y3C..        83
#>       [5]       149         0 | ST-E00600:137:H77Y3C..        99
#>       ...       ...       ... .                    ...       ...
#>    [9996]       122         0 | ST-E00600:137:H77Y3C..        83
#>    [9997]       143         0 | ST-E00600:137:H77Y3C..       163
#>    [9998]        82         0 | ST-E00600:137:H77Y3C..        83
#>    [9999]      1237         2 | ST-E00600:137:H77Y3C..       163
#>   [10000]       270         1 | ST-E00600:137:H77Y3C..        83
#>   -------
#>   seqinfo: 25 sequences from an unspecified genome
```

### Sequence and quality extraction

For `seq` and `qual`, `BamScaleR` supports both a
compatibility-preserving mode and an optimized compact mode. The
compatible mode (the default) returns the sequences and qualities as
[`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
and
[`Biostrings::PhredQuality`](https://rdrr.io/pkg/Biostrings/man/XStringQuality-class.html)
objects, identical to those from
[`Rsamtools::scanBam()`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html).
A `DataFrame` holds such columns naturally; a base `data.frame` can
carry them, but cannot print them.

``` r

seqqual_compatible <- bam_read(
    file = bam,
    what = c("qname", "seq", "qual"),
    as = "DataFrame",
    seqqual_mode = "compatible",
    threads = 2
)

seqqual_compatible
#> DataFrame with 10000 rows and 3 columns
#>                        qname                     seq                    qual
#>                  <character>          <DNAStringSet>          <PhredQuality>
#> 1     ST-E00600:137:H77Y3C.. NGACCGGCGA...CCTGCTCCTT #AA-AAFJFJ...)7777<))--
#> 2     ST-E00600:137:H77Y3C.. TGCCGGCACA...GCCAGAAAAN -FAAA)-AJF...JJF-AFAAA#
#> 3     ST-E00600:137:H77Y3C.. NGGAGCTTCG...GATGGAGGAG #AA-AFFJFJ...7<--7---A7
#> 4     ST-E00600:137:H77Y3C.. TGGGAAGAGC...AGATGGTGAN AJA-FFJJJF...JJF<FFAAA#
#> 5     ST-E00600:137:H77Y3C.. NGCATCTCTA...ACATTATTCG #AAFFF-FJJ...JJJJ-<<FJJ
#> ...                      ...                     ...                     ...
#> 9996  ST-E00600:137:H77Y3C.. TTTGATGTTT...GTGCGGTTAT JJA-<FJJJJ...JJJ7JFFF<A
#> 9997  ST-E00600:137:H77Y3C.. CGTAAATTGT...TTTAGTTGAA A<A-AAFJ<7...)<-77<7-7<
#> 9998  ST-E00600:137:H77Y3C.. CAGTTTTCTG...TGAGCCTCTA AA<<<-<)A<...JJAFAFFFAA
#> 9999  ST-E00600:137:H77Y3C.. AAGCGGCTGG...CAATCTGGCA A<AAFFJJ7F...FJJ-<---AF
#> 10000 ST-E00600:137:H77Y3C.. TTCTTGGACA...TTGAATCGAC FJ7A<<<-7<...JJJJAFFFAA
```

Compact mode returns lower-level raw vectors for throughput-oriented
workflows. These values can be decoded back to ordinary strings
explicitly when needed.

``` r

seqqual_compact <- bam_read(
    file = bam,
    what = c("qname", "qwidth", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compact",
    threads = 2
)

head(seqqual_compact)
#>                                       qname qwidth          seq         qual
#> 1  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309    150 f4, 12, .... 02, 20, ....
#> 2  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309    150 84, 22, .... 0c, 25, ....
#> 3 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309    150 f4, 41, .... 02, 20, ....
#> 4 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309    150 84, 44, .... 20, 29, ....
#> 5 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309    150 f4, 21, .... 02, 20, ....
#> 6 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309    150 11, 21, .... 16, 20, ....

seqqual_decoded <- decode_seqqual_compact(seqqual_compact)
head(seqqual_decoded)
#>                                       qname qwidth
#> 1  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309    150
#> 2  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309    150
#> 3 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309    150
#> 4 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309    150
#> 5 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309    150
#> 6 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309    150
#>                                                                                                                                                      seq
#> 1 NGACCGGCGAGGAATAGGAATCATGGCGGCTGCGCTGTTCGTGCTGCTGGGATTCGCGCTGCTGGGCACCCACGGAGCCTCCGGGGCTGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCGAAGATACTCCTCACCTGCTCCTT
#> 2 TGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCCAAGATACTCCTCACCTGCTCCTTGAATGACAGCGCCACAGAGGTCACAGGGCACCGCTGGCTGAAGGCGGGCGTGGTGCTGAAGGAGGACGCGCTGCCAGGCCAGAAAAN
#> 3 NGGAGCTTCGAGGTGGTATATATGACCGAGCCCATTAATGAGTACTGTGTGCAGCAGCGTAAGCAATTTGATGGGAAGAGCGTGGTCTCAGTTAGCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGGAGGAGCTGAAGAAGATGGAGGAG
#> 4 TGGGAAGAGCCTGGTCTCAGTTACCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGTAGGAGAAGAAGAAGATGGAAGAGAGCAAGGTAAAGTTTGAGAACCTCTGCAAGCTCATGAAAGAAATCTTAGATAAGAAGGTCGAGATGGTGAN
#> 5 NGCATCTCTACCCCTACTGTCCAGTAGGTGGGATGTGGCTGGGCTGGACAGTCCAGATTATACGAACTGGCAACTCTGAACAAACACCCTCCCTGGAAACAATATTATATTTGATGGTTAGATTCTTTAGCAAACCTATTACATTATTCG
#> 6 AACAAAAACCCTCCCTAGCAACAATATTATATTTGATCGTTAGATTCTTTAGCAAACCTATTACCTTATTCGATGTCAGCTAACCTCTTGGTTTGATCATCTTTTCCAGCTGCTCTAGGGGCCTAACCACCTTGGATAACCTGGTCTCTN
#>                                                                                                                                                     qual
#> 1 #AA-AAFJFJJJJJJJJJAAJJJJJJJJJJJJJF<-FJF-AAFJJJJFFJJFFF<FFJ-<<AFAFF7-AFJFJJJJJ77FJAAJJ7AA777AJJAJFJ-A<AA<AAF7FFA-FAFA<-77<FJA--)77A7AA7--7AF<)7777<))--
#> 2 -FAAA)-AJFF7JJFFAF-FJAF-<A-7--JJJFA-F-A7-<-F7-FFJFFJAFJFJJAJFJ<FFJJJFAF7JJFFJFFJ7<-<AA-<AF-JJA<7<777A7<JA77-7JJJJJFJJJFJJJFFJJFJJJJJJJJJF7-JJJF-AFAAA#
#> 3 #AA-AFFJFJJAAFFJAJJJ<AJJJFJJJJJA<<F7-F-<-<FFAFFJFAJJF<<FJ----<<-7<J-FF-AJJJF<J<JF-AFAA7-7AJJJF-7777A-7AA---777---7FJ--AFA7-7AAJFAA)7))7---A<7<--7---A7
#> 4 AJA-FFJJJFF<7FF-AA<7F<AF<-7FJFJJJJJFF-<7-JFJJJJFFAJFJ<A-JF-JJJAFFFFJJJJFJ7<--F-A-<7JF<-FF<<---<-77<-A-FJJFJJFJFJJJJJJJJJJJ<JJJJJJJJJJJJJJ<7JJJF<FFAAA#
#> 5 #AAFFF-FJJJA-AFJJJJJJJJJJJJJJJJJJJJJJJJFJF7-7FJJJFJJJ<<F-FF-FJJFJFJ-F<JA7-7-7AF-JJJFJJJFJJJ7FJJJ<FJJJJJJJJJJJJFJJJJJ7FFJ-7-7-A-AFFF-F--7F7AJJJJJ-<<FJJ
#> 6 7A7-7--<A-)<)AA7-7-A<<7-FF<7<F<7-7A-F-JJAAFJJJJJJA<JF<A7JJJFAA-J<JJFJJA-7<JJJJFFJJAJJA7<<<7<FJJAJJJFJJFJJJA-JJFJJFJJJJJJJJFJJJJJJ<F-JJJJJFFJJFJJFFAAA#
```

### Fast count summaries

[`bam_count()`](https://cparsania.github.io/BamScaleR/reference/bam_count.md)
provides chromosome-level count summaries using the same BAM filtering
model.

``` r

counts <- bam_count(bam, threads = 2)
counts
#>    seqname seqlength count
#> 1        1 248956422  1308
#> 2       10 133797422   308
#> 3       11 135086622   600
#> 4       12 133275309   648
#> 5       13 114364328   162
#> 6       14 107043718   230
#> 7       15 101991189   294
#> 8       16  90338345   334
#> 9       17  83257441   570
#> 10      18  80373285    78
#> 11      19  58617616   600
#> 12       2 242193529   586
#> 13      20  64444167   168
#> 14      21  46709983    76
#> 15      22  50818468   200
#> 16       3 198295559   454
#> 17       4 190214555   252
#> 18       5 181538259   516
#> 19       6 170805979   824
#> 20       7 159345973   412
#> 21       8 145138636   452
#> 22       9 138394717   372
#> 23      MT     16569   288
#> 24       X 156040895   254
#> 25       Y  57227415    14
#> 26       *        NA     0
```

## Summaries computed inside the reader

Many analyses do not need the alignments themselves, only a summary of
them: a fragment-size distribution for ATAC-seq or ChIP-seq quality
control, a mapping-quality profile, or per-base coverage for a
genome-browser track. The usual route reads every alignment into R first
and then summarises it. The four functions in this section instead
compute the summary inside the multithreaded reader. Each thread keeps
its own running totals while it decodes its share of the file, and only
the finished summary is returned to R. The checks below compare each
result with the standard Bioconductor route.

All four functions accept a `ScanBamParam`, of which the `flag` and
`mapqFilter` settings are applied to each read. (They cannot restrict
the summary to regions, so setting `which` gives a warning.) The
following filter (mapped, primary, passing QC) is the usual one for
fragment-size QC:

``` r

qc_flag <- Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE,
    isNotPassingQualityControls = FALSE
)
qc_param <- Rsamtools::ScanBamParam(flag = qc_flag)
```

### Fragment-size distribution: `fragment_sizes()`

[`fragment_sizes()`](https://cparsania.github.io/BamScaleR/reference/fragment_sizes.md)
tabulates the absolute template length (`abs(TLEN)`, the `isize` field
in `Rsamtools`) of every read that passes the filter. It returns a
`data.frame` with one row per observed fragment size, in increasing
order, and the number of reads with that size. Reads that are unmapped,
or whose mate is unmapped, have no defined fragment size and are
excluded by default (`drop_mate_unmapped = TRUE`), just as `scanBam()`
reports `NA` for them.

``` r

fs <- fragment_sizes(bam, param = qc_param, threads = 2)
str(fs)
#> 'data.frame':    2047 obs. of  2 variables:
#>  $ fragment_size: int  59 76 78 79 84 85 86 89 91 93 ...
#>  $ count        : num  2 2 2 2 4 2 2 4 4 2 ...

## the most common fragment size
fs[which.max(fs$count), ]
#>     fragment_size count
#> 122           212    62
```

``` r

short <- fs[fs$fragment_size <= 1000, ]
plot(short$fragment_size, short$count, type = "h",
     xlab = "Fragment size (bp)", ylab = "Reads")
```

![Fragment-size distribution of the example BAM (fragments up to 1
kb).](BamScaleR-intro_files/figure-html/fragment-sizes-plot-1.png)

Fragment-size distribution of the example BAM (fragments up to 1 kb).

The same table computed the standard way, by reading every `isize` into
R and tabulating it, is identical:

``` r

isize <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "isize", flag = qc_flag)
)[[1]]$isize
tab <- table(abs(isize))
ref <- data.frame(fragment_size = as.integer(names(tab)),
                  count = as.numeric(tab))
identical(fs, ref)
#> [1] TRUE
```

### Mapping-quality distribution: `mapq_dist()`

[`mapq_dist()`](https://cparsania.github.io/BamScaleR/reference/mapq_dist.md)
returns the number of reads at each mapping quality (MAPQ, 0 to 255) as
a `data.frame` with columns `mapq` and `count`. Unmapped reads have no
mapping quality (`scanBam()` reports `NA` for them) and are not counted.

``` r

mq <- mapq_dist(bam, param = qc_param, threads = 2)
mq
#>   mapq count
#> 1  255 10000
```

It is identical to tabulating the `mapq` values read with `scanBam()`:

``` r

mapq <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "mapq", flag = qc_flag)
)[[1]]$mapq
tab <- table(mapq)
identical(mq, data.frame(mapq = as.integer(names(tab)),
                         count = as.numeric(tab)))
#> [1] TRUE
```

The example BAM was aligned with STAR, which gives every uniquely mapped
read a MAPQ of 255, so the distribution has a single value. On a BAM
from an aligner such as Bowtie 2 or BWA, the table shows the full spread
of MAPQ values, and setting `mapqFilter` in `param` removes reads below
that value from the distribution.

### Per-base coverage: `bam_coverage()`

[`bam_coverage()`](https://cparsania.github.io/BamScaleR/reference/bam_coverage.md)
returns per-base read depth as an `RleList`: one run-length encoded
integer vector per sequence in the BAM header, in header order, each as
long as the sequence. Sequences with no reads are a single run of zeros.
CIGAR operations are counted in the same way as
[`GenomicAlignments::coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html):
`M`, `=`, `X` and `D` add depth, `N` (a splice junction) does not, and
`I`, `S`, `H` and `P` are ignored.

``` r

cov <- bam_coverage(bam, threads = 2)
cov[c("1", "MT")]
#> RleList of length 2
#> $`1`
#> integer-Rle of length 248956422 with 3132 runs
#>   Lengths:  630211      19       4       4 ...      26     123      26 1857649
#>   Values :       0       1       2       3 ...       1       2       1       0
#> 
#> $MT
#> integer-Rle of length 16569 with 518 runs
#>   Lengths: 140   6  24  34  44  42   4  19 ...  32  41  82 520 127  19 131 684
#>   Values :   0   1   2   4   5   6   5   4 ...   1   2   1   0   1   2   1   0
```

The result is an ordinary `RleList`, so the usual `IRanges` tools apply,
for example per-chromosome summaries or the regions above a depth
threshold:

``` r

head(data.frame(
    seqname = names(cov),
    covered_bp = sum(cov > 0),
    max_depth = max(cov)
))
#>    seqname covered_bp max_depth
#> 1        1     108800        21
#> 10      10      28831        10
#> 11      11      52762         9
#> 12      12      52669        12
#> 13      13      15105        11
#> 14      14      20147        13

## regions of the mitochondrial genome covered by at least 10 reads
IRanges::slice(cov[["MT"]], lower = 10)
#> Views on a 16569-length Rle subject
#> 
#> views:
#>      start  end width
#>  [1]  2502 2597    96 [12 12 12 13 13 13 13 13 13 13 13 13 13 13 12 12 12 ...]
#>  [2]  6909 6988    80 [10 10 10 10 10 10 10 10 11 11 11 11 11 11 11 11 11 ...]
#>  [3]  7046 7048     3 [10 10 10]
#>  [4]  7107 7124    18 [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 ...]
#>  [5]  7158 7277   120 [10 10 10 10 10 10 11 11 11 11 11 11 12 12 12 12 12 ...]
#>  [6]  7301 7313    13 [10 10 10 10 10 10 10 10 10 10 10 10 10]
#>  [7]  7345 7352     8 [10 10 10 10 10 10 10 10]
#>  [8]  7819 7829    11 [10 10 10 10 10 10 10 10 10 10 10]
#>  [9]  8095 8115    21 [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 ...]
#> [10]  8146 8203    58 [10 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 ...]
#> [11]  9766 9820    55 [10 10 10 10 10 10 11 11 11 11 11 11 11 11 10 10 10 ...]
#> [12]  9826 9932   107 [10 10 10 10 10 10 11 11 11 11 13 13 13 13 13 13 13 ...]
```

It is identical to the coverage computed from a `GAlignments` object:

``` r

ref_cov <- GenomicAlignments::coverage(
    GenomicAlignments::readGAlignments(bam)
)
identical(cov, ref_cov)
#> [1] TRUE
```

### Coverage straight to a bigWig file: `bam_coverage_bigwig()`

To make a genome-browser track,
[`bam_coverage_bigwig()`](https://cparsania.github.io/BamScaleR/reference/bam_coverage_bigwig.md)
computes the same coverage and writes it directly to a bigWig file,
without creating the coverage `RleList` in R. It returns the path of the
file it wrote (invisibly). Only non-zero runs are written, which is the
bigWig convention; uncovered bases read back as zero.

``` r

bw <- tempfile(fileext = ".bw")
bam_coverage_bigwig(bam, outfile = bw, threads = 2)
file.size(bw)
#> [1] 271289
```

The file can be loaded into a genome browser or read back into R with
`rtracklayer`. The values match
[`bam_coverage()`](https://cparsania.github.io/BamScaleR/reference/bam_coverage.md)
on every sequence (a bigWig stores them as floating-point numbers, so
they are compared as numbers):

``` r

imported <- rtracklayer::import(bw, as = "RleList")
all(vapply(names(cov), function(chr) {
    isTRUE(all.equal(as.numeric(imported[[chr]]), as.numeric(cov[[chr]])))
}, logical(1)))
#> [1] TRUE
```

### Several files at once

Every function accepts several BAM files, as a character vector or an
[`Rsamtools::BamFileList`](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html),
and returns a named list with one result per file. Files are processed
in parallel according to `BPPARAM`, and `threads` sets the number of
threads used within each file. The `BPPARAM` object is used as supplied
and is not modified, so it can be reused for other work afterwards. (On
Linux, GCC’s OpenMP runtime cannot start threads in a forked
`MulticoreParam` worker once the R session has used OpenMP threads, and
the worker would wait forever. `BamScaleR` releases the runtime’s idle
threads before each `fork()`, so forked workers use `threads` as
requested. This needs libgomp 9 or later; with an older runtime, use
`SnowParam()` instead.)

``` r

two <- fragment_sizes(
    c(rep1 = bam, rep2 = bam),
    param = qc_param,
    threads = 1,
    BPPARAM = BiocParallel::SerialParam()
)
vapply(two, function(x) sum(x$count), numeric(1))
#>  rep1  rep2 
#> 10000 10000
```

## Relationship to existing Bioconductor tools

`BamScaleR` is intended to complement, not replace, existing
Bioconductor packages.

- `Rsamtools` remains the canonical low-level BAM access layer in R and
  defines the filtering idioms that `BamScaleR` reuses through
  `ScanBamParam`.
- `GenomicAlignments` remains the standard package for alignment-centric
  downstream workflows, and `BamScaleR` supports direct generation of
  compatible alignment objects.
- `BiocParallel` remains the standard mechanism for file-level
  parallelism, and `BamScaleR` adds a second axis of parallelism through
  per-file OpenMP threads.

The main difference is therefore performance-oriented: `BamScaleR`
accelerates the traversal step itself while staying close to existing
Bioconductor usage patterns.

### BamScaleR and ompBAM

`ompBAM` and `BamScaleR` work at different levels.

`ompBAM` is a header-only C++ library for package developers. It opens a
BAM file, decompresses its compressed (BGZF) blocks on several OpenMP
threads, and gives each thread a share of the raw alignment records,
with C++ functions to read the fields of one record (flag, position,
CIGAR, tags and so on). It has no R function that reads a BAM file into
R; its R functions are helpers for developers, such as `use_ompBAM()`
and `example_BAM()`. To do anything with the records, a package must
write its own C++ loop over them.

`BamScaleR` is that loop, together with what is needed to use its output
in Bioconductor workflows:

- **R functions with `Rsamtools` semantics.**
  [`bam_read()`](https://cparsania.github.io/BamScaleR/reference/bam_read.md)
  takes the `what` fields of
  [`Rsamtools::scanBam()`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html)
  (except `groupid` and `mate_status`), its `tag` fields and
  `ScanBamParam` settings (`flag`, `mapqFilter`, `which`), and applies
  them to each record in C++. (`which` selects the alignments that
  overlap a range, as in `Rsamtools`;
  [`?bam_read`](https://cparsania.github.io/BamScaleR/reference/bam_read.md)
  lists the remaining differences.)
  [`bam_count()`](https://cparsania.github.io/BamScaleR/reference/bam_count.md)
  and the summary functions apply the filters listed on their help
  pages, and warn about settings they cannot apply.
- **Bioconductor output objects.** Results from all threads are merged
  in file order and converted into a `data.frame` or `DataFrame`, a
  `scanBam()`-style list (with `DNAStringSet` and `PhredQuality` for
  sequences and qualities), a `GAlignments` or a `GAlignmentPairs`
  object; the columns and the `GAlignments` slots are built in compiled
  code. The package tests compare fields, sequences and qualities with
  [`Rsamtools::scanBam()`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html),
  and `GAlignments` objects with `readGAlignments()`, using
  [`identical()`](https://rdrr.io/r/base/identical.html). They also
  check that, without `which`, `GAlignmentPairs` output contains the
  same pairs as `readGAlignmentPairs()`, including for secondary and
  supplementary alignments.
- **Summaries inside the reader.** Per-chromosome counts, fragment-size
  and mapping-quality distributions, per-base coverage, and a bigWig
  writer (a bundled copy of libBigWig, extended to compress blocks in
  parallel).
- **Two levels of parallelism.** OpenMP threads within a file and
  `BiocParallel` workers across files. With `auto_threads = TRUE` the
  two are balanced against the number of available cores.
- **Portable builds.** The `configure` script described under
  [Installation](#installation) finds a working OpenMP setup on macOS.

## Session information

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.32.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] BamScaleR_0.99.5 BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] SummarizedExperiment_1.42.0 rjson_0.2.23               
#>  [3] xfun_0.62                   bslib_0.12.0               
#>  [5] Biobase_2.72.0              lattice_0.23-1             
#>  [7] tools_4.6.2                 bitops_1.1-0               
#>  [9] generics_0.1.4              stats4_4.6.2               
#> [11] curl_8.1.0                  parallel_4.6.2             
#> [13] Matrix_1.7-6                desc_1.4.3                 
#> [15] S4Vectors_0.50.3            cigarillo_1.2.1            
#> [17] lifecycle_1.0.5             compiler_4.6.2             
#> [19] Rsamtools_2.28.0            textshaping_1.1.0          
#> [21] Biostrings_2.80.2           Seqinfo_1.2.0              
#> [23] codetools_0.2-20            GenomeInfoDb_1.48.0        
#> [25] htmltools_0.5.9             sass_0.5.0                 
#> [27] RCurl_1.98-1.20             yaml_2.3.12                
#> [29] pkgdown_2.2.1               crayon_1.5.3               
#> [31] jquerylib_0.1.4             BiocParallel_1.46.0        
#> [33] cachem_1.1.0                DelayedArray_0.38.2        
#> [35] abind_1.4-8                 ompBAM_1.16.0              
#> [37] digest_0.6.39               restfulr_0.0.17            
#> [39] bookdown_0.48               fastmap_1.2.0              
#> [41] grid_4.6.2                  cli_3.6.6                  
#> [43] SparseArray_1.12.3          S4Arrays_1.12.1            
#> [45] XML_3.99-0.25               UCSC.utils_1.8.0           
#> [47] rmarkdown_2.32              XVector_0.52.0             
#> [49] httr_1.4.9                  matrixStats_1.5.0          
#> [51] otel_0.3.0                  ragg_1.5.2                 
#> [53] evaluate_1.0.5              knitr_1.52                 
#> [55] GenomicRanges_1.64.0        IRanges_2.46.0             
#> [57] BiocIO_1.22.0               rtracklayer_1.72.0         
#> [59] rlang_1.3.0                 Rcpp_1.1.2                 
#> [61] BiocManager_1.30.27         BiocGenerics_0.58.1        
#> [63] jsonlite_2.0.0              R6_2.6.1                   
#> [65] MatrixGenerics_1.24.0       GenomicAlignments_1.48.0   
#> [67] systemfonts_1.3.4           fs_2.1.0
```
