# BamScaleR: Bioconductor-Friendly Multithreaded BAM Processing

BamScaleR reads BAM files using multiple threads per file and returns
the objects that Bioconductor workflows already use, so it can take the
place of Rsamtools::scanBam() and GenomicAlignments::readGAlignments()
in existing analyses. It builds on the ompBAM C++ library, which
decompresses BAM blocks in parallel, and adds readers that select fields
like scanBam(), apply the flag, mapping-quality and region filters of
Rsamtools::ScanBamParam, and return a data.frame, DataFrame,
scanBam-style list, GAlignments or GAlignmentPairs object. It also
computes common summaries inside the multithreaded reader without
creating per-read R objects: per-chromosome read counts, paired-end
fragment-size and mapping-quality distributions, per-base coverage as an
RleList, and coverage written directly to a bigWig file. Several BAM
files can additionally be processed in parallel with BiocParallel.

## See also

Useful links:

- <https://cparsania.github.io/BamScaleR/>

- Report bugs at <https://github.com/cparsania/BamScaleR/issues>

## Author

**Maintainer**: Chirag Parsania <chirag.parsania@gmail.com>

Authors:

- Chirag Parsania <chirag.parsania@gmail.com>
