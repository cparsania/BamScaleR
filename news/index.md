# Changelog

## BamScaleR 0.99.0

### Package rename

- The package has been renamed from **BamScale** to **BamScaleR**. The
  former name collided with the existing NCBI command-line tool BAMscale
  (<https://github.com/ncbi/BAMscale>; Pongor et al., *Epigenetics &
  Chromatin* 2020, <doi:10.1186/s13072-020-00343-x>), which performs a
  different task: normalised (“scaled”) coverage-track generation and
  peak quantification. The rename removes that ambiguity.
- No user-facing functionality changed. All exported functions keep
  their names, arguments and behaviour.
- The option `BamScale.ga_fastpath` is now `BamScaleR.ga_fastpath`.
- BamScale was accepted into Bioconductor devel on 4 May 2026 and was
  never part of a numbered Bioconductor release, so no deprecation cycle
  applies. BamScaleR is submitted as its renamed continuation and
  version numbering restarts at 0.99.0, as required for a new
  submission. The entries below document this package’s history under
  its former name.
