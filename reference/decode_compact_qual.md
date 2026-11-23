# Decode compact BamScaleR quality output

Decodes `qual` values returned by
`bam_read(..., seqqual_mode = "compact")` back to ASCII Phred-quality
strings.

## Usage

``` r
decode_compact_qual(qual)
```

## Arguments

- qual:

  A list (or list-column) of `raw` vectors produced by compact BamScaleR
  quality extraction.

## Value

A character vector containing decoded quality strings. Entries with
all-missing quality bytes (value 255) are returned as `"*"`, the SAM
notation for a missing quality. (The compatible mode of
[`bam_read()`](https://cparsania.github.io/BamScaleR/reference/bam_read.md)
instead returns
[`Rsamtools::scanBam()`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html)'s
`PhredQuality` for such reads.)

## See also

[`decode_compact_seq()`](https://cparsania.github.io/BamScaleR/reference/decode_compact_seq.md),
[`decode_seqqual_compact()`](https://cparsania.github.io/BamScaleR/reference/decode_seqqual_compact.md),
[`bam_read()`](https://cparsania.github.io/BamScaleR/reference/bam_read.md)

## Examples

``` r
decode_compact_qual(
  qual = list(as.raw(c(0L, 1L, 2L, 3L)))
)
#> [1] "!\"#$"
```
