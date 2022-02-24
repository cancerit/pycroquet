# CHANGES

## 1.5.1

- Improve outputs when `R2_R1` has been applied to ensure input fields are replicated "as is" in counts file
- Add warning to indicate `sgrna_strands` is minimally implemented within the manifest and not applied in the processing (see #13)

## 1.5.0

- Dual-guide count output is now same as single-guide
- All count output files are compressed with gzip
- Merging of count and statistics files implemented

## 1.4.1

Use classifications in the `*.query_class.tsv.gz` file.

## 1.4.0

- Adds `mean_count_per_guide` statistic to all modes.
- Correct bug in `--boundary` mode mapping (TinQ and QinT reversed).
- dual-guide: expand classifications to handle multiple guide-pair hits when one end unique.
- dual-guide: Adds the `*.query_class.tsv.gz` file for use in debugging (see #4).

## 1.3.0

- dual-guide: added library header item to deal with reversed read order when comparing against sgRNA.
- dual-guide: handled data multiplication issue in CRAM outputs, quicker.

## 1.2.1

Handle change to how bgzip data is reported by magic decode.

## 1.2.0

Add `--unique` option to `long-read` mode.  Exits as soon as unique read sequence counts are generated.
Data is added to outputs when run without this flag.

## 1.1.1

Correct license text, referencing child package

## 1.1.0

- First public release
- Includes dual-guide and long-read (single-guide, suitable for SGE)

## 1.0.0

- Exact flag dropped in preference for deducing from lack of rules
- Allow substring matching fall through to cope with length mismatch
- Add `--no-alignment` option

## 0.1.1

- Correction in dockerfile
- Internally handle rescaling phred scores when input is fastq

## 0.1.0

- Uses pySSA 0.3.0
- Initial release for user interaction
