# CHANGES

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
