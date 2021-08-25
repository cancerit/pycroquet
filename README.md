# pyCROQUET

python Crispr Read to Oligo QUantification Enhancement Tool

[![cancerit](https://circleci.com/gh/cancerit/pycroquet.svg?style=svg)](https://circleci.com/gh/cancerit/pycroquet)

- [Publications](#publications)
- [General](#general)
- [Subcommands](#subcommands)
- [Options](#options)
  - [`guidelib`](#guidelib)
  - [`queries`](#queries)
  - [`chunks`](#chunks)
  - [`rules`](#rules)
- [Output files](#output-files)
  - [CRAM](#cram)
    - [Reads mapping to a `sgrna_id`](#reads-mapping-to-a-sgrna_id)
    - [Reads assigned to a guide](#reads-assigned-to-a-guide)
- [Dual guide](#dual-guide)
  - [Statistic file extension](#statistic-file-extension)
- [Boundary mode details](#boundary-mode-details)
  - [`exact`](#exact)
  - [`TinQ` - target in query](#tinq---target-in-query)
  - [`QinT` - query in target](#qint---query-in-target)
  - [`any`](#any)
- [Viewing alignments](#viewing-alignments)
- [Installation](#installation)
  - [Pypi](#pypi)
  - [Docker and Singularity](#docker-and-singularity)
- [Development](#development)
  - [Linux](#linux)
  - [Mac](#mac)
- [Testing](#testing)
  - [Local `venv` testing](#local-venv-testing)
  - [Local `pre-commit` hooks](#local-pre-commit-hooks)
  - [Docker testing](#docker-testing)
  - [CI tests](#ci-tests)
  - [Updating licence headers](#updating-licence-headers)
- [LICENSE](#license)

## Publications

Please contact the following for appropriate referencing methods:

- Keiran Raine (kr2@sanger.ac.uk)
- Emre Karakoc (ek11@sanger.ac.uk)
- Victoria Offord (vo1@sanger.ac.uk)

## General

Code in place to support read input from any of the following formats:

- fastq (also gzip compressed)
- sam
- bam
- cram

## Subcommands

- single-guide
  - Short single end read quantification.
- dual-guide
  - Paired end read quantification.
- long-read
  - Long single end read quantification
- guides-to-fa
  - Convert guides to fasta for use with `samtools tview`

## Options

All options are not applicable to all subcommands, however the majority are common.

### `guidelib`

Please see the [Guide library format][guide-format] for a description of this file.

### `queries`

Currently the `dual-guide` mode only supports SAM/BAM/CRAM as input.  Convert fastq to unmapped CRAM with:

```
# if data has casava read barcode/qc (text after space in read name) please add "-i"
samtools import  --output-fmt CRAM,no_ref=1 -@ 4 -1 $READ1 -2 $READ2 -o $OUTFILE.cram
```

### `chunks`

Chunks should be set to a value that allows all CPUs to be utilized.  The value is multiplied by the number of
CPUs requested and this give the number of unique read sequences held in memory during the mapping phase.

This has a direct impact on memory. The value is automatically reduced when too large to allow full use of requested CPUs.

### `rules`

For single-guide `--rules MM` (allow 2 mismatches in alignment) is a sensible value.  For other subcommands the decision
is dependent on the library protocol.

Rules have a direct impact on run time as they increase the time taken to abort an alignment, individual costs are as follows:

- `M` = 1
- `I` = 2 (single b.p.)
- `D` = 2 (single b.p.)

Performance is only impacted by the maximum penalty you allow.

Be aware if you with to allow up to `2 mismatch` or `1 mismatch + 1 b.p. insert` you must specify:

```
pycroquet ... --rules MM --rules MI
```

## Output files

### CRAM

For single-guide you have the option to suppress it via the `--no-alignment` (`-n`) option.  In dual-guide it is tightly
linked to the pairing code so not possible to disable.

Reads that map uniquely are written with `MAPQ>0` (score calculations have not been refined at this time).  There are some
differences in how to interpret the data depending on if you are processing single-guide, dual-guide or long-read.

#### Reads mapping to a `sgrna_id`

To get the reads that map uniquely to a guide element (sgrna_seq) use the `sgrna_id`.  This is primarily of use for `single-guide` and `long-read`:

```bash
samtools view -F 4 -q 1 result.cram $SGRNA_ID
```

To get a single instance of reads that map to a guide element but map equally well to others select for the `SA` tag (requires samtools>=1.12):

```bash
samtools view -F 4 -F 256 -d SA result.cram $SGRNA_ID
```

To get reads that failed to map:

```bash
samtools view -f 4 result.cram
```

#### Reads assigned to a guide

This is only applicable to `dual-guide`.

You can pull reads by the guide id using `samtools view`, this example counts R1 mapping to a guide (equivalent to the `*.counts.tsv`
result), replace/set `$GUIDE_ID` as required:

```bash
samtools view -F 4 -f 64 -c -d YG:$GUIDE_ID result.cram
```

To select all the reads mapped to this guide grouped by readname:

```bash
samtools view -u -F 4 -d YG:$GUIDE_ID result.cram | samtools sort -n - | samtools view -b - > $GUIDE_ID.bam
```

## Dual guide

FASTQ(.gz) input is not currently supported for dual guide, please prepare your data appropriately with `samtools import`:

```bash
samtools import -@ 4 -1 R1.fastq.gz -2 R2.fastq.gz -O BAM -o OUTPUT.bam
```

Please review the `import` options as casava information can be interpreted where appropriate.

### Statistic file extension

The dual guide output extends the standard json statistics file adding `pair_classifications`:

| Classification | Description                            |
|----------------|----------------------------------------|
| match          | same vector F/R                        |
| aberrant_match | same vector pair, aberrant orientation |
| f_multi_3p     | 5p mapped F, 3p multihit               |
| f_multi_5p     | 3p mapped F, 5p multihit               |
| r_multi_3p     | 5p mapped R, 3p multihit               |
| r_multi_5p     | 3p mapped R, 5p multihit               |
| f_open_3p      | 5p mapped F, 3p open (unmapped)        |
| f_open_5p      | 3p mapped F, 5p open (unmapped)        |
| r_open_3p      | 5p mapped R, 3p open (unmapped)        |
| r_open_5p      | 3p mapped R, 5p open (unmapped)        |
| swap           | multi vector, uniq mapped              |
| ambiguous      | both ends multi hit                    |
| no_match       | multi/unmapped either end              |

## Boundary mode details

The `-b/--boundary-mode` option controls how the guide and read are allowed to overlap.  Each section shows the types of
alignment allowed, to be valid they still need to pass rules and any minimum score.

In all cases `XXX` indicates original sequence.

#### `exact`

Boundary of sequence must be equal between target (guide) and query (read)

```
T: XXXXXX
Q: XXXXXX
```

#### `TinQ` - target in query

Like the name suggests, valid alignments include those via `exact` and:

```
T: XXXXX
Q: XXXXXX

T:  XXXX
Q: XXXXXX

T:  XXXXX
Q: XXXXXX
```

#### `QinT` - query in target

Reverse of `TinQ`, valid alignments include those via `exact` and:

```
T: XXXXXX
Q: XXXXX

T: XXXXXX
Q:  XXXX

T: XXXXXX
Q:  XXXXX
```

#### `any`

No boundary checks are performed, this allows more complex events, all alignments from `exact`, `TinQ`, `QinT` plus:

```
T:   XXXXXXX
Q: XXXXXXX

T: XXXXXXX
Q:   XXXXXXX
```

## Viewing alignments

You can use samtools tview to view the cram file, this is mainly useful when checking fuzzy matching or allowing all boundary types.

To make this more informative generate the fasta file for the sgrna elements:

```
pycroquet guides-to-fa --guidelib guide_library.tsv --fasta sgrna.fa
```

NOTE: a contig is the individual sgrna sequence, not the pair

Now use this with `samtool tview` to view your alignments, see command line help to jump directly to a contig of interest
or `?` when using interactively.

```
samtools tview result.cram sgrna.fa
```

Will give a full screen output like this (N is just padded to screen width for short contigs):

```
1         11
CTAGTTCAGATAAAACAACNNNNNN
...................
...................
...................
...................
...................
................C..
...................
...................
```

## Installation

### Pypi

```
pip install Cython
pip install pycroquet
```

### Docker and Singularity

There are pre-built images containing this codebase on [quay.io][quay-repo].  When pulling an image you must specify
the version there is no `latest`.

The docker images are known to work correctly after import into a singularity image.

## Development

Python 3.9 or better required.

### Linux

```bash
git clone git@github.com:cancerit/pycroquet.git
cd pycroquet
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python3 ./setup.py develop  # dynamic build

# see later, assumes global install of pre-commit
pre-commit install
```

### Mac

```bash
brew update
brew install python@3.9
brew install libmagic
git clone git@github.com:cancerit/pycroquet.git
cd pycroquet
python3.9 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python3 setup.sh develop
```

## Testing

There are 4 layers to testing and standards:

1. Local `venv` testing
1. Local `pre-commit` hooks
1. Tests embedded in `docker build`
1. `CI` tests

### Local `venv` testing

```bash
/tests/scripts/run_unit_tests.sh
```

### Local `pre-commit` hooks

This project additionally uses git pre-commit hooks via the [pre-commit tool](https://pre-commit.com/).  These are concerned
with file formats and standards, not the actual execution of code.  See `./.pre-commit-config.yaml`.

### Docker testing

The Docker build includes the unit tests, but removes many of the libraries before the final build stage.  Mainly for CI tests.

### CI tests

CI includes 2 additional tests, each based on the 2 datasets in the `./examples` directory.

### Updating licence headers

Please use [skywalking-eyes](https://github.com/apache/skywalking-eyes).

Expected workflow:

1. Check state before modifying `.licenserc.yaml`:
   - `docker run -it --rm -v $(pwd):/github/workspace apache/skywalking-eyes header check`
   - You should get some 'valid' here, those without a header as 'invalid'
1. Modify `.licenserc.yaml`
1. Apply the changes:
   - `docker run -it --rm -v $(pwd):/github/workspace apache/skywalking-eyes header fix`
1. Add/commit changes

This is executed in the CI pipeline.

*DO NOT* edit the header in the files, please modify the date component of `content` in `.licenserc.yaml`.  The only exception being:

- `README.md`

If you need to make more extensive changes to the license carefully test the pattern is functional.

## LICENSE

```
Copyright (c) 2021

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

This file is part of pygas.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’.
```

<!-- refs -->

[guide-format]: https://github.com/cancerit/pycroquet/wiki/Guide-library-format
[quay-repo]: https://quay.io/repository/wtsicgp/pycroquet?tab=tags
