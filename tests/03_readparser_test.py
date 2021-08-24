#
# Copyright (c) 2021
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of pycroquet.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
import io
import os

import pysam
import pytest

from pycroquet import readparser
from pycroquet.classes import Stats

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
)


@pytest.mark.parametrize(
    "file, result, info",
    [
        ("good_guide_01.tsv", False, "Not gzip"),
        ("fastq/unique.fq.gz", True, "Gzip"),
    ],
)
def test_01_readparser_gzip(file, result, info):
    assert readparser.is_gzip(os.path.join(DATA_DIR, file)) == result, info


@pytest.mark.parametrize(
    "header, result, info",
    [
        (
            "@HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436 1:N:0:GACGACGT",
            ("HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436", 1, False),
            "Casava 1.8",
        ),
        (
            "@A00471:89:HMTWVDMXX:1:1101:2871:1016 1:N:0:TTGGACGT+AGCACTTC",
            ("A00471:89:HMTWVDMXX:1:1101:2871:1016", 1, False),
            "Cassava dual-bc",
        ),
        (
            "@HS27_17643:2:2110:8108:93084#6/1",
            ("HS27_17643:2:2110:8108:93084#6", 1, False),
            "Illumina read 1",
        ),
        (
            "@HS27_17643:2:2110:8108:93084#6/2",
            ("HS27_17643:2:2110:8108:93084#6", 2, False),
            "Illumina read 2",
        ),
        (
            "@HS27_17643:2:2110:8108:93084#6",
            ("HS27_17643:2:2110:8108:93084#6", None, False),
            "Illumina single",
        ),
        (
            "@HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436 1:Y:0:GACGACGT",
            ("HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436", 1, False),
            "Casava 1.8 - failed filtering - keep",
        ),
        (
            "@HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436 1:Y:0:GACGACGT",
            ("HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436", 1, True),
            "Casava 1.8 - failed filtering - discard",
        ),
    ],
)
def test_02_readparser_parse_fq_header(header, result, info):
    readparser.parse_fq_header(header) == result, info


def test_03_readparser_parse_fq_header():
    with pytest.raises(ValueError, match="^Unsupported FastQ header format: ") as e_info:
        readparser.parse_fq_header("@HS27_17643:2:2110:8108:93084#6/9")


@pytest.mark.parametrize(
    "file, result, info",
    [
        (
            "fastq/unique.fq",
            (
                3,
                Stats(total_reads=3, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 1,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "unique",
        ),
        (
            "fastq/non-unique.fq",
            (
                3,
                Stats(total_reads=4, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 2,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "non-unique",
        ),
    ],
)
def test_04_readparser_parse_fastq(file, result, info):
    assert readparser.parse_fastq(open(os.path.join(DATA_DIR, file)), "bob") == result, info


@pytest.mark.parametrize(
    "file, result, info",
    [
        ("htsfile/unique.sam", "PD13371a", "single RG"),
        ("htsfile/multi_rg.sam", "PD13371a", "multi RG"),
        ("htsfile/rg_no_sm.sam", None, "single RG, no SM"),
    ],
)
def test_05_readparser_hts_sample(file, result, info):
    save = pysam.set_verbosity(0)
    sam = pysam.AlignmentFile(os.path.join(DATA_DIR, file), mode="r", require_index=False)
    pysam.set_verbosity(save)
    assert readparser._hts_sample(sam) == result, info


def test_06_readparser_hts_multi_sample():
    save = pysam.set_verbosity(0)
    sam = pysam.AlignmentFile(
        os.path.join(DATA_DIR, "htsfile", "multi_sample.sam"),
        mode="r",
        require_index=False,
    )
    pysam.set_verbosity(save)
    with pytest.raises(ValueError, match="^Multiple different sample names found in header") as e_info:
        readparser._hts_sample(sam)


@pytest.mark.parametrize(
    "file, exclude_qcf, result, info",
    [
        (
            "fastq/unique.fq",
            False,
            (
                3,
                Stats(total_reads=3, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 1,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "fa-unique",
        ),
        (
            "fastq/unique.fq.gz",
            False,
            (
                3,
                Stats(total_reads=3, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 1,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "fa.gz-unique",
        ),
        (
            "htsfile/unique.sam",
            False,
            (
                3,
                Stats(total_reads=4, vendor_failed_reads=1, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 1,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 2,
                },
                None,
            ),
            "sam-unique",
        ),
        (
            "htsfile/unique.sam",
            True,
            (
                3,
                Stats(total_reads=3, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 1,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "sam-unique-exclude-qc-fail",
        ),
        (
            "htsfile/unique.bam",
            False,
            (
                3,
                Stats(total_reads=3, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 1,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "bam-unique",
        ),
        (
            "htsfile/unique.cram",
            False,
            (
                3,
                Stats(total_reads=3, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 1,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "cram-unique",
        ),
        (
            "htsfile/non-unique.sam",
            False,
            (
                3,
                Stats(total_reads=5, vendor_failed_reads=1, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 3,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "sam-non-unique",
        ),
        (
            "htsfile/non-unique.bam",
            False,
            (
                3,
                Stats(total_reads=4, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 2,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "bam-non-unique",
        ),
        (
            "htsfile/non-unique.cram",
            False,
            (
                3,
                Stats(total_reads=4, sample_name="bob"),
                {
                    "TTCACCGAGCTTCCGGGAG": 2,
                    "TTCANCGAGCTTCCGGGAG": 1,
                    "TTCATCGAGCTTCCGGGAG": 1,
                },
                None,
            ),
            "cram-non-unique",
        ),
    ],
)
def test_10_readparser_parse_reads(file, exclude_qcf, result, info):
    assert (
        readparser.parse_reads(os.path.join(DATA_DIR, file), "bob", exclude_qcfail=exclude_qcf, cpus=1) == result
    ), info


@pytest.mark.parametrize(
    "str_file, exclude_qcf, result, info",
    [
        (
            "@HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436 1:N:0:GACGACGT\nACGT\n+\nMMMM\n",
            False,
            (1, Stats(total_reads=1, sample_name="bob"), {"ACGT": 1}, None),
            "Good read, (qc keep)",
        ),
        (
            "@HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436 1:N:0:GACGACGT\nACGT\n+\nMMMM\n",
            True,
            (1, Stats(total_reads=1, sample_name="bob"), {"ACGT": 1}, None),
            "Good read, (qc discard)",
        ),
        (
            "@HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436 1:Y:0:GACGACGT\nACGT\n+\nMMMM\n",
            False,
            (
                1,
                Stats(total_reads=1, vendor_failed_reads=1, sample_name="bob"),
                {"ACGT": 1},
                None,
            ),
            "Bad read, keep",
        ),
        (
            "@HISEQ2500-01:110:H7AGVADXX:1:1101:10737:10436 1:Y:0:GACGACGT\nACGT\n+\nMMMM\n",
            True,
            (0, Stats(total_reads=0, sample_name="bob"), {}, None),
            "Bad read, discard",
        ),
    ],
)
def test_11_readparser_qcreads(str_file, exclude_qcf, result, info):
    assert readparser.parse_fastq(io.StringIO(str_file), "bob", exclude_qcfail=exclude_qcf) == result, info
