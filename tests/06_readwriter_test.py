#
# Copyright (c) 2021-2022
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
import os
import tempfile
from array import array

import pytest

from pycroquet import readwriter
from pycroquet.classes import Seqread
from pycroquet.classes import Stats

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


@pytest.mark.parametrize(
    "folder, file, qcfail, exp_count, info",
    [
        ("htsfile", "unique.sam", False, 4, "SAM file, unique"),
        ("htsfile", "unique.bam", False, 3, "BAM file, unique"),
        ("htsfile", "unique.cram", False, 3, "CRAM file, unique"),
        ("htsfile", "non-unique.sam", False, 5, "SAM file, non-unique"),
        ("htsfile", "non-unique.sam", True, 4, "SAM file, non-unique, exclude qcf"),
    ],
)
def test_01_readwriter_iter_hts(folder, file, qcfail, exp_count, info):
    reads = 0
    for i in readwriter._hts_iter(
        os.path.join(DATA_DIR, folder, file),
        "rb",
        exclude_qcfail=qcfail,
        reverse=False,
        reference=None,
    ):
        reads += 1
    assert reads == exp_count, info


@pytest.mark.parametrize(
    "folder, file, qual_offset, qcfail, exp_count, info",
    [
        ("fastq", "unique.fq", 33, False, 3, "FASTQ file, unique"),
        ("fastq", "non-unique.fq", 33, False, 4, "FASTQ file, non-unique"),
        (
            "fastq",
            "unique_qcfail_casava.fq",
            None,
            True,
            3,
            "FASTQ file, qcfail casava",
        ),
    ],
)
def test_02_readwriter_iter_fastq(folder, file, qual_offset, qcfail, exp_count, info):
    ifh = open(os.path.join(DATA_DIR, folder, file), "rt")
    reads = 0
    for i in readwriter._fq_iter(ifh, exclude_qcfail=qcfail, reverse=False, offset_override=qual_offset):
        reads += 1
    assert reads == exp_count, info


@pytest.mark.parametrize(
    "folder, file, qual_offset, reverse, sequence, quals, info",
    [
        (
            "fastq",
            "unique.fq",
            33,
            False,
            "TTCACCGAGCTTCCGGGAG",
            array(
                "i",
                [
                    33,
                    33,
                    33,
                    33,
                    33,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                ],
            ),
            "FASTQ file, unique",
        ),
        (
            "fastq",
            "unique.fq.gz",
            33,
            False,
            "TTCACCGAGCTTCCGGGAG",
            array(
                "i",
                [
                    33,
                    33,
                    33,
                    33,
                    33,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                ],
            ),
            "FASTQ file, unique.gz",
        ),
        (
            "fastq",
            "unique.fq",
            33,
            True,
            "CTCCCGGAAGCTCGGTGAA",
            array(
                "i",
                [
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    33,
                    33,
                    33,
                    33,
                    33,
                ],
            ),
            "FASTQ file, unique reversed",
        ),
        (
            "htsfile",
            "unique.sam",
            33,
            False,
            "TTCACCGAGCTTCCGGGAG",
            array(
                "i",
                [
                    33,
                    33,
                    33,
                    33,
                    33,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                ],
            ),
            "SAM file, unique",
        ),
        (
            "htsfile",
            "unique.sam",
            33,
            True,
            "CTCCCGGAAGCTCGGTGAA",
            array(
                "i",
                [
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    33,
                    33,
                    33,
                    33,
                    33,
                ],
            ),
            "SAM file, unique reversed",
        ),
        (
            "htsfile",
            "unique.bam",
            33,
            False,
            "TTCACCGAGCTTCCGGGAG",
            array(
                "i",
                [
                    33,
                    33,
                    33,
                    33,
                    33,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                ],
            ),
            "BAM file, unique",
        ),
        (
            "htsfile",
            "unique.bam",
            33,
            True,
            "CTCCCGGAAGCTCGGTGAA",
            array(
                "i",
                [
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    33,
                    33,
                    33,
                    33,
                    33,
                ],
            ),
            "BAM file, unique reversed",
        ),
        (
            "htsfile",
            "unique.cram",
            33,
            False,
            "TTCACCGAGCTTCCGGGAG",
            array(
                "i",
                [
                    33,
                    33,
                    33,
                    33,
                    33,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                ],
            ),
            "CRAM file, unique",
        ),
        (
            "htsfile",
            "unique.cram",
            33,
            True,
            "CTCCCGGAAGCTCGGTGAA",
            array(
                "i",
                [
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    37,
                    33,
                    33,
                    33,
                    33,
                    33,
                ],
            ),
            "CRAM file, unique reversed",
        ),
    ],
)
def test_03_readwriter_read_iter(folder, file, qual_offset, reverse, sequence, quals, info):
    iter = readwriter.read_iter(os.path.join(DATA_DIR, folder, file), reverse, False, None, qual_offset)
    seqread = next(iter)
    assert seqread.sequence == sequence, f"{info} - seq"
    assert seqread.qual == quals, f"{info} - qual"
    del iter


@pytest.mark.parametrize(
    "seqread, qname, seq, reverse, qual_len, is_qcfail, is_reverse, is_paired, is_read1, is_read2, tags",
    [
        (
            Seqread(
                qname="qname",
                sequence="AAAA",
                qual=array("i", [33, 33, 33, 33]),
                rgid=1,
                qc_fail=False,
                member=None,
            ),
            "qname",
            "AAAA",
            False,
            4,
            False,
            False,
            False,
            False,
            False,
            [],
        ),
        (
            Seqread(
                qname="qname",
                sequence="AAAA",
                qual=array("i", [33, 33, 33, 33]),
                rgid=1,
                qc_fail=True,
                member=None,
            ),
            "qname",
            "AAAA",
            False,
            4,
            True,
            False,
            False,
            False,
            False,
            [],
        ),
        (
            Seqread(
                qname="qname",
                sequence="AAAA",
                qual=array("i", [33, 33, 33, 33]),
                rgid=1,
                qc_fail=False,
                member=1,
            ),
            "qname",
            "AAAA",
            False,
            4,
            False,
            False,
            True,
            True,
            False,
            [],
        ),
        (
            Seqread(
                qname="qname",
                sequence="AAAA",
                qual=array("i", [33, 33, 33, 33]),
                rgid=1,
                qc_fail=False,
                member=2,
            ),
            "qname",
            "AAAA",
            False,
            4,
            False,
            False,
            True,
            False,
            True,
            [("SA", "stuff")],
        ),
    ],
)
def test_04_readwriter_to_alignment(
    seqread,
    qname,
    seq,
    reverse,
    qual_len,
    is_qcfail,
    is_reverse,
    is_paired,
    is_read1,
    is_read2,
    tags,
):
    alignseg = readwriter.to_alignment(seqread, reverse, tags)
    assert alignseg.query_name == qname
    assert alignseg.is_qcfail == is_qcfail
    assert alignseg.query_sequence == seq
    assert len(alignseg.query_qualities) == qual_len
    assert alignseg.is_paired == is_paired
    assert alignseg.is_reverse == is_reverse
    assert alignseg.is_paired == is_paired
    assert alignseg.is_read1 == is_read1
    assert alignseg.is_read2 == is_read2


def test_05_readwriter_to_alignment_except():
    with pytest.raises(ValueError) as e_info:
        readwriter.to_alignment(
            Seqread(qname="q", sequence="A", qual=array("i", [33]), rgid="1", member=3),
            False,
            [],
        )


@pytest.mark.parametrize(
    "file, stats, exp_rgs, exp_pgs",
    [
        (
            os.path.join(DATA_DIR, "htsfile/unique.sam"),
            Stats(sample_name="bob", command="pycroquet single-guide ...."),
            1,
            1,
        ),
        (
            os.path.join(DATA_DIR, "htsfile/no_rg_1_pg.sam"),
            Stats(sample_name="bob", command="pycroquet single-guide ...."),
            1,
            2,
        ),
        (
            os.path.join(DATA_DIR, "fastq/unique.fq"),
            Stats(sample_name="bob", command="pycroquet single-guide ...."),
            1,
            1,
        ),
    ],
)
def test_06_readwriter_rg_pg(file, stats, exp_rgs, exp_pgs):
    (rgs, pgs) = readwriter.rg_pg(file, stats)
    assert len(rgs) == exp_rgs
    assert len(pgs) == exp_pgs
