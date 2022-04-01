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
import json
import os
import tempfile

import pytest

import pycroquet.merge as merge

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


def test_01_merge_header_line():
    h_line = merge.merge_header_line("command", "version", 1, "chksum")
    assert h_line == "##Count-col-#1: chksum; version; command"


def test_02_hash_file():
    assert merge.hash_file(f"{DATA_DIR}/merge/bob_1.counts.tsv", "md5") == "md5: 41aacd3ea72f328b1d8ce6427ae6f28a"
    assert (
        merge.hash_file(f"{DATA_DIR}/merge/bob_1.counts.tsv", "sha256")
        == "sha256: 05197cd393645e123b450b1fbda4fbaa4c49245f2d1ae3faf7f82269987b943a"
    )


@pytest.mark.parametrize(
    "inputs, low_count, compare_to",
    [
        (
            ["cli/output/bob_1.counts.tsv.gz", "cli/output/bob_2.counts.tsv.gz"],
            None,
            "merge/bob_1n2_nolow",
        ),
        (
            ["cli/output/bob_1.counts.tsv.gz", "cli/output/bob_2.counts.tsv.gz"],
            5,
            "merge/bob_1n2_lowcount",
        ),
        (
            ["cli/output/dual_low_count.counts.tsv.gz", "cli/output/dual_low_count_2.counts.tsv.gz"],
            None,
            "merge/dual_low",
        ),
    ],
)
def test_03_merge(inputs, low_count, compare_to):
    checksum = "md5"
    loglevel = "WARN"

    stats_new = None

    inputs = [os.path.join(DATA_DIR, input) for input in inputs]

    with tempfile.TemporaryDirectory() as tdir:
        output = os.path.join(tdir, "single")
        merge.merge_counts(output, inputs, low_count, checksum, loglevel)

        out_counts = f"{output}.merged-counts.tsv.gz"
        out_stats = f"{output}.merged-stats.json"

        assert os.path.exists(out_counts)
        assert os.path.exists(out_stats)

        with open(out_stats, "r") as sfp:
            stats_new = json.load(sfp)

    stats_old = None
    with open(os.path.join(DATA_DIR, f"{compare_to}.merged-stats.json"), "r") as sfp:
        stats_old = json.load(sfp)

    for i in ("command", "version"):
        del stats_old[i]
        del stats_new[i]
    assert stats_new == stats_old


def test_04_samefile(capsys):
    checksum = "md5"
    loglevel = "WARN"

    inputs = [os.path.join(DATA_DIR, "cli/output/bob_1.counts.tsv.gz") for i in range(0, 2)]

    with tempfile.TemporaryDirectory() as tdir:
        output = os.path.join(tdir, "single")
        with pytest.raises(SystemExit):
            merge.merge_counts(output, inputs, None, checksum, loglevel)
