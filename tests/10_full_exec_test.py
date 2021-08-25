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
import json
import os
import tempfile
from pprint import pprint

import pytest

from pycroquet import singleguide

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data/cli",
)

"""
This test set targets the pycroquet cli functions in order
to do more complex path validation.
"""


@pytest.mark.parametrize(
    "queries, sample, rules, qual_offset, no_alignment, stats, cram_expected",
    [
        (
            os.path.join(DATA_DIR, "input", "mini.fq.gz"),
            "bob",
            [],
            33,
            False,
            os.path.join(DATA_DIR, "output", "mini_exact.stats.json"),
            True,
        ),
        (
            os.path.join(DATA_DIR, "input", "mini.fq.gz"),
            "bob",
            ["M"],
            33,
            True,
            os.path.join(DATA_DIR, "output", "mini_M1.stats.json"),
            False,
        ),
    ],
)
def test_01_single_guide(queries, sample, rules, qual_offset, no_alignment, stats, cram_expected):
    guidelib = os.path.join(DATA_DIR, "input", "guides.tsv.gz")
    low_count = None
    minscore = 17
    cpus = 1
    chunks = 1000
    exact_mode = "all"
    reference = None
    excludeqcf = False
    stats_old = None
    stats_new = None

    with tempfile.TemporaryDirectory() as tdir:
        singleguide.run(
            guidelib,
            queries,
            sample,
            os.path.join(tdir, "result"),
            os.path.join(tdir, "workspace"),
            rules,
            low_count,
            minscore,
            qual_offset,
            cpus,
            reference,
            excludeqcf,
            chunks,
            no_alignment,
            exact_mode,
            "CRITICAL",
        )

        with open(os.path.join(tdir, "result.stats.json"), "r") as sfp:
            stats_new = json.load(sfp)

        assert os.path.exists(os.path.join(tdir, "result.cram")) == cram_expected

    with open(stats, "r") as sfp:
        stats_old = json.load(sfp)

    for k in ("mapped_to_guide_reads", "total_guides", "total_reads", "unmapped_reads"):
        assert stats_new[k] == stats_old[k], k
