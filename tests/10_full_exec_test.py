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

from pycroquet import dualguide
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


"""
pycroquet dual-guide -g tests/data/cli/input/dual_lib.tsv -q tests/data/cli/input/dual_reads.bam -o tests/data/cli/output/dual_nolow -b exact -w work
pycroquet dual-guide -g tests/data/cli/input/dual_lib.tsv -q tests/data/cli/input/dual_reads.bam -o tests/data/cli/output/dual_low_count -b exact -w work --low_count 5
"""


@pytest.mark.parametrize(
    "low_count, compare_to",
    [
        (None, "dual_nolow"),
        (5, "dual_low_count"),
    ],
)
def test_02_dual_guide(low_count, compare_to):
    guidelib = os.path.join(DATA_DIR, "input", "dual_lib.tsv.gz")
    queries = os.path.join(DATA_DIR, "input", "dual_reads.bam")
    sample = None
    rules = []
    # low_count = None
    minscore = 15
    qual_offset = None
    cpus = 1
    reference = None
    excludeqcf = False
    boundary_mode = "TinQ"
    trimseq = 0
    chunks = 1000
    loglevel = "WARN"

    stats_new = None

    with tempfile.TemporaryDirectory() as tdir:
        output = os.path.join(tdir, "result")
        workspace = os.path.join(tdir, "workspace")
        dualguide.run(
            guidelib,
            queries,
            sample,
            output,
            workspace,
            rules,
            low_count,
            minscore,
            qual_offset,
            cpus,
            reference,
            excludeqcf,
            boundary_mode,
            trimseq,
            chunks,
            loglevel,
        )
        out_counts = f"{output}.counts.tsv.gz"
        out_stats = f"{output}.stats.json"
        out_cram = f"{output}.cram"
        out_crai = f"{output}.cram.crai"

        assert os.path.exists(out_counts)
        assert os.path.exists(out_stats)
        assert os.path.exists(out_cram)
        assert os.path.exists(out_crai)

        with open(out_stats, "r") as sfp:
            stats_new = json.load(sfp)

    stats_old = None
    with open(os.path.join(DATA_DIR, "output", f"{compare_to}.stats.json"), "r") as sfp:
        stats_old = json.load(sfp)

    for i in ("command", "version"):
        del stats_old[i]
        del stats_new[i]
    assert stats_new == stats_old
