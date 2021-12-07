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
import gzip
import os
import tempfile

import pytest

from pycroquet import countwriter
from pycroquet.classes import Stats
from pycroquet.libparser import load

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data", "countwriter")


@pytest.mark.parametrize(
    "guide_set, stats, exp_count, info",
    [
        ({"AAAAAACCCCCC": 1}, Stats(total_reads=1), 1, "forward"),
        ({"GGGGGGTTTTTT": 1}, Stats(total_reads=1), 0, "reversed"),
        ({}, Stats(), 0, "nohit"),
    ],
)
def test_01_countwriter_guide_counts_single(guide_set, stats, exp_count, info):
    library = load(os.path.join(DATA_DIR, "lib.tsv"))
    with tempfile.TemporaryDirectory() as tdir:
        stub = os.path.join(tdir, "result.tsv")
        (output, total_count) = countwriter.guide_counts_single(library, guide_set, stub, stats=stats)
        assert len(gzip.open(output, "rt").readlines()) == 4, info  # cmd, version, header, result
        assert total_count == exp_count, info
