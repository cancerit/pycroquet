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
import pytest

from pycroquet.classes import Guide
from pycroquet.classes import LibraryHeader
from pycroquet.targets import guides_to_targets

GUIDE_A = Guide(idx=0, sgrna_seqs=["AAAA"])
GUIDE_C = Guide(idx=0, sgrna_seqs=["CCCC"])
GUIDE_G = Guide(idx=0, sgrna_seqs=["GGGG"])
GUIDE_T = Guide(idx=0, sgrna_seqs=["TTTT"])
GUIDE_AC = Guide(idx=0, sgrna_seqs=["AAAA", "CCCC"])


@pytest.mark.parametrize(
    "guides, single, result, info",
    [
        ([GUIDE_A], True, ({"AAAA": [0]}, ["AAAA"]), "Basic validation"),
        (
            [GUIDE_A, GUIDE_C],
            True,
            ({"AAAA": [0], "CCCC": [1]}, ["AAAA", "CCCC"]),
            "2 unique records",
        ),
        ([GUIDE_A, GUIDE_A], True, ({"AAAA": [0, 1]}, ["AAAA"]), "Duplicate"),
        (
            [GUIDE_AC, GUIDE_C, GUIDE_T],
            False,
            ({"AAAA": [0], "CCCC": [0, 1], "TTTT": [2]}, ["AAAA", "CCCC", "TTTT"]),
            "Dual with clashes",
        ),
    ],
)
def test_01_targets_guides_to_targets(guides, single, result, info):
    assert guides_to_targets(guides, single) == result, info
