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

import pytest
from pygas.alignercpu import AlignerCpu
from pygas.classes import Backtrack

import pycroquet.tools as ctools
from pycroquet import dualguide
from pycroquet import libparser
from pycroquet import readparser
from pycroquet.classes import Classification
from pycroquet.classes import Library

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
    "dualguide",
)


@pytest.fixture()
def map_resource(request):
    library = libparser.load(os.path.join(DATA_DIR, "guides.tsv"))
    (unique, stats, reads, read_pairs) = readparser.parse_reads(
        os.path.join(DATA_DIR, "reads.sam"),
        sample="bob",
        cpus=1,
        paired=True,
    )
    aligner = AlignerCpu(
        targets=library.targets,
        rules=[],
        score_min=5,
        rev_comp=True,  # as some reads can be reversed in DG
    )
    # removed all the "at-scale" wrapping
    results = aligner.align_queries(reads.keys(), keep_matrix=False)
    # for code path, pickle
    tdir = tempfile.TemporaryDirectory()
    pickles = [ctools.pickle_this(tdir.name, "pre_matrix_{:05d}".format(1), [results])]
    (aligned_results, multi_map, unique_map, unmap) = dualguide.pickles_to_mapset(pickles, reads, aligner)
    yield (aligned_results, multi_map, unique_map, unmap, library, stats)


@pytest.mark.parametrize(
    "read_l, read_r, mtype_l, mtype_r, exp_class, info",
    [
        (
            "GGGGCCCCCC",
            "TTTTAAAAAA",
            "unique",
            "unique",
            Classification.match,
            "unique f/r",
        ),
        (
            "AAAAAAAAAA",
            "GGGGGGGGGG",
            "unique",
            "unique",
            Classification.aberrant_match,
            "unique f/f",
        ),
    ],
)
def test_01_read_pair_to_guides(map_resource, read_l, read_r, mtype_l, mtype_r, exp_class, info):
    (aligned_results, multi_map, unique_map, unmap, library, stats) = map_resource
    # assert stats.total_reads == 4
    # assert stats.total_pairs == 2
    assert aligned_results[read_l][0] == mtype_l
    if mtype_l == "unique":
        assert len(aligned_results[read_l][1]) == 1
    elif mtype_l == "multimap":
        assert len(aligned_results[read_l][1]) > 1
    else:  # unmapped
        print(aligned_results[read_l])
        assert 1 == 2
    assert aligned_results[read_r][0] == mtype_r
    if mtype_r == "unique":
        assert len(aligned_results[read_r][1]) == 1
    elif mtype_r == "multimap":
        assert len(aligned_results[read_r][1]) > 1
    else:  # unmapped
        print(aligned_results[read_r])
        assert 1 == 2
    (classified, gidx, bt_l, br_r, orig_l, orig_r) = dualguide.classify_read_pair(
        aligned_results[read_l], aligned_results[read_r], library
    )
    assert classified == exp_class
