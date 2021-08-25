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
import logging
from typing import Dict
from typing import List
from typing import Tuple

from pygas.matrix import revcomp

from pycroquet.classes import Guide

"""
Handles minimizing target sequences and the de-convolution back to the guide sets they came from
Overall this is is easy for singe guide but more complex for multi-guide
"""


def guides_to_targets(guides: List[Guide], single: bool = False) -> Tuple[Dict[str, List[int]], List[str]]:
    # this links a target sequence to the guide library position(s) in the list
    # target orientation is that of first added, we then check for equivalent targets
    target_to_guides = {}  # Dict[str, List[int]]
    total_dups = 0
    for g_idx, guide in enumerate(guides):
        for target in guide.sgrna_seqs:
            if target in target_to_guides:
                if single:
                    logging.debug(f"Guide duplicate: {guide.id}")
                    total_dups += 1
                    # this only works for single guide
                    guide.unique = False
                    guides[target_to_guides[target][0]].unique = False
                target_to_guides[target].append(g_idx)
            else:
                target_to_guides[target] = [g_idx]

    # now we need the unique list of targets, in a stable order
    uniq_targets = list(target_to_guides.keys())
    uniq_targets.sort()
    if single:
        logging.info(f"Number of duplicate guides: {total_dups}")
    logging.info(f"Total unique guides: {len(uniq_targets)}")
    return (target_to_guides, uniq_targets)
