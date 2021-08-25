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
import os
import sys
from array import array
from dataclasses import dataclass
from dataclasses import field
from typing import Dict
from typing import Final
from typing import List

import pkg_resources


@dataclass
class Classification:
    match: Final = "match"  # same vector F/R
    aberrant_match: Final = "aberrant_match"  # same vector pair, aberrant orientation
    f_multi_3p: Final = "f_multi_5p"  # 5p mapped F, 3p multihit
    f_multi_5p: Final = "f_multi_5p"  # 3p mapped F, 5p multihit
    r_multi_3p: Final = "r_multi_5p"  # 5p mapped R, 3p multihit
    r_multi_5p: Final = "r_multi_5p"  # 3p mapped R, 5p multihit
    f_open_3p: Final = "f_open_3p"  # 5p mapped F, 3p open (unmapped)
    f_open_5p: Final = "f_open_5p"  # 3p mapped F, 5p open (unmapped)
    r_open_3p: Final = "f_open_3p"  # 5p mapped R, 3p open (unmapped)
    r_open_5p: Final = "f_open_5p"  # 3p mapped R, 5p open (unmapped)
    swap: Final = "swap"  # multi vector, uniq mapped
    ambiguous: Final = "ambiguous"  # both ends multi hit
    no_match: Final = "no_match"  # multi/unmapped either end


"""
Feature count total_number
MATCH: 18379264 20617796  // The match is perfect R1 F(Forward) and R2 R(reverse)
MATCHFF: 0 20617796.      // Match but R1 is F and R2 is F
MATCHRF: 0 20617796       // Match but R1 is R and R2 is F
MATCHRR: 0 20617796       // Match but R1 is R and R2 is R
SWAPFF: 0 20617796        // R1 is mapping a vector F and R2 is mapping another vector F
SWAPFR: 404446 20617796   // R1 is mapping a vector F and R2 is mapping another vector R
SWAPRF: 0 20617796        // R1 is mapping a vector R and R2 is mapping another vector F
SWAPRR: 0 20617796        // R1 is mapping a vector R and R2 is mapping another vector R
OEM1F: 1063485 20617796  //  R1 is mapping to a vector F and R2 is not mapping
OEM1R: 0 20617796        //  R1 is mapping to a vector R and R2 is not mapping
OEM2F: 0 20617796        //  R2 is mapping to a vector F and R1 is not mapping
OEM2R: 722427 20617796   //  R2 is mapping to a vector R and R1 is not mapping
NOMATCH: 48174 20617796 // Both R1 and R2 are not mapping
"""


@dataclass
class LibraryHeader:
    column_map: dict
    column_list: List[str]
    column_separators: dict
    required_cols: List[str]
    config: dict
    info_items: dict = None
    is_single: bool = False


@dataclass
class Guide:
    idx: int
    id: str = None
    sgrna_ids: List[str] = None
    sgrna_seqs: List[str] = None
    gene_pair_id: str = None
    sgrna_strands: List[int] = None
    other: dict[str, str] = field(default_factory=dict)
    unique: bool = True
    count: int = 0


@dataclass
class Library:
    """
    header: the header object
    guides: guide detail objects
    targets: unique guide sequences
    target_to_guides: mappings of target sequences back to guides
    """

    header: LibraryHeader
    guides: List[Guide]
    targets: List[str]
    target_to_guides: Dict[str, List[int]]
    _sgrna_ids_by_seq: Dict[str, int] = None
    _guide_by_sgrna_set: Dict[str, int] = None

    def min_target_len(self) -> int:
        return len(min(self.targets, key=len))

    def sgrna_ids_by_seq(self, target_seq):
        """
        get the target by the 0,1... position in a guide
        """
        if self._sgrna_ids_by_seq is None:
            sgrna_ids_by_seq = {}
            frags = len(self.guides[0].sgrna_ids)
            for g in self.guides:
                for i in range(0, frags):
                    g_seq = g.sgrna_seqs[i]
                    g_id = g.sgrna_ids[i]
                    if g_seq not in sgrna_ids_by_seq:
                        sgrna_ids_by_seq[g_seq] = {}
                    sgrna_ids_by_seq[g_seq][g_id] = None
            # reduce to list of keys, note data type changes from dict to list
            for g_seq in sgrna_ids_by_seq:
                sgrna_ids_by_seq[g_seq] = sorted(sgrna_ids_by_seq[g_seq].keys())
            self._sgrna_ids_by_seq = sgrna_ids_by_seq
        return self._sgrna_ids_by_seq[target_seq]

    def guide_by_sgrna_set(self, seq_l, seq_r):
        match = f"{seq_l}|{seq_r}"
        if self._guide_by_sgrna_set is None:
            data = {}
            for i, g in enumerate(self.guides):
                data["|".join(g.sgrna_seqs)] = i
            self._guide_by_sgrna_set = data
        return self._guide_by_sgrna_set.get(match)


@dataclass
class Stats:
    """
    Read and guide statistics
    """

    total_reads: int = 0
    total_pairs: int = 0
    vendor_failed_reads: int = 0
    mapped_to_guide_reads: int = 0
    # reads that map equally well to multiple guides, mainly for SG data
    multimap_reads: int = 0
    unmapped_reads: int = 0
    total_guides: int = 0
    zero_count_guides: int = 0
    low_count_guides_lt_15: int = 0
    low_count_guides_lt_30: int = 0
    length_excluded_reads: int = None
    low_count_guides_user: Dict[str, int] = None
    reversed_reads: bool = False
    version: str = None
    command: str = None
    sample_name: str = None
    pair_classifications: Dict[str, int] = None

    def __post_init__(self):
        self.command = ""
        for i, e in enumerate(sys.argv):
            if i == 0:
                self.command += f"{os.path.basename(e)}"
                continue
            self.command += f" {e}"
        self.version = pkg_resources.require(__name__.split(".")[0])[0].version


@dataclass
class Seqread:
    """
    Raw sequence info required for mapping outputs
    """

    qname: str
    sequence: str
    qual: array
    rgid: str
    qc_fail: bool = False
    member: int = None
