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
import gzip
import json
import logging
import os
import sys
from typing import Dict
from typing import Tuple

from pycroquet.classes import Guide
from pycroquet.classes import Library
from pycroquet.classes import Stats
from pycroquet.constants import COLS_REQ


def _fmt_counts(guide: Guide, inc_unique: bool = True, reverse_sgrna_seqs: bool = False) -> str:
    to_join = []
    for c in COLS_REQ:
        attr = getattr(guide, c)
        if type(attr) is list:
            if c == "sgrna_seqs" and reverse_sgrna_seqs:
                attr = "|".join(reversed(attr))
            else:
                attr = "|".join(attr)
        to_join.append(attr)

    if inc_unique:
        to_join.append(str(int(guide.unique)))  # better for the output

    to_join.append(str(guide.count))
    return "\t".join(to_join)


def _header(sample: str, inc_unique: bool = True) -> str:
    header = "#"
    for c in COLS_REQ:
        header += f"{c}\t"
    if inc_unique:
        header += "unique_guide\t"
    header += f"reads_{sample}"
    return header


def guide_counts_single(
    library: Library,
    guide_results: Dict[str, int],
    output: str,
    stats: Stats,
    low_count: int = None,
) -> Tuple[str, int]:
    """
    Generates the primary result file, reads hitting guides
    """
    if low_count is True:
        stats.low_count_guides_user = {"lt": low_count, "count": 0}
    count_output = f"{output}.counts.tsv.gz"
    logging.info(f"Writing counts file: {count_output}")
    count_total = 0
    # target_to_guides = library.target_to_guides  # can't remember why I included this
    with gzip.open(count_output, "wt") as cout:
        print("##Command: " + stats.command, file=cout)
        print("##Version: " + stats.version, file=cout)
        print(_header(stats.sample_name), file=cout)
        # guide_results = {seq : count,...}
        for guide in library.guides:
            # need to output columns in a stable order, core fields then alpha?
            guide_seq = guide.sgrna_seqs[0]
            count = guide_results.get(guide_seq, 0)
            if count < 30:
                stats.low_count_guides_lt_30 += 1
                if count < 15:
                    stats.low_count_guides_lt_15 += 1
                if count == 0:
                    stats.zero_count_guides += 1
            if low_count is True and count < low_count:
                stats.low_count_guides_user["count"] += 1
            guide.count = count
            print(_fmt_counts(guide), file=cout)
            count_total += count

    stats.mean_count_per_guide = round(count_total / len(library.guides), 2)

    stats_output = f"{output}.stats.json"
    logging.info(f"Writing statistics file: {stats_output}")
    with open(stats_output, "wt") as jout:
        print(stats.as_json(), file=jout)
    return (count_output, count_total)


def query_counts(
    query_dict: Dict[str, int],
    stats: Stats,
    output: str,
):
    """
    This generates a file with the number of incidents of the same query sequence
    """
    count_output = f"{output}.query_counts.tsv.gz"
    logging.info(f"Writing query counts file: {count_output}")
    with gzip.open(count_output, "wt") as cout:
        print("##Command: " + stats.command, file=cout)
        print("##Version: " + stats.version, file=cout)
        print("#QUERY\tCOUNT", file=cout)
        for k in sorted(query_dict.keys()):
            print(f"{k}\t{query_dict[k]}", file=cout)
