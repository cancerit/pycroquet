#!/usr/bin/env python3
#
# Copyright (c) 2022
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
import re
import sys

INCREMENT = 1

sgrna_id_pattern = re.compile("([^:]+):(\d+)-(\d+)_([+-])")


def sgrna_id_inc(sgrna_id: str):
    (chr, start, stop, strand) = re.search(sgrna_id_pattern, sgrna_id).group(1, 2, 3, 4)
    start = int(start) + INCREMENT
    stop = int(stop) + INCREMENT
    return f"{chr}:{str(start)}-{str(stop)}_{strand}"


if len(sys.argv) < 3:
    print(f"USAGE: {sys.argv[0]} library.tsv L|D", file=sys.stderr)
    print(f"\tWhere:", file=sys.stderr)
    print(f"\t\t L = Output list of bad ids", file=sys.stderr)
    print(f"\t\t D = Output data for bad ids", file=sys.stderr)
    sys.exit(1)

mode = sys.argv[2]

lookup = {}
with open(sys.argv[1]) as ifh:
    for line in ifh:
        if line.startswith("#"):
            continue
        data = line.strip().split("\t")
        lookup[data[1]] = data

all_bad_paralogues = {}

for sgrna_ids in lookup.keys():
    # get the item and try +1 on each end to identify "mates"
    (lhs, rhs) = sgrna_ids.split("|")
    related = []
    alt_lhs = f"{sgrna_id_inc(lhs)}|{rhs}"
    if alt_lhs in lookup:
        related.append((lookup[alt_lhs][3], lookup[alt_lhs][0], alt_lhs))
    alt_rhs = f"{lhs}|{sgrna_id_inc(rhs)}"
    if alt_rhs in lookup:
        related.append((lookup[alt_rhs][3], lookup[alt_rhs][0], alt_rhs))
    if len(related):
        if mode == "D":
            print(f"{lookup[sgrna_ids][3]}\t{lookup[sgrna_ids][0]}\t{sgrna_ids}", end="")
        all_bad_paralogues[lookup[sgrna_ids][0]] = 1
        for t in related:
            all_bad_paralogues[t[0]] = 1
            if mode == "D":
                print(f"\t{t[0]}\t{t[1]}\t{t[2]}", end="")
        if mode == "D":
            print()

if mode == "L":
    print("\n".join(sorted(all_bad_paralogues.keys())))
