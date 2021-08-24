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
import re
from typing import Dict
from typing import List
from typing import TextIO
from typing import Tuple

import yaml
from pkg_resources import resource_string

from pycroquet.classes import Guide
from pycroquet.classes import Library
from pycroquet.classes import LibraryHeader
from pycroquet.targets import guides_to_targets

"""
See: https://docs.google.com/document/d/18G5kT0Z6gmrUwEq6JZ-sd0XuZf7iG5WsDLBY2Bb_uLE/edit

This document needs pulling into the README.

Flexible file definition, some core fields and then additional columns... but that means we need a dynamic parser.
"""

H_INFO_RE = re.compile("^##([^:]+): ?(.+)$")
ACGT_ONLY = re.compile("[ACGT]+")
LIBRARY_TYPES = ("single", "dual", "other")
LIBRARY_STRAND_DEFAULTS = {"single": (0), "dual": (0, 1)}


def _load_config():
    return yaml.safe_load(resource_string(__name__, f"resources/library.yaml").decode("utf-8", "strict"))


def parse_header(ifh: TextIO, line: str) -> LibraryHeader:
    """
    Expects file handle, and first line of file (to cope with open vs gz.open).
    Capture any info header items '##', validate any we decide to be required.
    Finally pass the column '#' row for parsing, do not read another line following this.
    """

    config = _load_config()

    info_items = {}
    lib_type = None
    while line and line.startswith("##"):
        # these are the info lines
        line = line.rstrip()
        m_res = H_INFO_RE.match(line)
        if m_res is None:
            raise ValueError("Header info line is not of expected format '^##([^:]+): ?(.+)$': " + line)
        (tag, value) = m_res.group(1, 2)
        if tag in info_items:
            raise ValueError(f"Duplicate key '{tag}' found for '##' header line")
        if tag == "library-type":
            if value not in LIBRARY_TYPES:
                raise ValueError(f"Value for 'library-type' ({value}) is not valid, choose from: " + str(LIBRARY_TYPES))
            lib_type = value
        info_items[tag] = value
        line = ifh.readline()
    (cols, col_idx_map, col_sep) = columns_from_header(config, line, lib_type)
    if "library-type" not in info_items:
        info_items["library-type"] = "UNKNOWN"
    return LibraryHeader(
        column_map=col_idx_map,
        column_list=cols,
        info_items=info_items,
        required_cols=config["columns"]["required"],
        config=config,
        column_separators=col_sep,
    )


def columns_from_header(config, line: str, lib_type: str = None) -> Tuple[List, Dict[int, str], Dict[str, str]]:
    """
    this should have single leading #
    """
    line = line.rstrip()
    if not line.startswith("#"):
        raise ValueError("Failed to find the column headers in the guide library file")
    # map columns to index positions so easy to map data into object/dicts
    idx_to_col = {}
    # also capture the list of headings seen
    cols = []
    # also remove any left whitespace after the # just to minimise errors
    for i, c in enumerate(line.lstrip("#").lstrip().split("\t")):
        idx_to_col[i] = c
        cols.append(c)

    conf_cols = config["columns"]
    required_keys = list(conf_cols["required"].keys())
    optional_keys = list(conf_cols["optional"].keys())

    for required in required_keys:
        if required not in cols:
            raise ValueError(f"Required column ({required}) missing from column headers")

    for c in cols:
        if c not in required_keys + optional_keys:
            raise ValueError(f"Column ({c}) is not an expected column")

    if "sgrna_strands" not in cols and (lib_type is None or lib_type not in ("single", "dual")):
        raise ValueError("Either column 'sgrna_strands' or info header '##library-type' (single/dual) must be defined")

    split_cols = {}
    for c in cols:
        if c in required_keys and conf_cols["required"][c] is not None and "separator" in conf_cols["required"][c]:
            split_cols[c] = conf_cols["required"][c]["separator"]
            continue
        if c in optional_keys and conf_cols["optional"][c] is not None and "separator" in conf_cols["optional"][c]:
            split_cols[c] = conf_cols["optional"][c]["separator"]

    return (cols, idx_to_col, split_cols)


def _validate_strands(lib_type, this_no, strand_no, line):
    if lib_type != "UNKNOWN":
        raise ValueError(
            f"'sgrna_strands' column elements ({this_no}) mismatch vs rule for '##library-type: {lib_type}' ({strand_no})"
        )
    else:
        raise ValueError(
            f"'sgrna_strands' column elements ({this_no}) mismatch against previous rows ({strand_no}), line: {line}"
        )


def _strand_info(lh: LibraryHeader) -> Tuple[int, Tuple[int]]:
    strand_no = None
    strand_set = None
    if lh.info_items["library-type"] == "single":
        strand_no = 1
        strand_set = LIBRARY_STRAND_DEFAULTS[lh.info_items["library-type"]]
    elif lh.info_items["library-type"] == "dual":
        strand_no = 2
        strand_set = LIBRARY_STRAND_DEFAULTS[lh.info_items["library-type"]]
    return (strand_no, strand_set)


def parse_data_rows(lh: LibraryHeader, ifh: TextIO) -> List[Guide]:
    col_list = lh.column_list
    col_no = len(col_list)
    col_map = lh.column_map
    required_cols = lh.required_cols
    (strand_no, strand_set) = _strand_info(lh)
    idx = 0
    guides = []

    cols_to_split = lh.column_separators

    while True:
        line = ifh.readline()
        if line == "":
            break
        elements = line.strip().split("\t")
        if len(elements) != col_no:
            raise ValueError(
                f"Column header indicates {col_no} columns, but record has only {len(elements)} (all fields require a value, '.' to omit)\n> {line}"
            )
        guide = Guide(idx=idx)
        for pos, key in col_map.items():
            value = elements[pos]
            if key in cols_to_split:
                # autosplit
                value = value.split(cols_to_split[key])
            if key in required_cols or key == "sgrna_strands":
                setattr(guide, key, value)
            else:
                guide.other[key] = value

        if guide.sgrna_strands is None:
            guide.sgrna_strands = strand_set
        else:
            this_no = len(guide.sgrna_strands)
            if strand_no is None:
                strand_no = this_no
            elif strand_no != this_no:
                _validate_strands(lh.info_items["library-type"], this_no, strand_no, line)

        expect_len = None
        initial_col = None
        for col in cols_to_split:
            if col in guide.other:
                split_len = len(guide.other[col])
            else:
                split_len = len(getattr(guide, col))
            if expect_len is None:
                expect_len = split_len
                initial_col = col
            elif expect_len != split_len:
                raise ValueError(
                    f"All columns that can be split should have the same number of elements.  Column {col} ({split_len}), differs from {initial_col} ({expect_len})"
                )

        for seq in guide.sgrna_seqs:
            if ACGT_ONLY.fullmatch(seq) is None:
                print(guide)
                raise ValueError(
                    f"'sgrna_seqs' can only contain ACGT and the separator character '|'.  Got '{seq}' after splitting"
                )

        guides.append(guide)
        idx += 1

    if strand_no == 1 and lh.info_items["library-type"] == "UNKNOWN":
        lh.info_items["library-type"] = "single"

    if lh.info_items["library-type"] == "single":
        lh.is_single = True

    return guides


def load(library_file: str) -> Library:
    i_fh = None
    line = None
    try:
        i_fh = gzip.open(library_file, "rt")
        line = i_fh.readline()
    except gzip.BadGzipFile:
        i_fh = open(library_file, "rt")
        line = i_fh.readline()
    lh = parse_header(i_fh, line)
    guides = parse_data_rows(lh, i_fh)
    i_fh.close()

    (target_to_guides, uniq_targets) = guides_to_targets(guides, lh.is_single)
    return Library(
        header=lh,
        guides=guides,
        targets=uniq_targets,
        target_to_guides=target_to_guides,
    )
