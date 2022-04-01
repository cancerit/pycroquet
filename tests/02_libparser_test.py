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

import pytest

from pycroquet import libparser
from pycroquet.classes import Library
from pycroquet.classes import LibraryHeader

DATA_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
)
CONF_ROOT_KEY_CHECK = ("columns", "headers")


@pytest.mark.parametrize(
    "config, line, lib_type, expected_err",
    [
        ({}, "", None, "Failed to find the column headers in the guide library file"),
        (
            {"columns": {"required": {"id": None}, "optional": {}}},
            "#id",
            None,
            r"Either column 'sgrna_strands' or info header '##library-type' \(single/dual\) must be defined",
        ),
        (
            {"columns": {"required": {"id": None}, "optional": {}}},
            "#",
            None,
            r"^Required column \(.+\) missing from column headers$",
        ),
        # confirm missing works for multiples
        (
            {
                "columns": {
                    "required": {"id": None, "sgrna_ids": None, "sgrna_seqs": None},
                    "optional": {},
                }
            },
            "#id\tsgrna_ids",
            None,
            r"^Required column \(.+\) missing from column headers$",
        ),
        (
            {"columns": {"required": {"id": None}, "optional": {}}},
            "#id\tbob",
            None,
            r"^Column (.+) is not an expected column$",
        ),
    ],
)
def test_01_libparser_columns_from_header_raises(config, line, lib_type, expected_err):
    with pytest.raises(ValueError, match=expected_err) as e_info:
        libparser.columns_from_header(config, line, lib_type)


@pytest.mark.parametrize(
    "config, line, lib_type, e_map, e_cols, e_sep",
    [
        # basic check for cols starting #id, although techincally order is unimportant
        (
            {"columns": {"required": {"id": None}, "optional": {}}},
            "#id",
            "single",
            {0: "id"},
            ["id"],
            {},
        ),
        # confirm we've allowed "# id"
        (
            {"columns": {"required": {"id": None}, "optional": {}}},
            "# id",
            "single",
            {0: "id"},
            ["id"],
            {},
        ),
        # confirm multiple work as expected
        (
            {"columns": {"required": {"id": None, "sgrna_ids": None}, "optional": {}}},
            "#id\tsgrna_ids",
            "single",
            {0: "id", 1: "sgrna_ids"},
            ["id", "sgrna_ids"],
            {},
        ),
        # confirm optional work as expected
        (
            {"columns": {"required": {"id": None}, "optional": {"sgrna_symbols": {}}}},
            "#id\tsgrna_symbols",
            "single",
            {0: "id", 1: "sgrna_symbols"},
            ["id", "sgrna_symbols"],
            {},
        ),
        # separators work
        (
            {
                "columns": {
                    "required": {"id": None, "sgrna_ids": {"separator": "|"}},
                    "optional": {"sgrna_symbols": {"separator": "|"}},
                }
            },
            "#id\tsgrna_ids\tsgrna_symbols",
            "single",
            {0: "id", 1: "sgrna_ids", 2: "sgrna_symbols"},
            ["id", "sgrna_ids", "sgrna_symbols"],
            {"sgrna_ids": "|", "sgrna_symbols": "|"},
        ),
    ],
)
def test_02_libparser_columns_from_header(config, line, lib_type, e_map, e_cols, e_sep):
    # expected format for start of column line
    (cols, map, col_sep) = libparser.columns_from_header(config, line, lib_type)
    assert map == e_map
    assert cols == e_cols
    assert col_sep == e_sep


def test_03_libparser_load_resource_config():
    conf = libparser._load_config()
    assert type(conf) is dict, "Expected data type from resource config yaml"
    for k in CONF_ROOT_KEY_CHECK:
        assert k in conf, "Root config entries found"


@pytest.mark.parametrize(
    "file, expected_err",
    [
        ("bad_library_info_header.tsv", r"^Header info line is not of expected format"),
        (
            "bad_library_type.tsv",
            r"^Value for 'library-type' \(.+\) is not valid, choose from: ",
        ),
        (
            "bad_library_info_header_duplicates.tsv",
            r"^Duplicate key '.+' found for '##' header line$",
        ),
    ],
)
def test_04_libparser_parse_header_raises(file, expected_err):
    with pytest.raises(ValueError, match=expected_err) as e_info:
        with open(os.path.join(DATA_DIR, file)) as ifh:
            line = ifh.readline()
            libparser.parse_header(ifh, line)


def test_05_libparser_parse_header():
    with open(os.path.join(DATA_DIR, "good_library_header.tsv")) as ifh:
        line = ifh.readline()
        lh = libparser.parse_header(ifh, line)
        assert type(lh) is LibraryHeader


@pytest.mark.parametrize(
    "file, expected_err",
    [
        (
            "bad_guide_01.tsv",
            r"^Column header indicates \d+ columns, but record has only \d+",
        ),
        (
            "bad_guide_02.tsv",
            r"^'sgrna_strands' column elements \(\d+\) mismatch vs rule for '##library-type:",
        ),
        (
            "bad_guide_03.tsv",
            r"^'sgrna_strands' column elements \(\d+\) mismatch against previous rows \(\d+\)",
        ),
    ],
)
def test_06_libparser_parse_data_rows_raises(file, expected_err):
    with pytest.raises(ValueError, match=expected_err) as e_info:
        with open(os.path.join(DATA_DIR, file)) as ifh:
            line = ifh.readline()
            lh = libparser.parse_header(ifh, line)
            libparser.parse_data_rows(lh, ifh)


@pytest.mark.parametrize(
    "file, e_guides",
    [
        (
            "good_guide_01.tsv",
            2,
        ),
        (
            "good_guide_02.tsv",
            2,
        ),
        (
            "good_guide_03.tsv",
            2,
        ),
        (
            "good_guide_04.tsv",
            2,
        ),
    ],
)
def test_07_libparser_parse_data_rows(file, e_guides):
    with open(os.path.join(DATA_DIR, file)) as ifh:
        line = ifh.readline()
        lh = libparser.parse_header(ifh, line)
        guides = libparser.parse_data_rows(lh, ifh)
        assert len(guides) == e_guides


def test_08_libparser_parse_data_rows_other():
    with open(os.path.join(DATA_DIR, "good_guide_05.tsv")) as ifh:
        line = ifh.readline()
        lh = libparser.parse_header(ifh, line)
        guides = libparser.parse_data_rows(lh, ifh)
        assert len(guides) == 2
        assert "custom_annotation" in guides[0].other
        assert guides[0].other["custom_annotation"] == "comment 1"
        assert "custom_annotation" in guides[1].other
        assert guides[1].other["custom_annotation"] == "comment 2"
        assert type(guides[0].sgrna_seqs) is list
        assert type(guides[0].sgrna_strands) is list
        assert type(guides[1].sgrna_seqs) is list
        assert type(guides[1].sgrna_strands) is list


def test_09_libparser_parse_data_rows_badsplits():
    with pytest.raises(
        ValueError,
        match="^All columns that can be split should have the same number of elements.  Column.+",
    ) as e_info:
        with open(os.path.join(DATA_DIR, "bad_guide_04.tsv")) as ifh:
            line = ifh.readline()
            lh = libparser.parse_header(ifh, line)
            guides = libparser.parse_data_rows(lh, ifh)


@pytest.mark.parametrize(
    "file, error",
    [
        (
            "bad_guide_05.tsv",
            r"'sgrna_seqs' can only contain ACGT and the separator character '\|'.  Got 'acgt' after splitting",
        ),
        (
            "bad_guide_06.tsv",
            r"'sgrna_seqs' can only contain ACGT and the separator character '\|'.  Got '' after splitting",
        ),
    ],
)
def test_10_libparser_parse_data_rows_badseq(file, error):
    with pytest.raises(ValueError, match=error) as e_info:
        with open(os.path.join(DATA_DIR, file)) as ifh:
            line = ifh.readline()
            lh = libparser.parse_header(ifh, line)
            guides = libparser.parse_data_rows(lh, ifh)


@pytest.mark.parametrize(
    "file, info",
    [
        ("good_guide_01.tsv", "text file"),
        ("good_guide_01.tsv.gz", "gzip text file"),
    ],
)
def test_20_libparser_load(file, info):
    assert type(libparser.load(os.path.join(DATA_DIR, file))) is Library, info
