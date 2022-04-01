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
import pathlib
import shutil
import tempfile
import types

import magic

from pycroquet import tools


def test_01_tools_chunks():
    newlist = tools.chunks([1, 2, 3], 2)
    assert type(newlist) is types.GeneratorType
    assert list(newlist) == [[1, 2], [3]]


def test_02_tools_pickle():
    original_data = ["hello", "world"]
    with tempfile.TemporaryDirectory() as tdir:
        path = tools.pickle_this(tdir, "1", original_data)
        assert magic.from_file(path).startswith("gzip compressed data")
        assert tools.unpickle(path) == original_data
