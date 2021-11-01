#!/usr/bin/env python3
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
from setuptools import setup

config = {
    "name": "pycroquet",
    "description": "python Crispr Read to Oligo QUantification Enhancement Tool",
    "long_description": open("README.md").read(),
    "long_description_content_type": "text/markdown",
    "author": "Keiran M Raine",
    "url": "https://github.com/cancerit/pycroquet",
    "author_email": "cgphelp@sanger.ac.uk",
    "version": "1.2.1",
    "license": "AGPL-3.0",
    "python_requires": ">= 3.9",
    "install_requires": ["click", "click-option-group", "python-magic", "pysam", "pygas", "PyYAML"],
    "packages": ["pycroquet"],
    "setup_requires": ["click"],
    "test_suite": "tests",
    "tests_require": ["pytest"],
    "entry_points": {
        "console_scripts": ["pycroquet=pycroquet.cli:cli"],
    },
    "package_data": {
        "pycroquet": [
            "resources/library.yaml",
        ],
    },
}

setup(**config)
