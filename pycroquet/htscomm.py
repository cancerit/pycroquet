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
"""
For common hts file actions
"""
import logging
import os
from tempfile import TemporaryDirectory

import pysam


def hts_reader(seq_file, mode, cpus, reference=None) -> pysam.AlignmentFile:
    hts_cpus = cpus if cpus < 4 else 4
    save = pysam.set_verbosity(0)
    sam = pysam.AlignmentFile(
        seq_file,
        mode=mode,
        reference_filename=reference,
        require_index=False,
        threads=hts_cpus,
        check_sq=False,
    )
    pysam.set_verbosity(save)
    return sam


def hts_sort_n_index(unsorted, guide_fa, output_stub, workspace, cpus=1):
    sam_threads = 4
    if cpus < sam_threads:
        sam_threads = cpus
    aligned_cram = f"{output_stub}.cram"
    logging.info(f"Sorting alignments to: {aligned_cram}")
    with TemporaryDirectory(prefix="sort_to_cram", dir=workspace) as td:
        pysam.sort(
            "--no-PG",
            "-@",
            str(sam_threads),
            "--reference",
            guide_fa,
            "--output-fmt",
            "CRAM,no_ref=1",
            "-T",
            os.path.join(td, "samsort"),
            "-o",
            aligned_cram,
            unsorted,
        )
    logging.info(f"Generating alignment index: {aligned_cram}.crai")
    pysam.index(aligned_cram)
