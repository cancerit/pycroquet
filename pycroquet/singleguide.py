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
import shutil

from pycroquet import cli
from pycroquet import countwriter
from pycroquet import libparser
from pycroquet import main
from pycroquet import readwriter


def run(
    guidelib,
    queries,
    sample,
    output,
    workspace,
    rules,
    low_count,
    minscore,
    qual_offset,
    cpus,
    reference,
    excludeqcf,
    chunks,
    no_alignment,
    boundary_mode,
    loglevel,
):
    (usable_cpu, work_tmp, workspace, boundary_mode) = cli.common_setup(
        loglevel, cpus, workspace, output, boundary_mode
    )

    library = libparser.load(guidelib)
    reverse = False
    (_, guide_results, aligned_results, stats) = main.process_reads(
        library,
        queries,
        workspace,
        rules,
        minscore,
        usable_cpu,
        chunks,
        sample=sample,
        reference=reference,
        exclude_qcfail=excludeqcf,
        reverse=reverse,
        boundary_mode=boundary_mode,
    )

    countwriter.guide_counts_single(library, guide_results, output, stats, low_count)
    if no_alignment is False:
        readwriter.reads_to_hts(
            library,
            aligned_results,
            queries,
            qual_offset,
            workspace,
            stats,
            output,
            usable_cpu,
            exclude_qcfail=excludeqcf,
            reference=reference,
            reverse=reverse,
        )

    # TODO: write everything out and then if mapped count is a very low fraction repeat.  If the revcomp result is a higher
    # fraction then replace data

    if not work_tmp:
        shutil.rmtree(workspace)
