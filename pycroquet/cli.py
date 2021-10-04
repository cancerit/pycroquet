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
import os
import sys
import tempfile
from functools import wraps

import click
import pkg_resources  # part of setuptools
from click_option_group import OptionGroup

from pycroquet import dualguide
from pycroquet import libparser
from pycroquet import main
from pycroquet import readwriter
from pycroquet import sge as pysge
from pycroquet import singleguide
from pycroquet import tools as ctools

LOG_LEVELS = ("WARNING", "INFO", "DEBUG")

HELP_TRIMSEQ = "Trim reads back to use first N bases only. This is a destructive process, the alignment CRAM will not include the full read sequence."
HELP_GUIDELIB = "Expanded guide library definition tsv file with optional headers (common format for single/dual/other)"
HELP_QUERIES = "Query sequence file (fastq[.gz], sam, bam, cram)"
HELP_SAMPLE = (
    "Sample name to apply to count column, required for fastq, interrogate header for others when not defined."
)
HELP_OUTPUT = "Final output to this filename prefix"
HELP_WORKSPACE = "Workspace for temporary files"
HELP_RULES = "Rules for allowing fuzzy matching (no rules = fast exact/substr match only), e.g. MM=2 mismatched bases, MDI allows 1 base of mismatch, deletion and insertion."
HELP_LOW_COUNT = (
    "*.stats.json includes low_count_guides_lt_{15,30}, this option allow specification of an additional cut-off."
)
HELP_QUAL_OFFSET = "Specify phread offset (for fastq) if detection by readname fails"
HELP_CPUS = "CPUs to use (0 to detect)"
HELP_SGE_UNIQUE = "Only generate the unique sequence counts file, then exit (--guide can be omitted)"
HELP_CHUNKS = "Reads per mapping block"
HELP_MINSCORE = "Minimum score to retain, regardless of rule penalties.  Perfect match has score equal to query length."
HELP_REFERENCE = "Required for cram"
HELP_EXCLUDEQCF = "Exclude reads that fail QC (sam/bam/cram, CASAVA fastq only)"
HELP_NO_ALIGNMENT = "Do not output cram alignments"
HELP_BOUNDARY = "Control boundary matching types, see end of options"
HELP_FASTA = "Write fasta to this file"

HELP_EPILOG = """
Additional option info:

-b/--boundary-mode

\b
  all   = all combinations
  exact = requires all guides and queries to be same length
  TinQ  = target in query
  QinT  = query in target
"""


def _log_setup(loglevel):  # pragma: no cover
    sh_info = logging.StreamHandler(stream=sys.stdout)
    sh_info.setLevel(logging.DEBUG)
    sh_error = logging.StreamHandler(stream=sys.stderr)
    sh_error.setLevel(logging.ERROR)

    logging.basicConfig(
        level=getattr(logging, loglevel.upper()), format="%(levelname)s: %(message)s", handlers=[sh_info, sh_error]
    )


def common_setup(loglevel, cpus, workspace, output, boundary_mode):
    _log_setup(loglevel)
    usable_cpu = cpus if cpus > 0 else len(os.sched_getaffinity(0))
    work_tmp = None
    if workspace:
        os.makedirs(workspace, exist_ok=True)
    else:
        work_tmp = tempfile.TemporaryDirectory()
        workspace = work_tmp.name
    outfolder = os.path.dirname(os.path.abspath(output))
    if not os.path.exists(outfolder):
        os.makedirs(outfolder, exist_ok=True)
    boundary_mode = ctools.boundary_mode(boundary_mode)
    return (usable_cpu, work_tmp, workspace, boundary_mode)


def _file_exists():
    return click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    )


def chunk_default(f):
    @click.option(
        "--chunks", required=False, type=int, default=main.READ_CHUNK_INT, show_default=True, help=HELP_CHUNKS
    )
    @wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def sge_extra(f):
    @click.option("--unique", "unique_only", required=False, type=bool, is_flag=True, help=HELP_SGE_UNIQUE)
    @click.option(
        "--chunks", required=False, type=int, default=main.READ_CHUNK_SGE_INT, show_default=True, help=HELP_CHUNKS
    )
    @click.option(
        "-n",
        "--no-alignment",
        required=False,
        default=False,
        type=bool,
        help=HELP_NO_ALIGNMENT,
        show_default=True,
        is_flag=True,
    )
    @wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def common_params(f):
    @click.option("-g", "--guidelib", required=False, default=None, type=_file_exists(), help=HELP_GUIDELIB)
    @click.option("-q", "--queries", required=True, type=_file_exists(), help=HELP_QUERIES)
    @click.option("-s", "--sample", required=False, type=str, help=HELP_SAMPLE)
    @click.option(
        "-o",
        "--output",
        required=True,
        type=click.Path(exists=False, file_okay=True, resolve_path=True),
        help=HELP_OUTPUT,
    )
    @click.option(
        "-w",
        "--workspace",
        required=False,
        type=click.Path(exists=False, dir_okay=True, resolve_path=True),
        help=HELP_WORKSPACE,
    )
    @click.option(
        "--rules",
        required=False,
        type=str,
        default=[],
        multiple=True,
        help=HELP_RULES,
        show_default=True,
    )
    @click.option("--low_count", required=False, type=int, default=None, help=HELP_LOW_COUNT, show_default=True)
    @click.option(
        "--qual_offset",
        required=False,
        type=click.Choice(["33", "64"]),
        default=None,
        help=HELP_QUAL_OFFSET,
        show_default=True,
    )
    @click.option("-c", "--cpus", required=False, type=int, default=1, show_default=True, help=HELP_CPUS)
    @click.option("-m", "--minscore", required=False, type=int, default=15, help=HELP_MINSCORE, show_default=True)
    @click.option("-r", "--reference", required=False, type=_file_exists(), help=HELP_REFERENCE)
    @click.option(
        "-e",
        "--excludeqcf",
        required=False,
        default=False,
        type=bool,
        help=HELP_EXCLUDEQCF,
        show_default=True,
        is_flag=True,
    )
    @click.option(
        "-b",
        "--boundary-mode",
        required=False,
        default="all",
        help=HELP_BOUNDARY,
        show_default=True,
        type=click.Choice(("all", "exact", "TinQ", "QinT"), case_sensitive=False),
    )
    @wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


optgroup_debug = OptionGroup("\nDebug options", help="Options specific to troubleshooting, testing and debugging")


def debug_params(f):
    @optgroup_debug.option(
        "-l",
        "--loglevel",
        required=False,
        default="INFO",
        show_default=True,
        type=click.Choice(LOG_LEVELS, case_sensitive=False),
        help="Set logging verbosity",
    )
    @wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


@click.group()
@click.version_option(pkg_resources.require(__name__.split(".")[0])[0].version)
def cli():  # pragma: no cover
    pass


@cli.command(epilog=HELP_EPILOG)
@common_params
@chunk_default
@click.option(
    "-n",
    "--no-alignment",
    required=False,
    default=False,
    type=bool,
    help=HELP_NO_ALIGNMENT,
    show_default=True,
    is_flag=True,
)
@debug_params
def single_guide(*args, **kwargs):
    """
    Single guide - map read file to library guides and output counts, statistics and alignments files.
    """
    singleguide.run(*args, **kwargs)


@cli.command()
@common_params
@chunk_default
@click.option(
    "-t",
    "--trimseq",
    required=False,
    type=int,
    default=0,
    show_default=True,
    help=HELP_TRIMSEQ,
)
@debug_params
def dual_guide(*args, **kwargs):
    """
    Dual guide - map read file to library guides and output counts, statistics and alignments files.
    """
    dualguide.run(*args, **kwargs)


@cli.command()
@common_params
@sge_extra
@debug_params
def long_read(*args, **kwargs):
    """
    log-read - map read file to library guides and output counts, statistics and alignments files.
    Different to single-guide due to read filtering and other defaults
    """
    pysge.run(*args, **kwargs)


@cli.command()
@click.option("-g", "--guidelib", required=True, type=_file_exists(), help=HELP_GUIDELIB)
@click.option(
    "-f",
    "--fasta",
    required=True,
    type=click.Path(exists=False, file_okay=True, resolve_path=True),
    help=HELP_FASTA,
)
@debug_params
def guides_to_fa(guidelib, fasta, loglevel):
    """
    Convert guide library to fasta file, mainly for debug use via "samtools tview".
    """
    _log_setup(loglevel)
    readwriter.guide_fasta(libparser.load(guidelib), fasta, index=True)
