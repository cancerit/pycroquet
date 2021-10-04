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
import multiprocessing as mp
import random
from functools import partial
from time import time
from typing import Dict
from typing import List
from typing import Tuple

from pygas.alignercpu import AlignerCpu
from pygas.classes import AlignmentBatch
from pygas.classes import Backtrack

import pycroquet.tools as ctools
from pycroquet import readparser
from pycroquet.classes import Library
from pycroquet.classes import Stats
from pycroquet.countwriter import guide_counts_single

READ_CHUNK_INT = 20000
READ_CHUNK_SGE_INT = 1000


def map_thread(query_seqs: List[str], aligner: AlignerCpu) -> AlignmentBatch:
    return aligner.align_queries(query_seqs, keep_matrix=False)


def mapping_by_chunk(aligner: AlignerCpu, query_seqs, read_chunk, workspace, usable_cpu) -> List[str]:
    # randomizes read order to distribute harder tasks, result is still reproducible
    random.Random().shuffle(query_seqs)
    # library.targets  # the targets for pygas
    # query_seqs  # the queries for pygas
    seq_sets = ctools.chunks(query_seqs, read_chunk * usable_cpu)
    del query_seqs  # make sure that it is removed as it's been duplicated

    # setup for multi-process

    if __name__ == "__main__" and mp.get_start_method(allow_none=True) is None:
        mp.set_start_method("spawn")

    start = time()
    was = start
    pickled_files = []
    for seq_set in seq_sets:
        this_seq_set = list(ctools.chunks(seq_set, read_chunk))
        results = None
        if usable_cpu == 1 or len(this_seq_set) == 1:  # save overhead
            results = [map_thread(this_seq_set[0], aligner)]
        else:
            with mp.Pool(processes=usable_cpu) as pool:
                results = pool.map(partial(map_thread, aligner=aligner), this_seq_set)
        now = time()
        logging.info(f"{usable_cpu} CPUs processed {len(seq_set)} reads in {int(now - was)}s (wall)")
        was = now

        pickled_files.append(ctools.pickle_this(workspace, "pre_matrix_{:05d}".format(len(pickled_files) + 1), results))
        del results
        # logging.warning("exiting loop early for development purposes")
        # break
    return pickled_files


def map_reads(
    aligner: AlignerCpu,
    unique: int,
    read_chunk: int,
    cpus: int,
    workspace: str,
    query_seqs: List[str],
):
    if unique < read_chunk * cpus:
        new_chunk = int(unique / cpus) + 1
        logging.warning(f"--chunks value {read_chunk} rescaled to {new_chunk} to utilise all CPUs")
        read_chunk = new_chunk

    pickles = mapping_by_chunk(aligner, query_seqs, read_chunk, workspace, cpus)
    return pickles


def sg_select_alignment(hits: List[Backtrack], rules: List[str]) -> List[Backtrack]:
    """
    Single guide is very straight forward for selecting the mapping
    - best score that fulfils the rules
    """
    best = []
    best_score = 0
    for bt in hits:
        sm = bt.sm
        if sm.score > best_score:
            best = []
            best_score = sm.score
        elif sm.score < best_score:
            continue
        if bt.pass_rules(rules):
            best.append(bt)
    return best


def process_reads(
    library: Library,
    seqfile: str,
    workspace,
    rules: List[str],
    minscore: int,
    cpus: int,
    read_chunk: int,
    sample=None,
    reference=None,
    exclude_qcfail=None,
    reverse=False,
    exclude_by_len=None,
    boundary_mode=3,
) -> Tuple[Dict[str, int], Dict[str, Tuple[str, List[Backtrack]]], Stats]:
    (unique, stats, query_dict, _) = readparser.parse_reads(
        seqfile,
        sample=sample,
        cpus=cpus,
        reference=reference,
        exclude_qcfail=exclude_qcfail,
        exclude_by_len=exclude_by_len,
    )
    aligner = AlignerCpu(
        targets=library.targets,
        rules=rules,
        score_min=minscore,
        rev_comp=reverse,
        match_type=boundary_mode,
    )

    pickles = map_reads(aligner, unique, read_chunk, cpus, workspace, list(query_dict.keys()))

    # here we are collecting the results into a dict so we can assess them as we pass over the read file again
    aligned_results = {}
    guide_results = {}
    (mapped, multimap, unmapped) = (0, 0, 0)
    for p in pickles:
        logging.info(f"Collating data from {p}")
        batches = ctools.unpickle(p)
        ab: AlignmentBatch
        for ab in batches:
            # unmapped is a simple list
            for s in ab.unmapped:
                aligned_results[s] = ("unmapped", None)
                unmapped += query_dict[s]
            for hits in ab.mapped:
                # function will need to be split out to work via:
                #  library.header.is_single
                best_bt = sg_select_alignment(hits, aligner.rules)
                # for ease of access, common to all hits
                original_seq = hits[0].sm.original_seq
                if len(best_bt) == 0:
                    aligned_results[original_seq] = ("unmapped", None)
                    unmapped += query_dict[original_seq]
                    continue
                elif len(best_bt) > 1:
                    aligned_results[original_seq] = ("multimap", best_bt)
                    multimap += query_dict[original_seq]
                    continue

                # don't need to check for existence on this one
                aligned_results[original_seq] = ("unique", best_bt)
                best_align = best_bt[0]
                # but do here as read with mm/d/i can map to the guide
                if best_align.sm.target not in guide_results:
                    guide_results[best_align.sm.target] = 0
                guide_results[best_align.sm.target] += query_dict[original_seq]
                mapped += query_dict[original_seq]

    logging.info(f"Mapped: {mapped}, Multimap: {multimap} , Unmapped: {unmapped}")
    stats.mapped_to_guide_reads = mapped
    stats.multimap_reads = multimap
    stats.unmapped_reads = unmapped
    stats.total_guides = len(library.guides)
    return (query_dict, guide_results, aligned_results, stats)
