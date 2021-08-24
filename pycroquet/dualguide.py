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
import json
import logging
import os
import shutil
from time import time
from typing import Dict
from typing import Final
from typing import List
from typing import Tuple

import pysam
from pygas.alignercpu import AlignerCpu
from pygas.classes import Backtrack

import pycroquet.tools as ctools
from pycroquet import cli
from pycroquet import libparser
from pycroquet import readparser
from pycroquet.classes import Classification
from pycroquet.classes import Library
from pycroquet.classes import Stats
from pycroquet.htscomm import hts_sort_n_index
from pycroquet.main import map_reads
from pycroquet.main import sg_select_alignment
from pycroquet.readwriter import guide_header
from pycroquet.readwriter import read_iter
from pycroquet.readwriter import to_alignment
from pycroquet.readwriter import to_mapped_read

CLASSIFICATION: Final = Classification()


def best_multimatch_set(library: Library, bt_set_l: List[Backtrack], bt_set_r: List[Backtrack]):
    """
    try to identify a single guide pair when one or both ends are multi-map
    """
    guides = {}  # to hold matches
    last_key = None  # makes getting the single entry far lower impact
    for bt_l in bt_set_l:
        guide_idxs_l = set(library.target_to_guides[bt_l.sm.target])
        for bt_r in bt_set_r:
            guide_idxs_r = set(library.target_to_guides[bt_r.sm.target])
            guide_intersect = guide_idxs_l.intersection(guide_idxs_r)
            if len(guide_intersect) == 1:
                gidx = guide_intersect.pop()
                if bt_l.sm.reversed is False and bt_r.sm.reversed is True:
                    guides[gidx] = (gidx, bt_l, bt_r)
                    last_key = gidx
                    # on multi we don't try to deal with everything
    if len(guides) == 1:
        return guides[last_key]
    return None


def intersect_guide_targets(library, tseq_l, tseq_r):
    guide_idxs_l = set(library.target_to_guides[tseq_l])
    guide_idxs_r = set(library.target_to_guides[tseq_r])
    # gives the common guides
    return guide_idxs_l.intersection(guide_idxs_r)


def classify_read_pair(
    map_l: Tuple[str, List[Backtrack]],
    map_r: Tuple[str, List[Backtrack]],
    library: Library,
) -> Tuple[str, int, Backtrack, Backtrack, str, str]:  #  noqa R701
    """
    Returns:
        classification string (see classes.Classification)
        R1 backtrack
        R2 backtrack
        R1 guide index
        R2 guide index
    """
    (mtype_l, mdata_l) = map_l
    (mtype_r, mdata_r) = map_r
    if mtype_l == "unique" and mtype_l == mtype_r:
        (bt_l, bt_r) = (mdata_l[0], mdata_r[0])
        (tseq_l, tseq_r) = (bt_l.sm.target, bt_r.sm.target)
        gidx = library.guide_by_sgrna_set(tseq_l, tseq_r)
        if gidx is not None and bt_l.sm.reversed is False and bt_r.sm.reversed is True:
            return (CLASSIFICATION.match, gidx, bt_l, bt_r, mtype_l, mtype_r)

        # gives the common guides
        guide_intersect = intersect_guide_targets(library, tseq_l, tseq_r)
        # but now need to know if any give the appropriate orientations
        g_size = len(guide_intersect)
        if g_size == 0:
            return (CLASSIFICATION.swap, -1, bt_l, bt_r, mtype_l, mtype_r)
        if g_size == 1:
            gidx = guide_intersect.pop()
            if bt_l.sm.reversed is False and bt_r.sm.reversed is True:
                return (CLASSIFICATION.match, gidx, bt_l, bt_r, mtype_l, mtype_r)
            else:
                return (CLASSIFICATION.aberrant_match, -1, bt_l, bt_l, mtype_l, mtype_r)
        # g_size > 1
        if bt_l.sm.reversed is False and bt_r.sm.reversed is True:
            for gidx in guide_intersect:
                guide = library.guides[gidx]
                if guide.sgrna_seqs[0] == tseq_l and guide.sgrna_seqs[1] == tseq_r:
                    return (
                        CLASSIFICATION.match,
                        gidx,
                        bt_l,
                        bt_r,
                        mtype_l,
                        mtype_r,
                    )
        return (CLASSIFICATION.aberrant_match, -1, bt_l, bt_l, mtype_l, mtype_r)

    if mtype_l == "unique":
        if mtype_r == "unmapped":
            bt_l = mdata_l[0]
            if bt_l.sm.reversed:
                return (CLASSIFICATION.r_open_3p, -1, bt_l, None, mtype_l, mtype_r)
            else:
                return (CLASSIFICATION.f_open_3p, -1, bt_l, None, mtype_l, mtype_r)
        # then multimap
        best_pair = best_multimatch_set(library, mdata_l, mdata_r)
        if best_pair:
            (
                gidx,
                bt_l,
                bt_r,
            ) = best_pair  # this will impact how we write out alignments
            return (CLASSIFICATION.match, gidx, bt_l, bt_r, mtype_l, mtype_r)
        bt_l = mdata_l[0]
        if bt_l.sm.reversed:
            return (CLASSIFICATION.r_multi_3p, -1, bt_l, None, mtype_l, mtype_r)
        else:
            return (CLASSIFICATION.f_multi_3p, -1, bt_l, None, mtype_l, mtype_r)
    if mtype_r == "unique":
        if mtype_l == "unmapped":
            bt_r = mdata_r[0]
            if bt_r.sm.reversed:
                return (CLASSIFICATION.r_open_5p, -1, None, bt_r, mtype_l, mtype_r)
            else:
                return (CLASSIFICATION.f_open_5p, -1, None, bt_r, mtype_l, mtype_r)
        # then multi
        best_pair = best_multimatch_set(library, mdata_l, mdata_r)
        if best_pair:
            (
                gidx,
                bt_l,
                bt_r,
            ) = best_pair  # this will impact how we write out alignments
            return (CLASSIFICATION.match, gidx, bt_l, bt_r, mtype_l, mtype_r)
        bt_r = mdata_r[0]
        if bt_r.sm.reversed:
            return (CLASSIFICATION.r_multi_5p, -1, None, bt_r, mtype_l, mtype_r)
        else:
            return (CLASSIFICATION.f_multi_5p, -1, None, bt_r, mtype_l, mtype_r)
    if mtype_l == "multimap" and mtype_l == mtype_r:
        best_pair = best_multimatch_set(library, mdata_l, mdata_r)
        if best_pair:
            (
                gidx,
                bt_l,
                bt_r,
            ) = best_pair  # this will impact how we write out alignments
            return (CLASSIFICATION.match, gidx, bt_l, bt_r, mtype_l, mtype_r)
        return (CLASSIFICATION.no_match, -1, None, None, mtype_l, mtype_r)

    # all that remains are the combination of unmapped/multimap on both ends
    return (CLASSIFICATION.no_match, -1, None, None, mtype_l, mtype_r)


def order_hits(hits: List[Backtrack], first_bt: Backtrack):
    if hits is None or first_bt is None:
        return []
    if len(hits) == 1:
        return hits
    reordered = [first_bt]
    for bt in hits:
        if bt.sm.target != first_bt.sm.target:
            reordered.append(bt)
    return reordered


def _init_class_counts():
    counts = {}
    for v in vars(CLASSIFICATION):
        counts[v] = 0
    return counts


def read_pairs_to_guides(
    workspace: str,
    aligned_results: Dict[str, Tuple[str, List[Backtrack]]],
    library: Library,
    seq_file,
    guide_fa,
    header,
    ref_ids,
    default_rgid,
    stats: Stats,
    cpus=1,
    trim_len=0,
):
    start = time()
    class_cache = {}
    # initialise counts
    counts = _init_class_counts()
    align_file = os.path.join(workspace, "tmp.cram")

    with pysam.AlignmentFile(align_file, "wb", header=header, reference_filename=guide_fa, threads=cpus) as af:
        iter = read_iter(seq_file, default_rgid=default_rgid, cpus=cpus, trim_len=trim_len)
        for seqread_l in iter:
            seqread_r = next(iter, None)
            if seqread_r is None:
                raise ValueError("Collated BAM exhausted between records")

            (read_l, read_r) = (seqread_l.sequence, seqread_r.sequence)
            pair_lookup = f"{read_l}|{read_r}"

            a_l, a_r = None, None

            (class_type, guide_idx, hits_l, hits_r, orig_l, orig_r) = (
                None,
                -1,
                [],
                [],
                None,
                None,
            )

            if pair_lookup not in class_cache:
                (map_l, map_r) = (aligned_results[read_l], aligned_results[read_r])
                # try to order from most likely to least
                (
                    class_type,
                    guide_idx,
                    bt_l,
                    bt_r,
                    orig_l,
                    orig_r,
                ) = classify_read_pair(map_l, map_r, library)

                class_cache[pair_lookup] = (
                    class_type,
                    guide_idx,
                    order_hits(map_l[1], bt_l),
                    order_hits(map_r[1], bt_r),
                    orig_l,
                    orig_r,
                )

            (class_type, guide_idx, hits_l, hits_r, orig_l, orig_r) = class_cache[pair_lookup]
            if guide_idx != -1:
                library.guides[guide_idx].count += 1

            if hits_l:
                (a_l, multi) = to_mapped_read(seqread_l, ref_ids, library, hits_l, guide_idx=guide_idx)
                if multi:
                    stats.multimap_reads += 1
                stats.mapped_to_guide_reads += 1
            else:
                a_l = to_alignment(seqread_l, False, [], unmapped=True)
                if class_type == CLASSIFICATION.r_multi_5p:
                    stats.multimap_reads += 1
                else:
                    if orig_l == "multimap":
                        stats.multimap_reads += 1
                    else:
                        stats.unmapped_reads += 1
            if hits_r:
                (a_r, multi) = to_mapped_read(seqread_r, ref_ids, library, hits_r, guide_idx=guide_idx)
                if multi:
                    stats.multimap_reads += 1
                stats.mapped_to_guide_reads += 1
            else:
                a_r = to_alignment(seqread_r, False, [], unmapped=True)
                if class_type == CLASSIFICATION.r_multi_3p:
                    stats.multimap_reads += 1
                else:
                    if orig_r == "multimap":
                        stats.multimap_reads += 1
                    else:
                        stats.unmapped_reads += 1

            af.write(a_l)
            af.write(a_r)

            counts[class_type] += 1
    logging.info(f"Alignment grouping took: {int(time() - start)}s")
    return (counts, align_file)


def pickles_to_mapset(pickles: List[str], reads: Dict[str, int], aligner: AlignerCpu):
    (aligned_results, multi_map, unique_map, unmap) = ({}, 0, 0, 0)
    for p in pickles:
        logging.info(f"Collating data from {p}")
        batches = ctools.unpickle(p)
        ab: AlignmentBatch
        for ab in batches:
            # unmapped is a simple list
            for s in ab.unmapped:
                aligned_results[s] = ("unmapped", None)
                unmap += reads[s]
            for hits in ab.mapped:
                # function will need to be split out to work via:
                #  library.header.is_single
                best_bt = sg_select_alignment(hits, aligner.rules)
                # for ease of access, common to all hits
                original_seq = hits[0].sm.original_seq
                if len(best_bt) == 0:
                    aligned_results[original_seq] = ("unmapped", None)
                    unmap += reads[original_seq]
                    continue
                elif len(best_bt) > 1:
                    aligned_results[original_seq] = ("multimap", best_bt)
                    multi_map += reads[original_seq]
                    continue
                # don't need to check for existence on this one
                aligned_results[original_seq] = ("unique", best_bt)
                unique_map += reads[original_seq]
    return (aligned_results, multi_map, unique_map, unmap)


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
    boundary_mode,
    trimseq,
    chunks,
    loglevel,
):
    (usable_cpu, work_tmp, workspace, boundary_mode) = cli.common_setup(
        loglevel, cpus, workspace, output, boundary_mode
    )
    library = libparser.load(guidelib)
    # returning a list of all guides in a single list
    # this allows us to map reads to both orientations at the same time (set reverse_comp)

    ###
    # input data has to be samtools-collated if hts
    seq_file = readparser.collate(queries, workspace, usable_cpu)
    # read loading needs to return both the unique list of read sequences, and the unique pairs (r1|r2)
    (unique, stats, reads, read_pairs) = readparser.parse_reads(
        seq_file,
        sample=sample,
        cpus=usable_cpu,
        reference=reference,
        exclude_qcfail=excludeqcf,
        paired=True,
        trim_len=trimseq,
    )

    # map the uniq list or individual reads and then use the r1|r2 info to bring the events back together
    aligner = AlignerCpu(
        targets=library.targets,
        rules=rules,
        score_min=minscore,
        rev_comp=True,  # as some reads can be reversed in DG
        match_type=boundary_mode,
    )
    pickles = map_reads(aligner, unique, chunks, cpus, workspace, list(reads.keys()))
    """
    Need to convert alignment batches into a dict by sequence, containing the possible mappings
    """
    (aligned_results, multi_map, unique_map, unmap) = pickles_to_mapset(pickles, reads, aligner)
    logging.info(f"Unique: {unique_map}, Multimap: {multi_map}, Unmapped, {unmap}")
    # * generate the fasta for the guides in workspace
    (guide_fa, header, ref_ids, default_rgid) = guide_header(workspace, library, stats, seq_file)

    (raw_counts, unsorted) = read_pairs_to_guides(
        workspace,
        aligned_results,
        library,
        seq_file,
        guide_fa,
        header,
        ref_ids,
        default_rgid,
        stats,
        cpus=usable_cpu,
        trim_len=trimseq,
    )
    count_output = f"{output}.counts.tsv"

    # initialise
    if low_count is True:
        stats.low_count_guides_user = {"lt": low_count, "count": 0}
    stats.pair_classifications = raw_counts

    logging.info(f"Writing counts file: {count_output}")
    with open(count_output, "wt") as cout:
        print("##Command: " + stats.command, file=cout)
        print("##Version: " + stats.version, file=cout)
        print("#" + "\t".join(["id", stats.sample_name]), file=cout)
        for g in library.guides:
            print(f"{g.id}\t{g.count}", file=cout)
            if g.count == 0:
                stats.zero_count_guides += 1
            if g.count < 15:
                stats.low_count_guides_lt_15 += 1
            if g.count < 30:
                stats.low_count_guides_lt_30 += 1
            if low_count is True and g.count < low_count:
                stats.low_count_guides_user["count"] += 1
            stats.total_guides += 1

    stats_output = f"{output}.stats.json"
    logging.info(f"Writing statistics file: {stats_output}")
    with open(stats_output, "wt") as jout:
        print(json.dumps(stats.__dict__, sort_keys=True, indent=2), file=jout)

    hts_sort_n_index(unsorted, guide_fa, output, workspace, cpus=usable_cpu)

    if not work_tmp:
        shutil.rmtree(workspace)
