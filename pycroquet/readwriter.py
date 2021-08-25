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
import hashlib
import logging
import os
from array import array
from tempfile import TemporaryDirectory
from typing import Dict
from typing import Iterator
from typing import List
from typing import Tuple

import pysam
from pygas.classes import Backtrack
from pygas.matrix import revcomp

from pycroquet.classes import Guide
from pycroquet.classes import Library
from pycroquet.classes import Seqread
from pycroquet.classes import Stats
from pycroquet.constants import EXT_TO_HTS
from pycroquet.htscomm import hts_reader
from pycroquet.htscomm import hts_sort_n_index
from pycroquet.readparser import is_gzip
from pycroquet.readparser import parse_fq_header


def _hts_iter(
    seq_file,
    mode,
    exclude_qcfail,
    reverse,
    reference,
    cpus=1,
    default_rgid=None,
    trim_len=0,
) -> Iterator[Seqread]:
    sam = hts_reader(seq_file, mode, cpus, reference)
    for read in sam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        if read.is_qcfail and exclude_qcfail:
            continue
        seq = read.get_forward_sequence()
        qual = read.get_forward_qualities()
        if trim_len:
            seq = seq[0:trim_len]
            qual = qual[0:trim_len]
        if reverse:
            seq = revcomp(seq)
            qual = qual[::-1]
        member = (1 if read.is_read1 else 2) if read.is_paired else None
        rgid = read.get_tag("RG") if read.has_tag("RG") else default_rgid
        yield Seqread(
            qname=read.query_name,
            sequence=seq,
            member=member,
            qual=qual,
            qc_fail=read.is_qcfail,
            rgid=rgid,
        )
    sam.close()


def _fq_iter(ifh, exclude_qcfail, reverse, offset_override=None, default_rgid=None, trim_len=0) -> Iterator[Seqread]:
    header = ifh.readline()
    while header:
        seq = ifh.readline().strip()
        _ = ifh.readline()  # throw away separator
        qual_str = ifh.readline().strip()
        if trim_len:
            seq = seq[0:trim_len]
            qual_str = qual_str[0:trim_len]
        (qname, member, qc_fail, phred_offset) = parse_fq_header(header.strip())  # for validation only here
        if qc_fail and exclude_qcfail:
            header = ifh.readline()
            continue
        if offset_override:
            phred_offset = offset_override
        phred_arr = [i - phred_offset for i in map(ord, qual_str)]

        if reverse:
            phred_arr.reverse()
            seq = revcomp(seq)
        yield Seqread(
            qname=qname,
            member=member,
            sequence=seq,
            qual=array("i", phred_arr),
            qc_fail=qc_fail,
            rgid=default_rgid,
        )
        header = ifh.readline()
    ifh.close()


def read_iter(
    seq_file: str,
    reverse: bool = False,
    exclude_qcfail: bool = False,
    reference: str = None,
    offset=None,
    default_rgid=None,
    cpus=1,
    trim_len=0,
) -> Iterator[Seqread]:
    base_iter = None
    ext = os.path.splitext(seq_file)[1]
    if ext in EXT_TO_HTS:
        base_iter = _hts_iter(
            seq_file,
            EXT_TO_HTS[ext],
            exclude_qcfail,
            reverse,
            reference,
            default_rgid=default_rgid,
            cpus=cpus,
            trim_len=trim_len,
        )
    else:
        fq_fo = None
        if is_gzip(seq_file):
            fq_fo = gzip.open(seq_file, "rt")
        else:
            fq_fo = open(seq_file, "rt")
        base_iter = _fq_iter(
            fq_fo,
            exclude_qcfail,
            reverse,
            offset_override=offset,
            default_rgid=default_rgid,
            trim_len=trim_len,
        )
    return base_iter


def to_alignment(seqread: Seqread, reverse: bool, tags: List[Tuple[str, str]], unmapped=False) -> pysam.AlignedSegment:
    tags.append(("RG", seqread.rgid))
    a = pysam.AlignedSegment()
    if reverse:
        a.query_sequence = revcomp(seqread.sequence)
        seqread.qual.reverse()
        a.query_qualities = seqread.qual
        a.is_reverse = True
    else:
        a.query_sequence = seqread.sequence
        a.query_qualities = seqread.qual

    a.query_name = seqread.qname
    a.tags = tags
    a.is_qcfail = seqread.qc_fail
    a.is_unmapped = unmapped

    if seqread.member:
        a.is_paired = True
        if seqread.member == 1:
            a.is_read1 = True
            a.is_read2 = False
        elif seqread.member == 2:
            a.is_read1 = False
            a.is_read2 = True
        else:
            raise ValueError("Only, unpaired, or paired 1/2 reads are anticipated")
    else:
        a.is_paired = False
    return a


def rg_pg(seq_file: str, stats: Stats) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    ext = os.path.splitext(seq_file)[1]
    rg = {"ID": "1", "SM": stats.sample_name}
    prog = " ".join(stats.command.split()[0:2])
    pg = {
        "ID": "pycroquet",
        "PN": prog,
        "CL": stats.command,
        "VN": stats.version,
        "DS": "python Crispr Read to Oligo QUantification Enhancement Tool",
    }

    if ext not in EXT_TO_HTS:
        return ([rg], [pg])

    sam = hts_reader(seq_file, EXT_TO_HTS[ext], 1)
    header = sam.header.as_dict()

    pgs = []
    if "PG" in header:
        pg["PP"] = header["PG"][-1]["ID"]
        pgs = header["PG"]
    pgs.append(pg)

    rgs = []
    if "RG" not in header:
        rgs.append(rg)
    else:
        for r in header["RG"]:
            r["SM"] = stats.sample_name
            rgs.append(r)

    return (rgs, pgs)


def guide_fasta(library: Library, guide_fa: str, index: bool = False):
    sq_data = []
    ref_ids = {}
    with open(guide_fa, "wt") as gfa:
        ref_idx = 0
        for g in library.guides:
            for sg_idx, sgrna_id in enumerate(g.sgrna_ids):
                if sgrna_id in ref_ids:
                    continue
                print(f">{sgrna_id}", file=gfa)
                print(f"{g.sgrna_seqs[sg_idx]}", file=gfa)
                sq_data.append(
                    {
                        "LN": len(g.sgrna_seqs[sg_idx]),
                        "SN": sgrna_id,
                        "M5": hashlib.md5(g.sgrna_seqs[sg_idx].encode("utf-8")).hexdigest(),
                    }
                )
                ref_ids[sgrna_id] = ref_idx
                ref_idx += 1
    pysam.faidx(guide_fa)
    return (sq_data, ref_ids)


def guide_header(workspace: str, library: Library, stats: Stats, seq_file: str) -> Tuple[str, dict, Dict[str, int]]:
    guide_fa = os.path.join(workspace, "guides.fa")
    (sq_data, ref_ids) = guide_fasta(library, guide_fa)

    (rg, pg) = rg_pg(seq_file, stats)
    default_rgid = None
    if len(rg) == 1:
        default_rgid = rg[0]["ID"]
    header = {"HD": {"VN": "1.0"}, "SQ": sq_data, "RG": rg, "PG": pg}
    return (guide_fa, header, ref_ids, default_rgid)


def sm_from_guide_hit(library: Library, hit: Backtrack, strand: str, mapq: int = 60) -> Tuple[Dict[str, str], int]:
    sgrna_ids = library.sgrna_ids_by_seq(hit.sm.target)
    sa_set = {}
    for sgrna_id in sgrna_ids:
        sa_set[sgrna_id] = f"{sgrna_id},{hit.t_pos},{strand},{hit.cigar},{mapq},{hit.nm}"
    return sa_set


def to_mapped_read(
    seqread,
    ref_ids,
    library: Library,
    hits: List[Backtrack],
    reverse_in=None,
    strand_in=None,
    guide_idx=-1,
) -> Tuple[pysam.AlignedSegment, bool]:
    mapq = 60
    if len(hits) > 1:
        mapq = 0
    sa_set = {}
    strand = strand_in
    reverse = reverse_in
    for hit in hits:
        if strand_in is None:
            strand = "-" if hit.sm.reversed else "+"
        this_sa_set = sm_from_guide_hit(library, hit, strand, mapq)
        sa_set = sa_set | this_sa_set
    is_secondary = False
    for hit in hits:
        ref_names = library.sgrna_ids_by_seq(hit.sm.target)
        tags = [("NM", hit.nm), ("MD", hit.md), ("AS", hit.sm.score)]
        if guide_idx >= 0 and is_secondary is False:
            # only apply to the primary mapping
            tags.append(("YG", library.guides[guide_idx].id))
        for ref_name in ref_names:
            if len(sa_set) > 1:
                tags.append(
                    (
                        "SA",
                        ";".join([v for k, v in sa_set.items() if k != ref_name]),
                    )
                )
                tags.append(("NH", len(sa_set)))
            if reverse_in is None:
                reverse = hit.sm.reversed

            a = to_alignment(seqread, reverse, tags)
            # complete mapped info
            a.is_secondary = is_secondary
            a.reference_id = ref_ids[ref_name]
            a.reference_start = hit.t_pos - 1  # zero based
            a.cigarstring = hit.cigar
            a.mapping_quality = mapq
            is_secondary = True
    return (a, len(sa_set) > 1)


def write_alignments(
    align_file,
    header,
    guide_fa,
    cpus,
    seq_file,
    reverse,
    exclude_qcfail,
    reference,
    qual_offset,
    aligned_results,
    library,
    ref_ids,
    strand,
    default_rgid,
    skipped=None,
):
    """
    own function to allow try due to incorrect qual threshold use
    """
    with pysam.AlignmentFile(align_file, "wb", header=header, reference_filename=guide_fa, threads=cpus) as af:
        # get the common-iterator
        iter = read_iter(
            seq_file,
            reverse,
            exclude_qcfail,
            reference,
            offset=qual_offset,
            default_rgid=default_rgid,
        )
        skipped_reads = 0
        for seqread in iter:
            # * Pair the sequence to the mapped records
            if seqread.sequence not in aligned_results:
                if skipped is not None:
                    skipped_reads += 1
                    if skipped_reads > skipped:
                        raise ValueError(
                            "Source read not found in aligned data and number of reads skipped due to length has beed exceeded. Has input been modified during execution?"
                        )
                    a = to_alignment(seqread, False, [], unmapped=True)
                    af.write(a)
                    continue

                raise ValueError("Source read not found in aligned data, has input been modified during execution?")
            (hit_type, hits) = aligned_results[seqread.sequence]

            # add the basics to the alignment

            a = None
            if hit_type in ("unique", "multimap"):
                (a, _) = to_mapped_read(seqread, ref_ids, library, hits, reverse, strand_in=strand)
            elif hit_type == "unmapped":
                a = to_alignment(seqread, reverse, [], unmapped=True)
            af.write(a)


def reads_to_hts(
    library: Library,
    aligned_results: Dict[str, Tuple[str, List[Backtrack]]],
    seq_file: str,
    qual_offset: str,
    workspace: str,
    stats: Stats,
    output: str,
    cpus: int,
    reverse: bool = False,
    exclude_qcfail: bool = False,
    reference: str = None,
):
    if qual_offset:
        qual_offset = int(qual_offset)
    strand = "-" if reverse else "+"
    # * generate the fasta for the guides in workspace
    (guide_fa, header, ref_ids, default_rgid) = guide_header(workspace, library, stats, seq_file)

    # * Create an AlignmentFile with relevant "sequences"
    align_file = os.path.join(workspace, "temp.bam")

    logging.info(f"Writing alignment file: {align_file} (intermediate)")

    sam_threads = cpus if cpus < 4 else 4

    try:
        write_alignments(
            align_file,
            header,
            guide_fa,
            sam_threads,
            seq_file,
            reverse,
            exclude_qcfail,
            reference,
            qual_offset,
            aligned_results,
            library,
            ref_ids,
            strand,
            default_rgid,
            skipped=stats.length_excluded_reads,
        )
    except OverflowError as e:
        seq_extn = os.path.splitext(seq_file)[1].lower
        if str(e) != "unsigned byte integer is less than minimum" and seq_extn in (
            ".sam",
            ".bam",
            ".cram",
        ):
            raise e
        logging.error("FASTQ: phred offset not deduced correctly from readname, attempting to handle internally")
        write_alignments(
            align_file,
            header,
            guide_fa,
            sam_threads,
            seq_file,
            reverse,
            exclude_qcfail,
            reference,
            33,
            aligned_results,
            library,
            ref_ids,
            strand,
            default_rgid,
            skipped=stats.length_excluded_reads,
        )

    hts_sort_n_index(align_file, guide_fa, output, workspace, cpus=sam_threads)
