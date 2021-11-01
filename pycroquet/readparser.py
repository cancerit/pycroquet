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
import logging
import os
import re
from time import time
from typing import Dict
from typing import TextIO
from typing import Tuple

import magic
import pysam
from pygas.matrix import revcomp

from pycroquet.classes import Stats
from pycroquet.constants import EXT_TO_HTS
from pycroquet.htscomm import hts_reader


ILLUMINA_SINGLE_FASTQ_HEADER_PATTERN = re.compile(r"^@([^\s/]+)$")
ILLUMINA_FASTQ_HEADER_PATTERN = re.compile(r"^@(\S+)/([12])$")
CASAVA_FASTQ_HEADER_PATTERN = re.compile(r"^@(\S+)\s([012]):([YN])+:[\d+]+:\S+$")

LOAD_INFO_THRESHOLD = 1000000
OFFSET_ILLUMINA = 64
OFFSET_CASAVA = 33


def parse_fq_header(header: str):
    qc_fail = False
    phred_offset = OFFSET_ILLUMINA
    if (match := ILLUMINA_SINGLE_FASTQ_HEADER_PATTERN.match(header)) is not None:
        # illumina format, unparied
        groups = match.groups()
        name = groups[0]
        pair_member = None
    elif (match := ILLUMINA_FASTQ_HEADER_PATTERN.match(header)) is not None:
        # illumina format, paired or single end wit read identifier
        groups = match.groups()
        name = groups[0]
        pair_member = int(groups[1])
    elif (match := CASAVA_FASTQ_HEADER_PATTERN.match(header)) is not None:
        # casava1.8+ format
        phred_offset = OFFSET_CASAVA
        groups = match.groups()
        name = groups[0]
        pair_member = int(groups[1])
        if groups[2] == "Y":
            qc_fail = True
    else:
        raise ValueError(f"Unsupported FastQ header format: {header}")
    return (name, pair_member, qc_fail, phred_offset)


def is_gzip(seq_file):
    magic_types = magic.from_file(seq_file)
    if "gzip compressed data" in magic_types or "gzip compatible" in magic_types:
        return True
    return False


def parse_reads(
    seq_file: str,
    sample,
    cpus,
    reference=None,
    exclude_qcfail=False,
    reverse=False,
    exclude_by_len=None,
    paired=False,
    trim_len=0,
) -> Tuple[int, Stats, Dict[str, int]]:
    """
    This function is for the initial collation of unique read sequences in the original orientation only (hts will do revcomp).
    The is no chunking of data so this relies on large memory lookups at present.

    Selecting correct underlying parser is via file extension:
    - cram/bam/sam -> htslib processing
    - gz assume gzip fastq
    - anything else assume uncompressed fastq

    Returns:
    Tuple[int unique, int total, List[str seq]]
    """
    start = time()
    ext = os.path.splitext(seq_file)[1]

    response = None

    if ext in EXT_TO_HTS:
        logging.info(f"Sequence input detected as *{ext}")
        response = parse_htsfile(
            seq_file,
            EXT_TO_HTS[ext],
            sample,
            cpus,
            reference=reference,
            exclude_qcfail=exclude_qcfail,
            reverse=reverse,
            exclude_by_len=exclude_by_len,
            paired=paired,
            trim_len=trim_len,
        )
    else:
        fq_fh = None
        if is_gzip(seq_file):
            logging.info("Sequence input detected as gzip (assume fastq)")
            fq_fh = gzip.open(seq_file, "rt")
        else:
            logging.info("uncompressed data (assume fastq)")
            fq_fh = open(seq_file)
        response = parse_fastq(
            fq_fh,
            sample,
            exclude_qcfail=exclude_qcfail,
            reverse=reverse,
            exclude_by_len=exclude_by_len,
            paired=paired,
            trim_len=trim_len,
        )
    if paired:
        logging.info(f"Parsed {response[1].total_pairs} pairs, {len(response[3])} were unique...")
    else:
        logging.info(f"Parsed {response[1].total_reads} reads, {response[0]} were unique...")
    logging.info(f"Read parsing took: {int(time() - start)}s")
    return response


def _hts_sample(sam: pysam.AlignmentFile, sample=None):
    if sample is not None:
        return sample
    header = sam.header.as_dict()
    if "RG" not in header:
        return sample
    for rg in header["RG"]:
        if "SM" not in rg:
            continue
        if sample is None:
            sample = rg["SM"]
        elif sample != rg["SM"]:
            raise ValueError("Multiple different sample names found in header")
    return sample


def parse_htsfile(
    seq_file: str,
    mode: str,
    sample,
    cpus,
    reference=None,
    exclude_qcfail=False,
    reverse=False,
    exclude_by_len=None,
    paired=False,
    trim_len=0,
) -> Tuple[int, Stats, Dict[str, int]]:
    sam = hts_reader(seq_file, mode, cpus, reference)

    stats = Stats(sample_name=_hts_sample(sam, sample=sample))
    if stats.sample_name is None:
        raise ValueError("No sample name found in input file header, please provide via '--sample'")

    reads = {}
    (unique, total, len_ex, total_pairs, unique_pairs) = (0, 0, 0, 0, 0)
    last_read = None
    last_seq = None
    pairs = {} if paired else None
    for read in sam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        if read.is_qcfail:
            if exclude_qcfail:
                continue
            # only count them if we keep them
            stats.vendor_failed_reads += 1

        seq = read.get_forward_sequence()
        if trim_len:
            seq = seq[0:trim_len]
        if exclude_by_len and len(seq) < exclude_by_len:
            len_ex += 1
            continue
        if reverse:
            seq = revcomp(seq)
        if paired:
            if not read.is_paired:
                continue
            if read.is_read1:
                last_read = read.query_name
                last_seq = seq
            else:
                if last_read != read.query_name:
                    raise ValueError("Paired reads require collation before parsing")
                pair_seq = f"{last_seq}|{seq}"
                total_pairs += 1
                if pair_seq in pairs:
                    pairs[pair_seq] += 1
                else:
                    pairs[pair_seq] = 1
                    unique_pairs += 1
        if seq in reads:
            reads[seq] += 1
        else:
            reads[seq] = 1
            unique += 1
        total += 1
        if total % LOAD_INFO_THRESHOLD == 0:  # pragma: no cover
            if paired:
                logging.debug(f"Parsed {total} reads, {unique} were unique...")
            else:
                logging.debug(f"Parsed {total_pairs} pairs, {unique_pairs} were unique...")
    sam.close()
    stats.total_reads = total
    stats.total_pairs = total_pairs
    stats.reversed_reads = reverse
    if exclude_by_len:
        stats.length_excluded_reads = len_ex
    return (unique, stats, reads, pairs)


def parse_fastq(
    ifh: TextIO,
    sample,
    exclude_qcfail=False,
    reverse=False,
    exclude_by_len=None,
    paired=False,
    trim_len=0,
) -> Tuple[int, Stats, Dict[str, int]]:
    """
    Closes received file handle
    Only used for the reads seq minimization process
    """
    stats = Stats(sample_name=sample)
    if stats.sample_name is None:
        raise ValueError("--sample must be provided for fastq inputs")

    if paired:
        raise ValueError("fastq doesn't support paired processing yet")

    reads = {}
    (unique, total, len_ex) = (0, 0, 0)
    pairs = {} if paired else None
    header = ifh.readline()
    while header:
        seq = ifh.readline().strip()
        _ = ifh.readline()  # throw away separator
        _ = ifh.readline()  # throw away qual at this point
        (_, _, qc_fail, _) = parse_fq_header(header.strip())  # for validation only here
        # looks odd, but best way to handle end of file
        header = ifh.readline()
        if qc_fail:
            if exclude_qcfail:
                continue
            # only count them if we keep them
            stats.vendor_failed_reads += 1
        if trim_len:
            seq = seq[0:trim_len]
        if exclude_by_len and len(seq) < exclude_by_len:
            len_ex += 1
            continue
        if reverse:
            seq = revcomp(seq)
        if seq in reads:
            reads[seq] += 1
        else:
            reads[seq] = 1
            unique += 1
        total += 1
        if total % LOAD_INFO_THRESHOLD == 0:  # pragma: no cover
            logging.debug(f"Parsed {total} reads, {unique} were unique...")
    ifh.close()
    stats.total_reads = total
    stats.reversed_reads = reverse
    if exclude_by_len:
        stats.length_excluded_reads = len_ex
    return (unique, stats, reads, pairs)


def collate(seq_file, workspace, cpus):
    """
    collates reads into pairs IF input is a hts file
    """
    start = time()
    hts_cpus = cpus
    if cpus > 4:
        hts_cpus = 4
    ext = os.path.splitext(seq_file)[1]
    if ext not in EXT_TO_HTS:
        return seq_file
    tmp_hts = os.path.join(workspace, "collated.reads.bam")
    logging.info(f"Collating pairs from {seq_file} to {tmp_hts}")
    pysam.collate("-f", "-l", "1", "--no-PG", "-@", str(hts_cpus), "-o", tmp_hts, seq_file)
    logging.info(f"Collation of reads took: {int(time() - start)}s")
    return tmp_hts
