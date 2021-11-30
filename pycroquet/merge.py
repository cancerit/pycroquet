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
import json
import logging
import sys
from typing import Any
from typing import List

from pycroquet import cli
from pycroquet.classes import Stats

SINGLE_COLS = ("#id", "sgrna_ids", "sgrna_seqs", "gene_pair_id", "unique_guide")
STATS_SUMABLE = (
    "length_excluded_reads",
    "mapped_to_guide_reads",
    "multimap_reads",
    "total_pairs",
    "total_reads",
    "unmapped_reads",
    "vendor_failed_reads",
)


def merge_header_line(cmd: str, version: str, idx: int, chk_item):
    cmd = cmd.lstrip("#")
    version = version.lstrip("#")
    return f"##Count-col-#{idx}: {chk_item}; {version}; {cmd}"


def hash_file(count_file: str, chksum_type) -> str:
    file_hash = None
    if chksum_type == "md5":
        file_hash = hashlib.md5()
    elif chksum_type == "sha256":
        file_hash = hashlib.sha256()
    with open(count_file, "rb") as f:
        while chunk := f.read(8192):
            file_hash.update(chunk)
    return f"{chksum_type}: {file_hash.hexdigest()}"


def merge_single_stats(count_files: List[str]):
    """
    - loads each json into a stats object
    - sum relevant counts into new stats object
    - add original object into ordered list to correlate with countheader
    """
    new_stats = Stats()
    # initialise all numbers
    for k in STATS_SUMABLE:
        setattr(new_stats, k, 0)
    stat_set = []
    for cf in count_files:
        ## need to handle possibility that counts aren't ".gz" compressed
        sf = cf.replace("counts.tsv", "stats.json").replace(".gz", "")
        with open(sf, "rt") as jfp:
            j_data = json.load(jfp)
            this_stats = Stats(**j_data)
            for k, v in vars(this_stats).items():
                if v is None:
                    continue
                if k in STATS_SUMABLE:
                    setattr(new_stats, k, getattr(new_stats, k) + v)
            stat_set.append(this_stats)
    new_stats.merged_from = stat_set
    new_stats.total_guides = stat_set[0].total_guides
    return new_stats


def load_single_count(count_file: str, file_idx: int, chksum_type: str):
    """
    Existing Header will be captured into single line and numbered. Cannot handle adding file to existing merged file.
    Mergins will only capture the id and sample count.
    """
    chk_item = hash_file(count_file, chksum_type)
    try:
        i_fh = gzip.open(count_file, "rt")
        line = i_fh.readline()
    except gzip.BadGzipFile:
        i_fh = open(count_file, "rt")
        line = i_fh.readline()

    cmd = None
    version = None
    sample = None

    data_set = []

    while True:
        if line == "":
            break
        line = line.strip()
        if line.startswith("#"):
            if line.startswith("##Command"):
                if "pycroquet single-guide" not in line:
                    logging.critical(f"Command header does not indicate single-guide: {line}")
                    sys.exit(2)
                cmd = line
            elif line.startswith("##Version"):
                version = line
            elif line.startswith("#id"):
                sample = line.split("\t")[-1]
        else:
            # just capture core fields (id, sgrna_ids, sgrna_seqs, gene_pair_id, unique_guide) and count
            items = line.split("\t")
            data_set.append([*items[0:5], int(items[-1])])

        line = i_fh.readline()
    i_fh.close()

    return (merge_header_line(cmd, version, file_idx, chk_item), sample, data_set)


def merge_single_data(inputs: List[str], chksum_type: str):
    input_idx = 0
    exp_sample = None
    header_lines = []
    final_data = None
    for input in inputs:
        input_idx += 1
        (header_line, sample, input_set) = load_single_count(input, input_idx, chksum_type)
        header_lines.append(header_line)
        if exp_sample is None:
            exp_sample = sample
        if exp_sample != sample:
            logging.critical(f"Input file {input_idx} ({input}) is a different sample to previous files.")
            sys.exit(2)

        if final_data is None:
            final_data = input_set
            continue
        if len(final_data) != len(input_set):
            logging.critical(f"Input file {input_idx} ({input}) has a different number of data rows to previous files.")
            sys.exit(2)
        for idx, row in enumerate(final_data):
            for idx_c, col_name in enumerate(SINGLE_COLS):
                # check the ID matches then extend the row
                if row[idx_c] != input_set[idx][idx_c]:
                    logging.critical(
                        f"Input file {input_idx} ({input}) has a different '{col_name}'' ({row[idx_c]} vs {input_set[idx][idx_c]}) on data row {idx+1} (+ header rows)."
                    )
                    sys.exit(2)
            row.append(input_set[idx][-1])
    return (exp_sample, final_data, header_lines)


def output_merged(output: str, sample: str, merged_data: List[List[Any]], header_lines: List[str], merged_stats: Stats):
    new_cols = [*SINGLE_COLS, sample]
    merged_counts = f"{output}.merged-counts.tsv.gz"
    total_count = 0
    with gzip.open(merged_counts, "wt") as ofh:
        sample_idx = 0
        for hl in header_lines:
            sample_idx += 1
            print(hl, file=ofh)
            new_cols.append(f"#{sample_idx}")
        print("\t".join(new_cols), file=ofh)
        for row in merged_data:
            # @todo - need to assess 15/30/custom low guide counts
            summed = sum(row[5:])
            total_count += summed
            if summed == 0:
                merged_stats.zero_count_guides += 1
            if summed < 15:
                merged_stats.low_count_guides_lt_15 += 1
            if summed < 30:
                merged_stats.low_count_guides_lt_30 += 1
            print("\t".join([str(e) for e in [*row[0:5], summed, *row[5:]]]), file=ofh)
    merged_stats.mean_count_per_guide = round(total_count / merged_stats.total_guides, 2)
    stats_file = f"{output}.merged-stats.json"
    with open(stats_file, "wt") as jfh:
        print(merged_stats.as_json(), file=jfh)


def merge_single(output: str, inputs: List[str], checksum: str, loglevel: str):
    # command/version will be required in new header
    # need to validate
    #   All columns are exact match except counts, order is maintained
    #   header of final column must match as this is "sample"
    cli._log_setup(loglevel)
    if len(inputs) < 2:
        logging.critical("At least 2 count files must be provided")
        sys.exit(2)
    merged_stats = merge_single_stats(inputs)
    (sample, merged_data, header_lines) = merge_single_data(inputs, checksum)
    output_merged(output, sample, merged_data, header_lines, merged_stats)
