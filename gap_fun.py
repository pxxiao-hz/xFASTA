#!/usr/bin/env python
'''
time: 2023-04-17
author: pxxiao
version: 1.0
description: 处理基因组的gap信息
'''
import os
from collections import defaultdict
from util import *

import re
from Bio import SeqIO


def func_get_gap_location(file_input):

    # # Parse command-line arguments
    # parser = argparse.ArgumentParser()
    # parser.add_argument("fasta")
    # args = parser.parse_args()
    pfx = os.path.basename(file_input).split(".fa")[0]
    pfx_gff = pfx + ".gaps.gff"
    print("gap.gff: ", pfx_gff)
    pfx_bed = pfx + ".gaps.bed"
    print("gap.bed: ", pfx_bed)
    out_gff = open(pfx_gff, "w")
    out_bed = open(pfx_bed, "w")

    # GFF uses 1-based inclusive coordinates; BED uses 0-based half-open coordinates.
    with open(file_input) as handle:
        gap_index = 0
        for record in SeqIO.parse(handle, "fasta"):
            for match in re.finditer('(?i)N+', str(record.seq)):
                gap_index += 1
                gap_size = match.end() - match.start()
                gff_fields = [
                    record.id, ".", "gap", str(match.start() + 1), str(match.end()), ".", ".", ".",
                    "Name=gap" + str(gap_index) + ";size=" + str(gap_size)
                ]
                out_gff.write("\t".join(gff_fields) + "\n")
                out_bed.write(record.id + "\t" + str(match.start()) + "\t" + str(match.end()) + "\n")
    out_gff.close()
    out_bed.close()
    return
