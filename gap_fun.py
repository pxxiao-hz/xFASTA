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

import argparse
import re
from Bio import SeqIO


def func_get_gap_location(input):

    # # Parse command-line arguments
    # parser = argparse.ArgumentParser()
    # parser.add_argument("fasta")
    # args = parser.parse_args()
    pfx = os.path.basename(input).split(".fa")[0]
    pfx_gff = pfx + ".gaps.gff"
    print("gap.gff: ", pfx_gff)
    pfx_bed = pfx + ".gaps.bed"
    print("gap.bed: ", pfx_bed)
    out_gff = open(pfx_gff, "w")
    out_bed = open(pfx_bed, "w")

    # Open FASTA, search for masked regions, print in GFF3 format
    with open(input) as handle:
        i = 0
        for record in SeqIO.parse(handle, "fasta"):
            for match in re.finditer('(?i)N+', str(record.seq)):
                i = i + 1
                # print(record.id, ".", "gap", match.start() + 1, match.end(), ".", ".", ".",
                #       "Name=gap" + str(i) + ";size=" + str(match.end() - match.start()), sep='\t')
                list = [record.id, ".", "gap", str(match.start() + 1), str(match.end()), ".", ".", ".",
                      "Name=gap" + str(i) + ";size=" + str(match.end() - match.start())]
                out_gff.write("\t".join(list) + "\n")
                out_bed.write(record.id + "\t" + str(int(match.start())) + "\t" + str(match.end()+1) + "\n")
    out_gff.close()
    out_bed.close()
    return