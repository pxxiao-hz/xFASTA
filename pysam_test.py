#!/usr/bin/env python
'''
time: 2023-07-01
author: zcyu && pxxiao
version: 1.0
description: pysam
'''
import subprocess
import pysam
import sys
import os
import argparse


def x_get_base_info(cigar):
    '''
    根据 CIGAR 值，获得match、insertion、deletion、left clipped、right clipped、edit distance的信息
    :param cigar:
    :return:
    '''
    match_base = 0
    insertion_base = 0
    deletion_base = 0
    soft_clipped_base = 0
    hard_clipped_base = 0
    # edit_distance = 0
    for op, length in cigar:
        # print(cigar)
        if op == 0:
            match_base += length
        elif op == 1:
            # print(f'insertion: {length}')
            if length == 9016:
                print(cigar)
            insertion_base += length
        elif op == 2:
            deletion_base += length
        elif op == 4: # soft clip
            soft_clipped_base += length
        elif op == 5:
            hard_clipped_base += length
    if cigar[0][0] == 5:
        left_hard_clipped_base = cigar[0][1]
    else:
        left_hard_clipped_base = 0
    if cigar[-1][0] == 5:
        right_hard_clipped_base = cigar[-1][1]
    else:
        right_hard_clipped_base = 0
    # print(f'match_base: {match_base}')
    # print(f'insertion_base: {insertion_base}')
    # print(f'deletion_base: {deletion_base}')
    # print(f'soft_clipped_base: {soft_clipped_base}')
    # print(f'hard_clipped_base: {hard_clipped_base}')
    # print(f'left_hard_clipped_base: {left_hard_clipped_base}')
    # print(f'right_hard_clipped_base: {right_hard_clipped_base}')
    return match_base, insertion_base, deletion_base, soft_clipped_base, hard_clipped_base, \
           left_hard_clipped_base, right_hard_clipped_base


def x_get_edit_distance(record):
    try:
        edit_distance = record.get_tag('NM')
        return edit_distance
    except KeyError:
        return None


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=' 比对bam文件提取比对信息 ')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help="bam file")
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='output file')
    parser.add_argument('-m', '--mapping_quality', type=int, required=False, default=20,
                        help='Minium mapping quality for filtering [20]')
    parser.add_argument('-v', '--verbose', type=str, help='打印详细信息')
    # 解析命令行参数
    args = parser.parse_args()
    # 处理命令行参数
    file_bam = args.input
    output_file = args.output
    min_mapping_quality = args.mapping_quality
    verbose = args.verbose
    # 是否打印相关的文件信息：true or false
    if verbose == "True" or verbose == "T" or verbose == "true" or verbose == "t" or verbose == "1":
        print(f'输入文件路径：{file_bam}')
        print(f'输出文件路径：{output_file}')
        print('打印详细信息')

    ### 执行其他的命令
    out = open(output_file, 'w')
    out.write('#Query_name\tQuery_hard_length\tQuery_aligned_start\tQuery_aligned_end\tCoverage\tIdentity\t'
              'Ref\tRef_start\tRef_end\tFLAG\tMapping_quality\tCIGAR\n')
    # check bam index
    file_bam_idx = os.path.basename(file_bam) + '.bai'
    if not os.path.exists(file_bam_idx):
        cmd = f'samtools index {file_bam}'
        subprocess.run(cmd, shell=True, close_fds=True)
    # process bam file
    bf = pysam.AlignmentFile(file_bam, 'rb')
    for r in bf.fetch():
        # flag
        flag = r.flag
        # map quality
        mapping_quality = r.mapping_quality
        if mapping_quality < min_mapping_quality:
            continue
        # edit distance
        edit_distance = x_get_edit_distance(r)
        # match, insertion, deletion, soft clipped, hard clipped, left_clipped, right_clipped
        match_bases, insertion_bases, deletion_bases, soft_clipped, hard_clipped,\
            left_clipped_bases, right_clipped_bases = x_get_base_info(r.cigar)
        # query length (not including hard clipped)
        query_length = r.query_length   # the length of the query/read; includes soft-clipped bases; == len(query_sequence)
        # query length (including hard clipped)
        query_length_hard = r.infer_read_length()
        # read 比对的起点
        query_start = r.query_alignment_start # (0-based, first base in query_seq that is not soft-clipped)
        # read 比对的终点
        query_end = r.query_alignment_end
        # read 比对的长度 length of the aligned query sequence
        query_alignment_length = r.query_alignment_length   # == query_alignment_end - query_alignment_start
        # ref length: aligned length of the read on the reference genome
        ref_length = r.reference_length # == reference_end - reference_start

        ## coverage
        coverage = query_alignment_length / query_length_hard
        ## identity
        identity = (query_alignment_length - edit_distance) / query_alignment_length

        if flag == 0 or flag == 2048:
            query_start = r.query_alignment_start + left_clipped_bases  # include Hard clipped
            query_end = r.query_alignment_end + left_clipped_bases
        elif flag == 16 or flag == 2064:
            query_start = r.query_alignment_start + right_clipped_bases
            query_end = r.query_alignment_end + right_clipped_bases

        # print(len(r.query_sequence))
        # print(f'{r.query_name}\t{query_length_hard}\t{query_start}\t{query_end}\t'
        #       f'{coverage}\t{identity}\t{r.reference_name}\t{r.reference_start}\t{r.reference_end}\t'
        #       f'{r.flag}\t{r.mapping_quality}')

    print(r.cigarstring)


if __name__ == '__main__':
    main()