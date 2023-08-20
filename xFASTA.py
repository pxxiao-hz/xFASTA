#!/usr/bin/env python
'''
time: 2023-04-05
author: pxxiao
version: 1.0
description: 用于日常的 FASTA 文件处理，其他功能待定
'''

import argparse
import os

from fasta_fun import *
from gff_fun import *
from gap_fun import *
from util import *
from reads_fun import *
import textwrap

def main():# 定义命令行参数和选项
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='pxx 的分析小工具')
    subparsers = parser.add_subparsers(dest='function', metavar='')

    # 定义 func-one 命令的子命令和选项: fasta
    func_one_parser = subparsers.add_parser('fasta', formatter_class=argparse.RawTextHelpFormatter,
                                            help='处理 FASTA 文件')
    # func_one_parser.add_argument('-s', '--split', dest='split', required=True, choices=['length', 'number'], metavar='',
    #                              help='拆分 FASTA 文件，1按照bp拆分，2按照number拆分')
    func_one_parser.add_argument('-i', '--input', type=str, required=True,
                                 help='file input: .fasta or .fa')
    func_one_parser.add_argument('-m', '--model', dest='model',
                                 choices=['length', 'number', 'extractk', 'extractl', 'extractc', 'falen',
                                          'IDsimp', 'compare', 'toPhy', 'count'], metavar='',
                                 help='Process on .fasta or .fa \n'
                                      'length: Split File Based on Length \n'
                                      'number: Split File Based on Number \n'
                                      'extractk: Sequence Extraction Based on Keywords,\n'
                                      '\t -k required \n'
                                      'extractl: Sequence Extraction Based on seq length,\n'
                                      '\t -l required \n'
                                      'extractc: Sequence Extraction Based on coords, \n'
                                      '\t -k required, e.g.: chr1-10-20; chr1-10-, means 10 - end of chr1 \n'
                                      '\t Note: 1-based in this script \n'
                                      'falen: Length of each sequence \n'
                                      'IDsimp: Simplified ID Name,\n'
                                      '\t -s required \n'
                                      'compare: Compare the contents of two FASTA files for equality \n'
                                      'toPhy: Convert FASTA format to phy format \n'
                                      'count: Count the number of strings')

    func_one_parser.add_argument('-l', '--length', type=int, default=20000000, required=False,
                                 help='拆分成多少bp的小文件 or 根据序列的长度从基因组文件提取序列 [20000000]')
    func_one_parser.add_argument('-n', '--number', type=int, default=50, required=False,
                                 help='拆分成多个文件，每个文件包含number个序列 [50]')
    func_one_parser.add_argument('-k', '--keywords', type=str, required=False,
                                 help='Keywords in the Sequence to be Extracted or Count.\n'
                                      'For Extracted, Chr1 or Chr1,Chr2,Chr3')
    func_one_parser.add_argument('-g', '--greedy', type=str, required=False,
                                 help='Greedy match or not, True or False')
    func_one_parser.add_argument('-s', '--separator', type=str, required=False,
                                 help='Separator of sequence ID; "None" == "\TAB"')
    func_one_parser.add_argument('-x', '--newFASTA', type=str, required=False,
                                 help='File 2 of two FASTA files for compare')
    func_one_parser.add_argument('-p', '--toPhy', type=str, required=False,
                                 help='FASTA file to Phy')


    # 定义 func-two 命令的子命令和选项: stat
    func_stat_parser = subparsers.add_parser('stat', formatter_class=argparse.RawTextHelpFormatter,
                                             help='统计基因组的相关信息')
    func_stat_parser.add_argument('-i', '--input', type=str, required=True,
                                  help='输入文件，待处理的 FASTA 文件')
    func_stat_parser.add_argument('-n', '--n50', type=str, required=False,
                                  help='基因组组装结果的统计')

    # 定义 fun-three 命令的子命令和选项: gff
    func_gff_parser = subparsers.add_parser('gff', formatter_class=argparse.RawTextHelpFormatter,
                                            help='Perform operations on GFF file')
    func_gff_parser.add_argument('-i', '--input', type=str, required=True,
                                  help='file input: GFF file')
    func_gff_parser.add_argument('-m', '--model', type=str, required=True,
                                 choices=['sort', 'extract', 'GetLongTrans', "NewGff"], metavar='',
                                 help= 'GFF sub command, including sort, extract, GetLongTrans\n'
                                       'sort: based on chr and start \n'
                                       'extract: extract exon, intron information, output format: bed\n'
                                       'GetLongTrans: only retain longest information, output format: gff\n'
                                       'NewGff: if gff not have "gene" mark, please use NewGff')
    func_gff_parser.add_argument('-s', '--split', type=str, required=False, default="ID=",
                                 help = 'the 9th column of gff file, maybe startswith "ID=", or "Parent=')

    # 定义 fun-four 命令的子命令和选项：gap
    func_stat_parser = subparsers.add_parser('gap', formatter_class=argparse.RawTextHelpFormatter,
                                             help = 'Genome gap statistics')
    func_stat_parser.add_argument('-i', '--input', type=str, required=True,
                                  help='finle input: Genome fasta file')
    func_stat_parser.add_argument('-m', '--model', type=str, required=True,
                                  choices=['coord'], metavar='',
                                  help = 'gap sub coomand, including coord\n'
                                         'coord: Obtaining the location information of gaps, '
                                         'NNNN... or nnnn... was defined as gap')

    # 定义 fun-five 命令的子命令和选项：reads
    func_stat_parser = subparsers.add_parser('read', formatter_class=argparse.RawTextHelpFormatter,
                                             help = 'Reads information')
    func_stat_parser.add_argument('-i', '--input', type=str, required=True,
                                  help = 'The input FASTQ file, can be in .fastq or .fastq.gz format')
    func_stat_parser.add_argument('-m', '--model', type=str, required=True,
                                  choices=['phreads', 'deduplication', 'extract'], metavar='',
                                  help = 'read sub command, including phreads\n'
                                         'phreads: Parse the quality value of reads\n'
                                         'deduplication: remove reads based on reads ID, if reads ID not uniq\n'
                                         'extract: Extract reads sequences')
    func_stat_parser.add_argument('-k', '--readid', type=str, required=False, metavar='',
                                  help = 'Extract reads sequence based the reads ID')


    # 解析命令行参数
    args = parser.parse_args()

    # 根据所选的命令行选项调用相应的功能和子功能
    if args.function == 'fasta':
        model = args.model
        inputfile = args.input
        length = args.length
        number = args.number
        keywords = args.keywords
        separator = args.separator
        newFASTA = args.newFASTA
        greedy = args.greedy
        if model == 'length':
            split_fasta_based_bp(inputfile, length)
        elif model == 'number':
            split_fasta_based_number(inputfile, number)
        elif model == 'extractk':
            extract_based_keywords(inputfile, keywords, greedy)
        elif model == 'extractl':
            extract_based_length(inputfile, length)
        elif model == 'extractc':
            extract_based_coords(inputfile, keywords)
        elif model == "falen":
            func_get_fasta_length(inputfile)
        elif model == "IDsimp":
            func_fasta_IDsimplified(inputfile, separator)
        elif model == "compare":
            func_fasta_compare(inputfile, newFASTA)
        elif model == "toPhy":
            func_fasta_toPhy(inputfile)
        elif model == "count":
            func_fasta_count(inputfile, keywords)
        else:
            print("Invalid model choice. \nPlease set \"model\" parameter!")

    elif args.function == "stat":
        file_input = args.input
        fun_1 = args.n50
        func_stat(file_input, fun_1)

    elif args.function == "gff":
        input_gff = args.input
        model = args.model
        splits = args.split
        # 如果不是gff或则gff3后缀，直接报错
        print(input_gff)
        if not args.input.lower().endswith(('.gff', '.gff3')):
            raise ValueError('Error: Input file must have a .gff or .gff3 extension.')
        if model == "sort":
            None
        elif model == "extract":
            func_get_ExonIntronInfo(input_gff, splits)
        elif model == "GetLongTrans":
            # get long transcripts
            func_get_LongTranscript(input_gff)
        elif model == "NewGff":
            # 造一个新的GFF文件
            func_gff_NewGFF(input_gff)

    elif args.function == "gap":
        input = args.input
        model = args.model
        print("input genome file: ", input)
        print("model: ", model)
        if not args.input.lower().endswith(('.fasta', '.fa')):
            raise ValueError('Error: Input file must have a .fasta or .fa extension.')
        if model == "coord":
            func_get_gap_location(input)

    elif args.function == "read":
        input = args.input
        model = args.model
        readsID = args.readid
        print("Input FASTQ reads: ", input)
        print("model: ", model)
        if model == "phreads":
            print("Parse the quality value of reads.")
            func_get_reads_phreads(input)
        elif model == "deduplication":
            print("Deduplication based on reads ID, if reads ID not uniq")
            func_reads_deduplication_based_ID(input)
        elif model == 'extract':
            print('Extract reads based on read ID')
            func_reads_extract_based_ID(input, readsID)


if __name__ =="__main__":
    main()