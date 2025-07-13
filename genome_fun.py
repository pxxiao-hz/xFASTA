#!/usr/bin/env python
'''
time: 2023-04-07
author: pxxiao
version: 1.0
description: 和 genome 有关的函数
包括：基因组 survey
'''


# 脚本目的：二代 reads 进行 survey
import subprocess
import sys, os, argparse
from util import *
from fasta_fun import *
from scipy.signal import find_peaks


def get_file_list(work_dir):
    filelists = []
    for home, dirs, files in os.walk(work_dir):
        for file in files:
            # 文件名列表，包含完整路径
            filelists.append(os.path.join(home, file))
    return filelists


def find_peak(file):
    # 从文件中读取数据
    with open(file, "r") as file:
        data = [int(line.strip().split()[1]) for line in file]

    # 找到数据中的峰值
    peaks, _ = find_peaks(data)
    # 输出峰值的索引和对应的值
    for peak in peaks:
        print(f"峰值位置: {peak}, 峰值值: {data[peak]}")


def genome_survey(kmersize, length, count, threads, path, genome_poly):
    # parser = argparse.ArgumentParser(description=' genome survey based on Illumina')
    # parser.add_argument('-k', '--kmersize', type=int, default=21, help='kmer size [21]')
    # parser.add_argument('-l', '--length', type=int, default=150, help='reads length [150]')
    # parser.add_argument('-c', '--count', type=int, default=10000000, help='High count value of histogram [10000000]')
    # parser.add_argument('-t', '--threads', type=int, default=16, help='threads [16]')
    # parser.add_argument('-p', '--path', type=str, help='path of fastq reads, decompress format')
    # args = parser.parse_args()
    #
    # kmersize = args.kmersize
    # length = args.length
    # count = args.count
    # threads = args.threads
    # path = args.path

    ### file list
    file_list = get_file_list(path)
    file_str = ' '.join(file_list)

    ### jellyfish
    cmd_jellyfish_count = f'jellyfish count -C -m {kmersize} -s 1000000000 -t {threads} ' \
                          f'{file_str} -o reads.{kmersize}.jf'
    cmd_jellyfish_histo = f'jellyfish histo -h {count} -t {threads} ' \
                          f'-o reads.{kmersize}.histo -v reads.{kmersize}.jf'
    print(cmd_jellyfish_count)
    print(cmd_jellyfish_histo)
    print('Start jellyfish count')
    check_file_in_path(f'reads.{kmersize}.jf', cmd_jellyfish_count)
    print('Start jellyfish histo')
    check_file_in_path(f'reads.{kmersize}.histo', cmd_jellyfish_histo)

    ### genome scope v1
    cmd_genomescope1 = f'Rscript ' \
                      f'/home/pxxiao/tools/GenomeScope_1/genomescope/genomescope.R ' \
                      f'reads.{kmersize}.histo ' \
                      f'{kmersize} ' \
                      f'{length} ' \
                      f'{kmersize}_genomescope1_output'
    print('Start genomescope v1')
    print(cmd_genomescope1)
    check_file_in_path(f'{kmersize}_genomescope1_output/summary.txt', cmd_genomescope1)

    ### genome scope v2
    if genome_poly == 2:
        cmd_genomescope2 = f'Rscript ' \
                          f'/home/pxxiao/tools/GenomeScope_2/genomescope2.0/genomescope.R ' \
                          f'-i reads.{kmersize}.histo ' \
                          f'-o {kmersize}_genomescope2_output ' \
                          f'-k {kmersize}'
    else:
        cmd_genomescope2 = f'Rscript ' \
                           f'/home/pxxiao/tools/GenomeScope_2/genomescope2.0/genomescope.R ' \
                           f'-i reads.{kmersize}.histo ' \
                           f'-o {kmersize}_genomescope2_output ' \
                           f'-k {kmersize} ' \
                           f'-p {genome_poly}'
    print(f'Start genomescope v2, genome poly: {genome_poly}')
    print(cmd_genomescope2)
    check_file_in_path(f'{kmersize}_genomescope2_output/summary.txt', cmd_genomescope2)

    ### find peak
    # find_peak(f'reads.{kmersize}.histo')

    ### smudgeplot
    cmd_smudgeplot = f'px_shell_smudgeplot.sh reads.{kmersize}.histo reads.{kmersize}.jf kmer_pairs_{kmersize} {kmersize}'
    print('Start smudgeplot')
    print(cmd_smudgeplot)
    # check_file_in_path(f'{kmersize}_smudgeplot.png', cmd_smudgeplot)

    ### findGSE
    cmd_findGSE = f'px_shell_findGSE.sh reads.{kmersize}.histo {kmersize}_findGSE_output {kmersize}'
    print('Start findGSE')
    print(cmd_findGSE)
    check_file_in_path(f'{kmersize}_findGSE_output/v1.94.est.reads.{kmersize}.histo.genome.size.estimated.k{kmersize}to{kmersize}.fitted.txt', cmd_findGSE)


def genome_telomere(genome_fasta):
    '''
    使用脚本鉴定端粒；需要借助 外部脚本； 后期可以整合进下FASTA！！！
    :param genome_fasta: genome FASTA
    :return:
    '''

    cmd_telomere_identified = 'python ~/test/scripts/04_T2T/03_evaluate_tools/02_find_telomere.py ' \
                              f'-i {genome_fasta} ' \
                              '-l 150000 ' \
                              '-o temp.telo.info'
    subprocess.run(cmd_telomere_identified, shell=True, close_fds=True)
    return None

