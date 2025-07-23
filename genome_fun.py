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


# def genome_split_chr_by_gap(genome_fasta, gap_size):
#     '''
#     在 Gap 处切割基因组，获得contig 水平的基因组，Gap 大小默认 100 bp
#     :param genome_fasta: genome FASTA
#     :return:
#     '''
#
#     def split_sequence_by_gap(sequence, gap="N" * 100):
#         """ 根据给定的gap（默认100个N）将序列切割成contig """
#         contigs = sequence.split(gap)
#         return [contig for contig in contigs if len(contig) > 0]  # 去除空contig
#
#     def write_fasta(output_file, contigs, seq_id):
#         """ 将分割后的序列写入新的FASTA文件 """
#         for i, contig in enumerate(contigs, 1):
#             output_file.write(f">{seq_id}_contig{i}\n")
#             for j in range(0, len(contig), 60):
#                 output_file.write(contig[j:j + 60] + "\n")  # 每行写60个碱基
#
#     def process_fasta(input_fasta, output_fasta, gap_size=100):
#         """ 处理FASTA文件，将序列按gap切割并输出新的FASTA文件 """
#         sequences = fasta_read(input_fasta)
#         gap = "N" * gap_size
#         with open(output_fasta, "w") as output_file:
#             for seq_id, sequence in sequences.items():
#                 contigs = split_sequence_by_gap(sequence, gap)
#                 write_fasta(output_file, contigs, seq_id)
#     # if ".fasta" in os.path.basename(genome_fasta):
#     #     output_fasta = genome_fasta.replace(".fasta", ".contig.fa")
#     # else:
#     #     output_fasta = genome_fasta.replace(".fa", ".contig.fa")
#     # output_fasta = genome_fasta.replace(".fasta", ".contig.fa").replace(".fa", ".contig.fa")
#     base, ext = os.path.splitext(genome_fasta)
#     output_fasta = base + ".contig.fa"
#     process_fasta(genome_fasta, output_fasta, gap_size)


def genome_split_chr_by_gap(genome_fasta, gap_size=None):
    '''
    在 Gap 处切割基因组，获得 contig 水平的基因组
    - gap_size: 指定多少个连续 N 才算一个 gap；如果不指定，则认为任意连续 N 都是 gap。
    '''

    def split_sequence_by_gap(sequence, gap_size=None):
        """ 根据连续N的数量切割序列 """
        if gap_size is None:
            # 匹配任意 >=1 个连续N
            pattern = r"N+"
        else:
            # 匹配大于等于 gap_size 个连续N
            pattern = r"N{" + str(gap_size) + ",}"
        return [contig for contig in re.split(pattern, sequence) if len(contig) > 0]

    def write_fasta(output_file, contigs, seq_id):
        for i, contig in enumerate(contigs, 1):
            output_file.write(f">{seq_id}_contig{i}\n")
            for j in range(0, len(contig), 60):
                output_file.write(contig[j:j + 60] + "\n")

    def process_fasta(input_fasta, output_fasta, gap_size=None):
        sequences = fasta_read(input_fasta)
        with open(output_fasta, "w") as output_file:
            for seq_id, sequence in sequences.items():
                contigs = split_sequence_by_gap(sequence, gap_size)
                write_fasta(output_file, contigs, seq_id)

    base, ext = os.path.splitext(genome_fasta)
    output_fasta = base + ".contig.fa"
    process_fasta(genome_fasta, output_fasta, gap_size)