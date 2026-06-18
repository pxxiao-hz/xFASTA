#!/usr/bin/env python
'''
time: 2023-04-27
author: pxxiao
version: 1.0
description: 处理FASTQ文件的函数
'''
import os
import readline, subprocess
import gzip
import shlex
from util import *
from Bio import SeqIO


def _open_fastq_text(file_input):
    """Open plain or gzipped FASTQ input as text."""
    if file_input.endswith(".gz"):
        return gzip.open(file_input, "rt")
    return open(file_input, "r")


def _fastq_prefix(file_input):
    """Return the basename prefix used by FASTQ output files."""
    name = os.path.basename(file_input)
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return os.path.splitext(name)[0]


def fun_qv_percent(list, qv):
    '''
    获得大于qv的比例，例如qv20：87%
    :param list:
    :param qv:
    :return:
    '''
    greater_num = [i for i in list if i >= qv]
    greater_percent = (len(greater_num) / len(list)) * 100
    print("QV{}: ".format(qv), greater_percent)
    return None


def func_get_reads_phreads(file_input):

    pfx = _fastq_prefix(file_input)
    out = open(pfx + ".qv.txt", "w")     # 文件存储每个read的平均质量值

    count = 0
    total = 0
    min_qv = None
    max_qv = None
    threshold_counts = {20: 0, 25: 0, 30: 0, 35: 0, 40: 0}
    with _open_fastq_text(file_input) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            quality_scores = record.letter_annotations["phred_quality"]
            if not quality_scores:
                continue
            qv = sum(quality_scores) / len(quality_scores)
            out.write(str(qv) + "\n")
            count += 1
            total += qv
            min_qv = qv if min_qv is None else min(min_qv, qv)
            max_qv = qv if max_qv is None else max(max_qv, qv)
            for threshold in threshold_counts:
                if qv >= threshold:
                    threshold_counts[threshold] += 1
    out.close()

    if count == 0:
        print("Number of quality scores:", 0)
        return None

    for threshold in sorted(threshold_counts):
        print("QV{}: ".format(threshold), (threshold_counts[threshold] / count) * 100)
    print("Number of quality scores:", count)
    print("Minimum quality score:", min_qv)
    print("Maximum quality score:", max_qv)
    print("Average quality score:", total / count)
    return None


def func_reads_deduplication_based_ID(file_input):
    '''
    根据reads ID, 将具有重复的reads过滤掉
    :param file_input: fastq.gz, 压缩文件
    :return: 输出过滤之后的reads文件
    '''
    # !/user/bin/env python
    # encoding: utf-8
    '''
    Created on Jun 02, 2023
    @author: pxxiao
    version: v1
    '''
    seen_read_ids = set()
    out = _fastq_prefix(file_input) + ".redu.fq.gz"
    # 处理
    with _open_fastq_text(file_input) as handle, gzip.open(out, "wt") as output:
        for record in SeqIO.parse(handle, "fastq"):
            ID = record.id  # reads ID
            if ID not in seen_read_ids:
                seen_read_ids.add(ID)
                SeqIO.write(record, output, "fastq")
            else:
                continue
    return None


def func_reads_extract_based_ID(file_input, read_ID):
    '''
    根据reads ID，提取相应的fastq文件
    :param file_input: fastq.gz, 压缩文件
    :return: 输出提取的reads序列，还是fastq格式
    '''
    out = _fastq_prefix(file_input) + '.extract.fastq.gz'
    # 处理
    with _open_fastq_text(file_input) as handle, gzip.open(out, "wt") as output:
        for record in SeqIO.parse(handle, "fastq"):
            ID = record.id  # reads ID
            if ID == read_ID:
                SeqIO.write(record, output, "fastq")
                break
            else:
                continue
    return None


def modify_read_id(record):
    read_id_parts = record.id.split(' ')
    new_read_id = read_id_parts[0]

    # 修改 Read ID
    record.id = new_read_id
    record.name = new_read_id
    record.description = new_read_id
    return record


def func_reads_ID_simplified(file_input):
    '''
    对 fastq reads ID 进行简化
    :param file_input: reads.fastq.gz
    :return: reads.simplified.fastq.gz
    '''
    out = _fastq_prefix(file_input) + '.simplified.fastq.gz'
    # 处理
    with _open_fastq_text(file_input) as handle, gzip.open(out, 'wt') as output:
        for record in SeqIO.parse(handle, "fastq"):
            modify_record = modify_read_id(record)
            SeqIO.write(modify_record, output, "fastq")
    return None


def build_chopper_filter_command(file_input, quality, length, chopper="chopper"):
    '''
    Build the ONT read filtering shell pipeline.
    :param file_input: input FASTQ.gz file
    :param quality: minimum read quality
    :param length: minimum read length
    :param chopper: chopper executable or full path
    :return: (command, output_file)
    '''
    file_output_name = _fastq_prefix(file_input)
    output = f"{file_output_name}.q{quality}.l{length}.fastq.gz"
    command = (
        f"gunzip -c {shlex.quote(file_input)} | "
        f"{shlex.quote(chopper)} -q {quality} -l {length} | "
        f"gzip > {shlex.quote(output)}"
    )
    return command, output


def func_reads_filter_ONT_reads(file_input, quality, length, chopper="chopper"):
    '''
    对 ONTreads 进行过滤，
    :param file_input: 待输入的 ONT reads
    :param quality: 过滤的 reads 质量值，默认为 7
    :param length: 过滤的长度阈值
    :param chopper: chopper executable or full path
    :return: None
    '''
    cmd_chopper, _ = build_chopper_filter_command(file_input, quality, length, chopper)
    subprocess.run(cmd_chopper, shell=True, close_fds=True, check=True)
    return None


import os
import pandas as pd

def func_reads_stat_from_seqkit_stat(file_input):
    '''
    对 seqkit stat 结果进行整理，输出数据量和 reads 数量（按样本聚合）
    :param file_input: seqkit stats 输出文件（如 stat.all），需包含 file, num_seqs, sum_len 等列
    :return: pandas.DataFrame，包含每个样本的 reads 数和总数据量（Gb）
    '''
    # 读取数据
    df = pd.read_csv(file_input, sep=r"\s+", engine="python")

    # 提取样本名（去掉 _1/_2 和扩展名）
    df['sample'] = df['file'].str.replace(r'_([12])\.f(ast)?q(\.clean)?\.gz$', '', regex=True)

    # 去掉逗号，并将数值字段转为数字
    df['num_seqs'] = df['num_seqs'].str.replace(',', '').astype(int)
    df['sum_len'] = df['sum_len'].str.replace(',', '').astype(int)

    # 每个样本汇总：双端总reads数、总测序量（Gb）
    summary = df.groupby('sample').agg({
        'num_seqs': 'sum',
        'sum_len': 'sum'
    }).reset_index()

    # 转换为 Gb 并保留两位小数
    summary['sum_len_Gb'] = (summary['sum_len'] / 1e9).round(2)

    # 格式化输出
    summary = summary[['sample', 'num_seqs', 'sum_len_Gb']]
    summary.columns = ['Sample', 'Total_Reads', 'Total_Data(Gb)']

    # 直接输出到终端
    print(summary)

    return None
