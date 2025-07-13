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
from util import *
from Bio import SeqIO


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

    read_file = os.path.basename(file_input)
    pfx = read_file.replace('.gz', '')
    out = open(pfx+".qv.txt", "w")     # 文件存储每个read的平均质量值
    # # 文件可以处理压缩和未压缩
    # if pfx.endswith('.gz'):
    #     with gzip.open(pfx, 'rt') as f:
    #         fastq_lines = f.readlines()
    # else:
    #     with open(pfx, 'r') as f:
    #         fastq_lines = f.readlines()

    qualities = []      # 质量分数list
    with gzip.open(read_file, 'rt') as f:
        # 逐行读取文件内容
        for line in f:
            # 选择读取包含质量值的行
            if line.startswith('+'):
                qline = next(f).strip()
                quality_scores = [ord(q) - 33 for q in qline]  # 将 ASCII 码转换为实际质量分数
                qv = sum(quality_scores) / len(quality_scores)
                out.write(str(qv) + "\n")
                # print(qv)
                # print(quality_scores)
                qualities.append(qv)  # 添加到质量分数列表中
    fun_qv_percent(qualities, 20)
    fun_qv_percent(qualities, 25)
    fun_qv_percent(qualities, 30)
    fun_qv_percent(qualities, 35)
    fun_qv_percent(qualities, 40)
    print("Number of quality scores:", len(qualities))
    print("Minimum quality score:", min(qualities))
    print("Maximum quality score:", max(qualities))
    print("Average quality score:", sum(qualities) / len(qualities))
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
    list_readID = []
    # 输出文件名称
    name = os.path.basename(file_input)
    if "fq.gz" in name:
        file_output_name = name.replace("fq.gz", "")
    elif "fastq.gz" in name:
        file_output_name = name.replace("fastq.gz", "")
    out = file_output_name + "redu.fq.gz"
    # 处理
    with gzip.open(file_input, "rt") as handle, gzip.open(out, "wt") as output:
        for record in SeqIO.parse(handle, "fastq"):
            ID = record.id  # reads ID
            if ID not in list_readID:
                list_readID.append(ID)
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
    name = os.path.basename(file_input)
    if 'fq.gz' in name:
        file_output_name = name.replace('fq.gz', '')
    elif 'fastq.gz' in name:
        file_output_name = name.replace('fastq.gz', '')
    out = file_output_name + 'extract.fastq.gz'
    # 处理
    with gzip.open(file_input, "rt") as handle, gzip.open(out, "wt") as output:
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
    print(new_read_id)

    # 修改 Read ID
    record.id = new_read_id
    # record.id = record.id.split(' ')[0]
    # print(record)
    return record


def func_reads_ID_simplified(file_input):
    '''
    对 fastq reads ID 进行简化
    :param file_input: reads.fastq.gz
    :return: reads.simplified.fastq.gz
    '''
    name = os.path.basename(file_input)
    if 'fq.gz' in name:
        file_output_name = name.replace('fq.gz', '')
    elif 'fastq.gz' in name:
        file_output_name = name.replace('fastq.gz', '')
    out = file_output_name + 'extract.fastq.gz'
    # 处理
    records = []
    with gzip.open(file_input, 'rt') as handle, gzip.open(out, 'wt') as output:
        for record in SeqIO.parse(handle, "fastq"):
            print(record)
            modify_record = modify_read_id(record)
            records.append(modify_record)
            # print(modify_record)
            SeqIO.write(records, output, "fastq")
    return None


def func_reads_filter_ONT_reads(file_input, quality, length):
    '''
    对 ONTreads 进行过滤，
    :param file_input: 待输入的 ONT reads
    :param quality: 过滤的 reads 质量值，默认为 7
    :param length: 过滤的长度阈值
    :return: None
    '''
    name = os.path.basename(file_input)
    if 'fq.gz' in name:
        file_output_name = name.replace('fq.gz', '')
    elif 'fastq.gz' in name:
        file_output_name = name.replace('fastq.gz', '')
    ### 处理
    cmd_chopper = f'gunzip -c {file_input} | ~/tools/Anaconda3/envs/chopper/bin/chopper ' \
                  f'-q {quality} -l {length} | gzip > {file_output_name}q{quality}.l{length}.fastq.gz'
    subprocess.run(cmd_chopper, shell=True, close_fds=True)
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
