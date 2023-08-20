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