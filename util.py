#!/usr/bin/env python
'''
time: 2023-04-05
author: pxxiao
version: 1.0
description: 函数文件
'''
import os
import subprocess
from collections import defaultdict


### 判断文件是否存在
def check_file_in_path(file, cmd):
    '''
    Run a command only when the expected output file is missing.
    :param file: expected output file path
    :param cmd: shell command used to create the file
    '''
    if os.path.exists(file):
        print(f'[info] {file} exists, CMD: {cmd}; PASS!')
    else:
        subprocess.run(cmd, shell=True, check=True)
    return None


def fasta_read(file_input):
    '''
    解析 FASTA 文件，存入字典 d
    :param file_input: .fasta
    :return: d: id as kye, seq as value
    '''
    return dict(iter_fasta_records(file_input))


def iter_fasta_records(file_input):
    '''
    Stream FASTA records one at a time.
    :param file_input: .fasta/.fa file
    :return: yields (record_id, sequence)
    '''
    record_id = None
    seq_parts = []
    with open(file_input, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, "".join(seq_parts)
                record_id = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if record_id is not None:
            yield record_id, "".join(seq_parts)


def write_fasta_record(handle, record_id, sequence, line_width=None):
    '''
    Write one FASTA record with wrapped sequence lines.
    :param handle: writable file handle
    :param record_id: FASTA record id without leading ">"
    :param sequence: sequence string
    :param line_width: optional sequence wrap width; None preserves the previous one-line output
    '''
    handle.write(f">{record_id}\n")
    if line_width is None:
        handle.write(sequence + "\n")
    else:
        for i in range(0, len(sequence), line_width):
            handle.write(sequence[i:i + line_width] + "\n")


def split_fasta_based_bp(file_input, length):
    '''
    将 FASTA 文件分割成固定bp的小文件
    :param file_input: 待拆分的 FASTA 文件
    :param length: 按多少bp进行拆分
    :return: none
    '''
    if length <= 0:
        raise ValueError("length must be greater than 0")

    s = 1   # 多少条序列
    for key, values in iter_fasta_records(file_input):
        print("正在拆分第{}条染色体...".format(s))
        m = 1
        for i in range(0, len(values), length):
            seq = values[i:i + length]
            with open("{}_{}.fa".format(key.replace(">", ""), str(m)), "w") as out_split:
                write_fasta_record(out_split, key + "_" + str(m), seq)
            m += 1
        s += 1
    print("染色体拆分完成！")
    return None


def split_fasta_based_number(file_input, number):
    '''
    将 FASTA 文件分割成小文件，每个文件包含的序列个数为 number
    :param file_input: 待拆分的 FASTA 文件
    :param number: 每个文件包含的序列个数
    :return: none
    '''
    if number <= 0:
        raise ValueError("number must be greater than 0")

    output_index = 1
    records_in_output = 0
    out = None
    try:
        for key, values in iter_fasta_records(file_input):
            if out is None or records_in_output >= number:
                if out is not None:
                    out.close()
                out = open("out_" + str(output_index), "w")
                output_index += 1
                records_in_output = 0
            write_fasta_record(out, key, values)
            records_in_output += 1
    finally:
        if out is not None:
            out.close()
    return None

def get_Nnumber(BaseSum, Length, a):
    'get N 50,N 90, N 60...'
    return _get_nnumber_from_sorted(BaseSum, sorted(Length, reverse=True), a)


def _get_nnumber_from_sorted(base_sum, sorted_lengths, percentage):
    """Calculate Nx/Lx from a descending sequence-length list."""
    N_pos = base_sum / 100 * percentage
    n = 0
    ValueSum = 0
    for value in sorted_lengths:
        ValueSum += value
        n += 1
        if N_pos <= ValueSum:
            return value, n
    raise ValueError("Length must contain at least one sequence")


def func_n50(file_input):
    '''
    计算基因的N50，长度等信息
    :param file_input: .fata 文件
    :return: none
    '''
    # 将结果输出到一个新的文件当中
    prefix = os.path.basename(file_input)
    prefix = prefix.replace(".fasta", "")
    prefix = prefix.replace(".fa", "")
    file_output = prefix + ".stat.txt"  # 存储相关信息
    BaseSum, Length = 0, []
    no_ctg_50k, no_ctg_100k = 0, 0
    len_ctg_50k, len_ctg_100k = 0, 0
    no_ctg_lt_100k = 0
    len_ctg_lt_100k = 0

    for _, sequence in iter_fasta_records(file_input):
        seq_length = len(sequence)
        BaseSum += seq_length
        Length.append(seq_length)
        if seq_length > 50000:
            no_ctg_50k += 1
            len_ctg_50k += seq_length
        if seq_length > 100000:
            no_ctg_100k += 1
            len_ctg_100k += seq_length
        if seq_length < 100000:
            no_ctg_lt_100k += 1
            len_ctg_lt_100k += seq_length

    if not Length:
        raise ValueError("FASTA file contains no sequences")

    sorted_lengths = sorted(Length, reverse=True)
    with open(file_output, "w") as out:
        for i in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
            N_length = "N" + str(i)
            N_number = "L" + str(i)
            N_len, N_num = _get_nnumber_from_sorted(BaseSum, sorted_lengths, i)
            if i in [10, 50, 90]:
                print(N_length + "\t" + str(N_len))
                print(N_number + "\t" + str(N_num))
            out.write(N_length + "\t" + str(N_len) + "\n")
            out.write(N_number + "\t" + str(N_num) + "\n")
        print('Seq. Num.\t' + str(len(Length)))
        print('Seq. Min\t' + str(min(Length)))
        print('Seq. Max\t' + str(max(Length)))
        print('Seq. Total\t' + str(BaseSum))
        out.write('Seq. Num.\t' + str(len(Length)) + "\n")
        out.write('Seq. Min\t' + str(min(Length)) + "\n")
        out.write('Seq. Max\t' + str(max(Length)) + "\n")
        out.write('Seq. Total\t' + str(BaseSum) + "\n")
        print(f'Length < 100k Number\t{str(no_ctg_lt_100k)}')
        print(f'Length < 100k Length\t{str(len_ctg_lt_100k)}')
        out.write(f'Length < 100k Num.\t{str(no_ctg_lt_100k)}\n')
        out.write(f'Length < 100k Len.\t{str(len_ctg_lt_100k)}\n')
        print('Length > 50k Number\t' + str(no_ctg_50k))
        print(f'Length > 50k Length\t{str(len_ctg_50k)}')
        print('Length > 100k Number\t' + str(no_ctg_100k))
        print(f'Length > 100k Length\t{str(len_ctg_100k)}')
        out.write('Length > 50k Num.\t' + str(no_ctg_50k) + '\n')
        out.write(f'Length > 50k Len.\t{str(len_ctg_50k)}\n')
        out.write('Length > 100k Num.\t' + str(no_ctg_100k) + '\n')
        out.write(f'Length > 100k Len.\t{str(len_ctg_100k)}\n')
    return None


def gff_sort(file_input):
    return None


def gff_extract(file_input, splits):
    '''
    提取GFF的exon, intron, gene信息，存入字典里面
    :param file_input:
    :return:
    '''
    d_gene = defaultdict(list)
    d_exon = defaultdict(list)
    d_intron = defaultdict(list)
    for lines in open(file_input):
        # filter
        if lines.startswith("#"): continue  # 过滤#
        line = lines.strip().split()
        if len(line) < 8: continue
        if line[2] not in ['mRNA', 'exon', 'CDS', 'gene', 'five_prime_UTR', 'three_prime_UTR',
                           'UTR_3', 'UTR_5']:
            continue
        chrom = line[0]
        type = line[2]
        start = line[3]
        end = line[4]
        strand = line[6]
        phase = line[7]
        info = line[8]
        if type == "gene":
            info = info.split(";")
            gene_id = info[0].replace(splits, "")    # marked这里可能不一样呀
            # evm 不一致
            if "evm.TU" in gene_id:     # 如果是evm注释的话，有的项目没有改基因名称，默认的是：基因ID evm.TU.xx; exon: evm.model.xx
                gene_id = gene_id.replace("evm.TU", "evm.model")
            d_gene[gene_id] = [chrom, start, end, strand]  # chr, start, end, strand; 为什么要强调strand，因为涉及到基因上下游的问题。
        elif type == "exon":
            info = info.split(";")
            for x in info:
                if "Parent=" in x:
                    gene_id = x.replace("Parent=", "")
                    # evm 不一致
                    if "evm.TU" in gene_id:
                        gene_id = gene_id.replace("evm.TU", "evm.model")
                else:
                    continue
                    # gene_id = info[0].split(".1.exon")
                    # gene_id = gene_id[0].replace("ID=", "")
            gene_id = gene_id.replace("transcript:", "")
                    # gene_id = gene_id.replace("Parent=", "")
            gene_id = gene_id.split(".")[0]       # 修改的地方！！！
            # gene_id = gene_id.split("_")[0]       # 修改的地方！！！

            # gene_id = re.findall('(?<=ID=).+?(?=;)',line[8])
            if strand == "+":
                d_exon[gene_id].append(start)
                d_exon[gene_id].append(end)
            elif strand == "-":
                d_exon[gene_id].insert(0, end)  # insert：在列表任意索引位置插入元素
                d_exon[gene_id].insert(0, start)
    for key, values in d_exon.items():
        for i in range(1, len(values) - 1):
            d_intron[key].append(d_exon[key][i])
    return d_gene, d_exon, d_intron


def dic_sort(d):
    '''
    根据字典的值，进行排序
    :param d: 字典，key任意，值为数字
    :return: d_sort
    '''
    l_values = []
    l_key = []
    for key, values in d.items():
        l_key.append(key)
        l_values.append(values)
    return d_sort


def rev_seq(seq):
    '''
    序列反向互补
    :param seq:
    :return:
    '''
    complement = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn",
    )
    return seq.translate(complement)[::-1]
