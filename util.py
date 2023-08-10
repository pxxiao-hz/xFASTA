#!/usr/bin/env python
'''
time: 2023-04-05
author: pxxiao
version: 1.0
description: 函数文件
'''
import readline, os, re
from Bio import SeqIO
from collections import defaultdict


def fasta_read(file_input):
    '''
    解析 FASTA 文件，存入字典 d
    :param file_input: .fasta
    :return: d: id as kye, seq as value
    '''
    d = {}
    for lines in open(file_input, "r"):
        if lines.startswith(">"):
            id = lines.strip().replace(">", "")
            d[id] = []
        else:
            d[id].append(lines.strip())
    for key, values in d.items():
        d[key] = "".join(values)
    return d


def split_fasta_based_bp(file_input, length):
    '''
    将 FASTA 文件分割成固定bp的小文件
    :param file_input: 待拆分的 FASTA 文件
    :param length: 按多少bp进行拆分
    :return: none
    '''
    d_fa = fasta_read(file_input)
    s = 1   # 多少条序列
    for key,values in d_fa.items():
        print("正在拆分第{}条染色体...".format(s))
        m = 1
        for i in range(0, len(values), length):
            if i + length < len(values):
                # 没到末尾
                out_split = open("{}_{}.fa".format(key.replace(">",""), str(m)), "w")
                out_split.write(key + "_" + str(m) + "\n" + values[i:i + length] + "\n")
                m += 1
            else:
                # 到末尾
                out_split = open("{}_{}.fa".format(key.replace(">", ""), str(m)), "w")
                out_split.write(key + "_" + str(m) + "\n" + values[i:] + "\n")
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
    d_fa = fasta_read(file_input)
    n = 1  # 控制文件包含的pep id个数
    s = 1  # 控制输出文件的名字
    for key, values in d_fa.items():
        out = "out_" + str(s)
        out = open(out, "a+")
        out.write(key + "\n" + values + "\n")
        out.close()
        n += 1
        if n <= number:
            continue
        else:
            s += 1
            n = 1
    return None

def get_Nnumber(BaseSum, Length, a):
    'get N 50,N 90, N 60...'
    N_pos = BaseSum / 100 * a
    Length.sort()
    Length.reverse()
    n = 0
    ValueSum = 0
    for value in Length:
        ValueSum += value
        n += 1
        if N_pos <= ValueSum:
            N_length = value
            N_number = n
            break
    return N_length, N_number


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
    out = open(file_output, "w")

    BaseSum, Length = 0, []
    ValueSum, N50 = 0, 0
    no_c, no_g, no_a, no_t, no_n = 0, 0, 0, 0, 0
    no_ctg_10k, no_ctg_50k, no_ctg_100k = 0, 0, 0

    for record in SeqIO.parse(open(file_input), "fasta"):
        BaseSum += len(record.seq)
        Length.append(len(record.seq))
        seq = record.seq.lower()
        no_c += seq.count('c')
        no_g += seq.count('g')
        no_a += seq.count('a')
        no_t += seq.count('t')
        no_n += seq.count('n')

        if len(record.seq) > 10000:
            no_ctg_10k += 1
        if len(record.seq) > 50000:
            no_ctg_50k += 1
        if len(record.seq) > 100000:
            no_ctg_100k += 1

    for i in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
        N_length = "N" + str(i)
        N_number = "L" + str(i)
        N_len, N_num = get_Nnumber(BaseSum, Length, i)
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
    ### contig 长度大于某个值
    # print('Ctg gt 10k:\t' + str(no_ctg_10k))
    print('Length > 50k\t' + str(no_ctg_50k))
    print('Length > 100k\t' + str(no_ctg_100k))
    # out.write('Ctg gt 10k:\t' + str(no_ctg_10k) + '\n')
    out.write('Length > 50k\t' + str(no_ctg_50k) + '\n')
    out.write('Length > 100k\t' + str(no_ctg_100k) + '\n')
    out.close()
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