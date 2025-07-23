#!/usr/bin/env python
'''
time: 2023-04-05
author: pxxiao
version: 1.0
description: 处理FASTA文件的函数
'''
import os.path
import readline, subprocess
from util import *


def output_based_FileInput(file_input, sufx):
    '''
    根据file_input的名称，定义fileoutput
    :param file_input:
    :param sufx:
    :return:
    '''
    basename = os.path.basename(file_input)
    if ".fasta" in basename:
        basename = basename.replace(".fasta", sufx)
    else:
        basename = basename.replace(".fa", sufx)
    out = open(basename, "w")
    return out


def func_split(choice, file_input, length, number):
    print("用于处理 FASTA 文件.")
    if choice == "length":
        print("\t拆分 FASTA 文件：按照 Length")
    else:
        print("\t拆分 FASTA 文件：按照 Number")

    if choice == 'length':
        split_fasta_based_bp(file_input, length)
    elif choice == 'number':
        split_fasta_based_number(file_input, number)
    else:
        print("Invalid choice.")


def func_stat(file_input, fun):
    print("正在统计 {} 的相关信息.".format(file_input))
    # print("\t1. child function 1")
    # print("\t2. child function 2")

    # choice = input("Choose an option: ")

    if fun == 'n50':
        func_n50(file_input)
    # elif choice == '2':
    #     func_two_child_two()
    else:
        print("Invalid choice.")

def func_three():
    print("This is function three.")
    print("\t1. child function 1")
    print("\t2. child function 2")

    choice = input("Choose an option: ")

    if choice == '1':
        func_three_child_one()
    elif choice == '2':
        func_three_child_two()
    else:
        print("Invalid choice.")
    return None


def extract_based_keywords(file_input, keywords, greedy):
    '''
    根据关键字，提取序列
    :param file_input: 后缀是.fasta or .fa
    :param keywords: 关键字，可以是完整的序列ID，也可以是ID的一部分
    :return: none
    '''
    l_keywords = keywords.split(',')
    out = output_based_FileInput(file_input, ".extract.fasta")
    d = fasta_read(file_input)
    for key,values in d.items():
        for k in l_keywords:
            if greedy == "True":
                if k == key:
                    out.write(">" + key + "\n" + values + "\n")
                else:
                    continue
            elif greedy == "False":
                if k in key:
                    out.write(">" + key + "\n" + values + "\n")
                else:
                    continue

    out.close()
    return None


def extract_based_length(file_input, seqlength):
    '''
    根据序列长度，提取序列
    :param file_input: 后缀是.fasta or .fa；可以是chr genome，也可以是contig genome
    :param seqlength: 序列的长度
    :return:
    '''
    out = output_based_FileInput(file_input, '.extract.fasta')
    d = fasta_read(file_input)
    for key, values in d.items():
        if len(values) > seqlength:
            out.write('>' + key + '\n' + values + '\n')
        else:
            continue
    out.close()
    return None


def extract_based_coords(file_input, coords):
    '''
    根据坐标信息，提取序列;
    这里的位置信息是 1-based。
    >Chr1
    123456789
    >Chr1:4-8
    45678
    :param file_input: 后缀是.fasta or .fa
    :param coords: chr-start-end
    :return:
    '''
    out = output_based_FileInput(file_input, '.extract.fasta')
    d = fasta_read(file_input)
    coords_info = coords.split('-')
    if coords_info[2] != '':
        chr = coords_info[0]
        start = int(coords_info[1]) - 1
        end = int(coords_info[2])
    elif coords_info[2] == '':
        chr = coords_info[0]
        start = int(coords_info[1]) - 1
        end = len(d[chr])
        # print(end)
    for key, values in d.items():
        if key == chr:
            seq = values[start:end]
            out.write('>' + chr + ':' + str(start + 1) + '-' + str(end) + '\n' + seq + '\n')
    out.close()
    return None


def func_get_fasta_length(file_input):
    '''
    获得FASTA文件每一条序列的长度
    :param file_input:
    :return:
    '''
    print(f"对{file_input}进行序列长度统计")
    out = output_based_FileInput(file_input, ".falen")
    out_sort = output_based_FileInput(file_input, ".sort.falen")
    d_fa = fasta_read(file_input)
    d_fa_len = {}
    # 输出 序列ID 序列长度（没有排序）
    for key,values in d_fa.items():
        d_fa_len[key] = len(values)
        out.write(key + "\t" + str(len(values)) + "\n")
    out.close()
    # 输出 序列ID 序列长度（长度逆向排序)
    d_fa_len_sort = sorted(d_fa_len.items(), key=lambda item:item[1], reverse=True) # 按照value进行排序
    for i in d_fa_len_sort:
        out_sort.write(i[0] + "\t" + str(i[1]) + "\n")
    out_sort.close()
    return None


def func_fasta_IDsimplified(file_input, s):
    '''
    fasta文件ID简化；如果有TAB，则取第一列；如果有
    :param file_input:
    :param s:
    :return:
    '''

    d_fa = {}
    basename = os.path.basename(file_input)
    if ".fasta" in basename:
        basename = basename.replace(".fasta", ".IDsimp.fasta")
    else:
        basename = basename.replace(".fa", ".IDsimp.fa")
    out = open(basename, "w")
    for lines in open(file_input, "r"):
        if lines.startswith(">"):
            # 如果TAB存在
            if lines.strip().count('\t') > 0:
                line = lines.strip().split()
                id = line[0]
            elif s in lines:
                if s == "." or s == "_":
                    line = lines.strip().split(s)
                    id = line[0]
                elif s == "-":
                    line = lines.strip().split(s)
                    id = ".".join(line[:-1])
            elif s == "None":
                line = lines.strip().split()
                id = line[0]
            else:
                line = lines.strip()
                id = line
        else:
            if id not in d_fa.keys():
                d_fa[id] = []
                d_fa[id].append(lines.strip())
            else:
                d_fa[id].append(lines.strip())
    for key, values in d_fa.items():
        d_fa[key] = "".join(values)

    for key, values in d_fa.items():
        out.write(key + "\n" + values + "\n")
    out.close()

    # cmd1 = "seqkit seq -w 80 {} > temp".format(basename)
    # cmd2 = " mv temp {}".format(basename)
    # subprocess.run(cmd1, shell=True, close_fds=True)
    # subprocess.run(cmd2, shell=True, close_fds=True)

    # # 有的文件中会有“*”字符，这时候可以选择删掉
    # cmd = " sed -i 's/*//g' {}".format(basename)
    # subprocess.run(cmd, shell=True, close_fds=True)
    return None


def func_fasta_compare(newFA, oldFA):
    '''
    比较两个FASTA文件是否一致
    :param newFA: new FASTA file
    :param oldFA: old FASTA file
    :return: none
    '''
    d_new = fasta_read(newFA)
    d_old = fasta_read(oldFA)
    for key,value in d_new.items():
        if key in d_old.keys():
            if value == d_old[key]:
                print("{} is same".format(key))
            else:
                print("{} is not same".format(key))
        else:
            print("{} not in {}".format(key, oldFA))
            print("{} is not same".format(key))
    return None


def func_fasta_toPhy(file_input):
    '''
    将FASTA转换为phy格式
    :param file_input: FASTA file
    :return: none
    '''
    out = output_based_FileInput(file_input, ".phy")
    with open(file_input, 'r') as fin:
        sequences = [(m.group(1), ''.join(m.group(2).split()))
                     for m in re.finditer(r'(?m)^>([^ \n]+)[^\n]*([^>]*)', fin.read())]
    out.write('%d %d\n' % (len(sequences), len(sequences[0][1])))
    for item in sequences:
        out.write('%-20s %s\n' % item)
    return None


def func_fasta_count(file_input, keywords):
    '''
    统计一个文件中，keywords出现的频率
    :param file_input: 给定的文件
    :param keywords: 关键词
    :return: none
    '''
    frequency = 0
    for lines in open(file_input, "r"):
        line = lines.strip()
        length = len(line)
        # print(f"长度为：{length}")
        frequency += line.count(keywords)
        indices = [i for i in range(len(line)) if line.startswith(keywords, i)]
    print(f"\"{keywords}\" 出现的频率为：{frequency} 次。")
    # print(f"出现的索引值为：{indices}")
    return None


def genome_change_chr_name(genome_fasta, file_chr_list):
    '''
    根据提供的染色体列表，更改染色体名称
    :param genome_fasta: 将要更改的基因组文件
    :param file_chr_list: 染色体替换文件：第一列原始 ID  第二列新的 ID
    :return:
    '''
    dic_genome_fasta = fasta_read(genome_fasta)
    dic_ID = {}

    # 读取染色体替换表
    for line in open(file_chr_list):
        line = line.strip().split()
        if len(line) >= 2:  # 确保至少有两列
            dic_ID[line[0]] = line[1]

    # 输出新文件
    out = output_based_FileInput(genome_fasta, '.ID.change.fasta')
    for old_id, sequence in dic_genome_fasta.items():
        if old_id in dic_ID:
            # 写入新ID + 原序列
            out.write(f'>{dic_ID[old_id]}\n{sequence}\n')
        else:
            # 保留未匹配的染色体
            out.write(f'>{old_id}\n{sequence}\n')
    out.close()


    # out = output_based_FileInput(genome_fasta, '.ID.change.fasta')
    # for lines in open(file_chr_list):
    #     line = lines.strip().split()
    #     dic_ID[line[0]] = line[1]
    # for key, value in dic_ID.items():
    #     out.write('>' + dic_ID[key] + '\n' + dic_genome_fasta[key] + '\n')
    #
    # for key, value in dic_genome_fasta.items():
    #     if key in dic_ID:
    #         out.write('>' + dic_genome_fasta[key] + '\n')
    #     else:
    #         continue

    # for key, values in dic_genome_fasta.items():
    #     if key in dic_ID.keys():
    #         out.write('>' + key + '\n' + dic_genome_fasta[dic_ID[key]] + '\n')
    #     else:
    #         out.write('>' + key + '\n' + dic_genome_fasta[key] + '\n')
    # out.close()
    return None


def genome_reverse_some_chr(genome_fasta, file_chr_list):
    '''
    根据提供的染色体列表，反向互补这些染色体的序列
    :param genome_fasta: genome_fasta
    :param file_chr_list: 需要反向互补的染色体列表: each row each ID
    :return:
    '''

    dic_genome_fasta = fasta_read(genome_fasta)
    out = output_based_FileInput(genome_fasta, '.rev.genome.fasta')
    list_rev_ID = []
    for lines in open(file_chr_list):
        line = lines.strip().split()
        list_rev_ID.append(line[0])
    for key, values in dic_genome_fasta.items():
        if key not in list_rev_ID:
            out.write('>' + key + '\n' + dic_genome_fasta[key] + '\n')
        else:
            new_seq = rev_seq(dic_genome_fasta[key])
            out.write('>' + key + '\n' + new_seq + '\n')
    out.close()
    return None


def genome_karyotype(genome_fasta):
    '''
    生成核型文件（karyotype），格式为：染色体    1    长度
    :param genome_fasta: 输入基因组FASTA文件
    :return: None（输出文件）
    '''
    print(f"正在处理 {genome_fasta}，生成核型文件...")
    d_fa = fasta_read(genome_fasta)
    out_karyotype = output_based_FileInput(genome_fasta, ".karyotype.txt")

    for chrom, seq in d_fa.items():
        length = len(seq)
        out_karyotype.write(f"{chrom}\t1\t{length}\n")

    out_karyotype.close()
    print(f"输出完成")
    return None


def fasta_extract_based_length(genome_fasta, length_cutoff):
    '''
    根据 Length cutoff 提取相应的序列，长于 cutoff 提取
    :param genome_fasta:
    :param length_cutoff:
    :return:
    '''
    dic_genome_fasta = fasta_read(genome_fasta)
    out = output_based_FileInput(genome_fasta, f'.{length_cutoff}.fa')
    for key, values in dic_genome_fasta.items():
        if len(values) > length_cutoff:
            out.write('>' + key + '\n' + values[0] + '\n')
        else:
            continue
    out.close()
    return None