#!/usr/bin/env python
'''
time: 2023-04-07
author: pxxiao
version: 1.0
description: 处理 GFF 文件的函数
'''
import os
from collections import defaultdict
from util import *

def func_get_ExonIntronInfo(file_input, splits):
    '''
    sort, extract
    :param file_input: gff file
    :param model: extract, sort
    :param splits: gff文件第九列以什么开头
    :return:
    '''

    # 获得exon, intron, gene的坐标信息，并将它们以bed格式输出到文件
    d_gene, d_exon, d_intron = gff_extract(file_input, splits)
    # output gene
    out = open("gene.bed", "w")
    for key, values in d_gene.items():
        # print(key)
        # print(values)
        chr = values[0]
        start = str(int(values[1]) - 1)
        end = values[2]
        strand = values[3]
        out.write(chr + "\t" + start + "\t" + end + "\t" + strand + "\t" + key + "\n")
    out.close()

    # output exon
    # print(d_exon)
    out = open("exon.bed", "w")
    for key, values in d_exon.items():
        for i in range(0, len(values), 2):
            print(key)
            print(d_gene[key][0])
            print(values[i])
            out.write(d_gene[key][0] + "\t" + str(int(values[i]) - 1) + "\t" + str(int(values[i + 1]))
                      + "\t" + d_gene[key][3] + "\t" + key + "\n")
    out.close()

    # output intron
    out = open("intron.bed", "w")
    for key, values in d_intron.items():
        for i in range(0, len(values), 2):
            out.write(d_gene[key][0] + "\t" + str(int(values[i])) + "\t" + str(int(values[i + 1]) - 1)
                      + "\t" + d_gene[key][3] + "\t" + key + "\n")
    out.close()
    return None


def func_get_LongTranscript(file_input):

    prefix = os.path.basename(file_input)
    if ".gff3" in prefix:
        output = prefix.replace(".gff3", ".LongTranscript.gff")
    elif ".gff" in prefix:
        output = prefix.replace(".gff", ".LongTranscript.gff")
    out = open(output, "w")
    d_gene = {}
    for lines in open(file_input, "r"):
        # filter
        if lines.startswith("#"): continue  # 过滤#
        line = lines.strip().split()
        if len(line) < 8: continue
        if line[2] not in ['mRNA', 'exon', 'CDS', 'gene', 'five_prime_UTR', 'three_prime_UTR',
                           'UTR_3', 'UTR_5']:
            continue
        if line[2] == "gene":
            info = line[8].split(";")
            gene_ID = info[0].replace("gene:", "")    ### 一会要取消注释的地方 ###
            # gene_ID = info[1].replace("Parent=","")
            gene_ID = gene_ID.replace("ID=", "")
            line = [i.replace("gene:", "") if "gene:" in i else i for i in line]
            d_gene[gene_ID] = []
            d_gene[gene_ID].append(line)
            print(gene_ID)
        elif line[2] == "mRNA":
            info = line[8].split(";")
            tx_ID = info[0].replace("transcript:", "")
            tx_ID = tx_ID.replace("ID=", "")
            # tx_ID = tx_ID.split(".")[0]
            line = [i.replace("transcript:", "") if "transcript:" in i else i for i in line]
            if "_" in tx_ID:
                gene_ID = tx_ID.split(".")[0]
                gene_ID = tx_ID.split("_")[0]     # 此时的tx_ID为gene_ID
            for x in info:
                if "Parent=" in x:
                    gene_ID = x.replace("Parent=", "")
                    gene_ID = gene_ID.replace("gene:", "")
                else:
                    continue
            if len(d_gene[gene_ID]) == 1:
                # 第一个mRNA
                # print(gene_ID)
                tx_long_ID = tx_ID
                mRNA_first_length = int(line[4]) - int(line[3]) + 1
                d_gene[gene_ID].append(line)
            else:
                # if gene_ID == "Vitvi01g00048":
                #     print(mRNA_first_length)
                #     print("Find")
                #     print(d_gene[gene_ID])
                mRNA_second_length = int(line[4]) - int(line[3]) + 1
                if mRNA_second_length > mRNA_first_length:
                    # 第二个转录本比较长，则保留第二个；否则，保留第一个
                    mRNA_first_length = mRNA_second_length
                    d_gene[gene_ID] = [d_gene[gene_ID][0]]
                    d_gene[gene_ID].append(line)
                    tx_long_ID = tx_ID
                else:
                    continue
        elif line[2] == "CDS" or line[2] == "exon":
            line = lines.strip().split()
            # 这里判断是不是最长的转录本，如果是最长的转录本，才将CDS，exon等信息输出；否则不输出
            for i in line[8].split(";"):
                if "Parent=" in i:
                    parent_ID = i.replace("Parent=","")
                    parent_ID = parent_ID.replace("transcript:", "")
            # print(parent_ID, "\t", tx_ID)
            if parent_ID == tx_long_ID:
                line = [i.replace("CDS:", "").replace("transcript:", "") if ("CDS:" in i) or ("exon:") in i
                        else i for i in line]
                d_gene[gene_ID].append(line)
            else:
                continue
        else:
            continue
        # output
    # print(d_gene)
    for key,values in d_gene.items():
        # key: gene_ID values: [[gene], [mRNA], [exon], [CDS] ...]
        for value in values:
            out.write("\t".join(value) + "\n")
    out.close()
    return None


def func_gff_NewGFF(file_input):
    '''
    有的gff不完整，缺少gene这一项，或者缺少exon，或者缺少exon；
    目的是根据mRNA补齐gene，根据CDS补齐exon。
    # version1.0: 不对缺少的东西判断，直接暴力将gene和exon同时输出，最后再去除冗余
    # version2.0: ???? (如何去除冗余？) ---- awk '!a[$0]++' file.txt > file_new.txt
    :param file_input: gff.file
    :return: none
    '''
    prefix = os.path.basename(file_input).replace(".gff3", "")
    prefix = prefix.replace(".gff", "")
    file_output = prefix + ".new.gff"
    out = open(file_output, "w")
    for lines in open(file_input, "r"):
        # filter
        if lines.startswith("#"): continue  # 过滤#
        line = lines.strip().split()
        if len(line) < 8: continue
        if line[2] not in ['mRNA', 'exon', 'CDS', 'gene', 'five_prime_UTR', 'three_prime_UTR',
                           'UTR_3', 'UTR_5']:
            continue
        # output
        if line[2] == "mRNA":
            gene_list = line[0:2] + ["gene"] + line[3:]
            out.write("\t".join(gene_list) + "\n")
            out.write("\t".join(line) + "\n")
        elif line[2] == "CDS":
            exon_list = line[0:2] + ["exon"] + line[3:]
            out.write("\t".join(exon_list) + "\n")
            out.write("\t".join(line) + "\n")
        else:
            out.write("\t".join(line) + "\n")
    out.close()
    cmd_gff_filter = 'awk \'!a[$0]++\' {} > xpx_temp.txt '.format(file_output)
    os.system(cmd_gff_filter)
    os.system("mv xpx_temp.txt {}".format(file_output))
    return None