#!/usr/bin/env python
'''
time: 2023-04-27
author: pxxiao
version: 1.0
description: 处理FASTQ文件的函数
'''
import subprocess


def func_telomere_info(input, length, minrptnum, output):

    cmd = f'python ~/test/scripts/04_T2T/03_evaluate_tools/02_find_telomere.py -i {input} -l {length}' \
          f' -m {minrptnum} -o {output}'
    subprocess.run(cmd, shell=True, close_fds=True)
    return None