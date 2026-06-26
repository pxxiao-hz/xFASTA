#!/usr/bin/env python
'''
time: 2023-04-27
author: pxxiao
version: 1.0
description: 处理端粒信息的函数
'''
import subprocess
import sys


def build_telomere_command(file_input, length, minrptnum, output, telomere_script, python_executable=None):
    '''
    Build the external telomere finder command.
    :param file_input: input FASTA file
    :param length: sequence length searched at both chromosome ends
    :param minrptnum: minimum repeat number
    :param output: output file
    :param telomere_script: path to 02_find_telomere.py or compatible script
    :param python_executable: Python executable used to run the script
    :return: command list for subprocess.run
    '''
    if python_executable is None:
        python_executable = sys.executable
    return [
        python_executable,
        telomere_script,
        "-i", file_input,
        "-l", str(length),
        "-m", str(minrptnum),
        "-o", output,
    ]


def func_telomere_info(input, length, minrptnum, output, telomere_script="02_find_telomere.py"):
    cmd = build_telomere_command(input, length, minrptnum, output, telomere_script)
    subprocess.run(cmd, close_fds=True, check=True)
    return None
