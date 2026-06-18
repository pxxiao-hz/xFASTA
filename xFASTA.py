#!/usr/bin/env python3
'''
time: 2023-04-05
author: pxxiao
version: 1.0
description: 用于日常的 FASTA 文件处理，其他功能待定
'''

import argparse


def build_parser():
    """Build the CLI parser without importing heavy bioinformatics modules."""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='pxx 的分析小工具')
    subparsers = parser.add_subparsers(dest='function', metavar='')

    # FASTA subcommand.
    func_one_parser = subparsers.add_parser('fasta', formatter_class=argparse.RawTextHelpFormatter,
                                            help='处理 FASTA 文件')
    func_one_parser.add_argument('-i', '--input', type=str, required=True,
                                 help='file input: .fasta or .fa')
    func_one_parser.add_argument('-m', '--model', dest='model', required=True,
                                 choices=['length', 'number', 'extractk', 'extractl', 'extractc', 'falen',
                                          'IDsimp', 'compare', 'toPhy', 'count', 'changeid', 'reverseq',
                                          'karyotype'], metavar='MODEL',
                                 help='Process on .fasta or .fa \n'
                                      'length: Split File Based on Length \n'
                                      'number: Split File Based on Number \n'
                                      'extractk: Sequence Extraction Based on Keywords,\n'
                                      '\t -k required \n'
                                      'extractl: Sequence Extraction Based on seq length,\n'
                                      '\t -l required \n'
                                      'extractc: Sequence Extraction Based on coords, \n'
                                      '\t -k required, e.g.: chr1-10-20; chr1-10-, means 10 - end of chr1 \n'
                                      '\t Note: 1-based in this script \n'
                                      'falen: Length of each sequence \n'
                                      'IDsimp: Simplified ID Name,\n'
                                      '\t -s required \n'
                                      'compare: Compare the contents of two FASTA files for equality \n'
                                      'toPhy: Convert FASTA format to phy format \n'
                                      'count: Count the number of strings \n'
                                      'changeid: change fasta id based id list, -c required \n'
                                      'reverseq: reverse some chr seq based id list, -c required \n'
                                      'karyotype: Karyotype information, such as: chr01\t1\tlength')
    func_one_parser.add_argument('-l', '--length', type=int, default=20000000, required=False,
                                 help='拆分成多少bp的小文件 or 根据序列的长度从基因组文件提取序列 [20000000]')
    func_one_parser.add_argument('-n', '--number', type=int, default=50, required=False,
                                 help='拆分成多个文件，每个文件包含number个序列 [50]')
    func_one_parser.add_argument('-k', '--keywords', type=str, required=False,
                                 help='Keywords in the Sequence to be Extracted or Count.\n'
                                      'For Extracted, Chr1 or Chr1,Chr2,Chr3')
    func_one_parser.add_argument('-g', '--greedy', type=str, required=False, default='False',
                                 choices=['True', 'False'],
                                 help='Exact ID match for extractk, True or False [False]')
    func_one_parser.add_argument('-s', '--separator', type=str, required=False,
                                 help='Separator of sequence ID; "None" == "TAB"')
    func_one_parser.add_argument('-x', '--newFASTA', type=str, required=False,
                                 help='File 2 of two FASTA files for compare')
    func_one_parser.add_argument('-p', '--toPhy', type=str, required=False,
                                 help='FASTA file to Phy')
    func_one_parser.add_argument('-c', '--chrlist', type=str, required=False,
                                 help='File chr list, first col is oldID, second col is newID, TAB separated')

    # Genome assembly statistics.
    func_stat_parser = subparsers.add_parser('stat', formatter_class=argparse.RawTextHelpFormatter,
                                             help='统计基因组的相关信息')
    func_stat_parser.add_argument('-i', '--input', type=str, required=True,
                                  help='输入文件，待处理的 FASTA 文件')
    func_stat_parser.add_argument('-n', '--n50', type=str, required=False,
                                  help='基因组组装结果的统计')

    # GFF subcommand.
    func_gff_parser = subparsers.add_parser('gff', formatter_class=argparse.RawTextHelpFormatter,
                                            help='Perform operations on GFF file')
    func_gff_parser.add_argument('-i', '--input', type=str, required=True,
                                 help='file input: GFF file')
    func_gff_parser.add_argument('-m', '--model', type=str, required=True,
                                 choices=['sort', 'extract', 'GetLongTrans', 'NewGff'], metavar='MODEL',
                                 help='GFF sub command, including sort, extract, GetLongTrans\n'
                                      'sort: based on chr and start \n'
                                      'extract: extract exon, intron information, output format: bed\n'
                                      'GetLongTrans: only retain longest information, output format: gff\n'
                                      'NewGff: if gff not have "gene" mark, please use NewGff')
    func_gff_parser.add_argument('-s', '--split', type=str, required=False, default='ID=',
                                 help='the 9th column of gff file, maybe startswith "ID=", or "Parent=')

    # Gap subcommand.
    func_gap_parser = subparsers.add_parser('gap', formatter_class=argparse.RawTextHelpFormatter,
                                            help='Genome gap statistics')
    func_gap_parser.add_argument('-i', '--input', type=str, required=True,
                                 help='finle input: Genome fasta file')
    func_gap_parser.add_argument('-m', '--model', type=str, required=True,
                                 choices=['coord'], metavar='MODEL',
                                 help='gap sub coomand, including coord\n'
                                      'coord: Obtaining the location information of gaps, '
                                      'NNNN... or nnnn... was defined as gap')

    # FASTQ/read subcommand.
    func_read_parser = subparsers.add_parser('read', formatter_class=argparse.RawTextHelpFormatter,
                                             help='Reads information')
    func_read_parser.add_argument('-i', '--input', type=str, required=True,
                                  help='The input FASTQ file, can be in .fastq or .fastq.gz format')
    func_read_parser.add_argument('-m', '--model', type=str, required=True,
                                  choices=['phreads', 'deduplication', 'extract', 'IDsimp', 'filterONT', 'stat'],
                                  metavar='MODEL',
                                  help='read sub command, including phreads\n'
                                       'phreads: Parse the quality value of reads\n'
                                       'deduplication: remove reads based on reads ID, if reads ID not uniq\n'
                                       'extract: Extract reads sequences; extensions: fq.gz or fastq.gz\n'
                                       'IDsimp: ID simplified\n'
                                       'filterONT: Filter ONT reads based on quality or length [l 100000] [q 7]\n'
                                       'stat: the results of seqkit stat')
    func_read_parser.add_argument('-k', '--readid', type=str, required=False, metavar='READ_ID',
                                  help='Extract reads sequence based the reads ID')
    func_read_parser.add_argument('-q', '--ONTquality', type=int, required=False, metavar='QUALITY',
                                  help='Filter ONT reads based quality, e.g., 7')
    func_read_parser.add_argument('-l', '--ONTlength', type=int, required=False, metavar='LENGTH',
                                  help='Filter ONT reads based length, e.g., 100000')
    func_read_parser.add_argument('--chopper', type=str, required=False, default='chopper',
                                  help='chopper executable or full path for read filterONT [chopper]')

    # Telomere subcommand.
    func_telomere_parser = subparsers.add_parser('telomere', formatter_class=argparse.RawTextHelpFormatter,
                                                 help='Telomere information')
    func_telomere_parser.add_argument('-i', '--input', type=str, required=True,
                                      help='FASTA file')
    func_telomere_parser.add_argument('-l', '--length', type=int, required=False, default=150000,
                                      help='Length of both end chromosome [150000]')
    func_telomere_parser.add_argument('-m', '--minrptnum', type=int, required=False, default=100,
                                      help='Minimum repeat number of each telomere [100]')
    func_telomere_parser.add_argument('-o', '--output', type=str, required=True,
                                      help='output')
    func_telomere_parser.add_argument('--telomere-script', type=str, required=False,
                                      default='02_find_telomere.py',
                                      help='path to 02_find_telomere.py [02_find_telomere.py]')

    # Genome workflow subcommand.
    func_genome_parser = subparsers.add_parser('genome', formatter_class=argparse.RawTextHelpFormatter,
                                               help='Genome operation')
    func_genome_parser.add_argument('-m', '--model', type=str, required=True,
                                    choices=['survey', 'telomere', 'split_by_gap'],
                                    help='genome sub command, including: \n'
                                         '\tsurvey: genome survey based on Illumina reads, -p required \n'
                                         '\ttelomere: genome telomere identified, [-i/--input_genome] required \n'
                                         '\tsplit_by_gap: Cut at the gap, obtain contig-level genome, '
                                         '[-i/--input_genome] required \n')
    func_genome_parser.add_argument('-k', '--kmersize', type=int, required=False, default=21,
                                    help='kmer size [21]')
    func_genome_parser.add_argument('-c', '--count', type=int, required=False, default=10000000,
                                    help='High count value of histogram [10000000]')
    func_genome_parser.add_argument('-b', '--poly', type=int, required=False, default=2,
                                    help='poly [2]')
    func_genome_parser.add_argument('-t', '--threads', type=int, required=False, default=16,
                                    help='threads [16]')
    func_genome_parser.add_argument('-p', '--path', type=str, required=False,
                                    help='path of fastq reads, decompress format ')
    func_genome_parser.add_argument('-l', '--length', type=int, required=False, default=150,
                                    help='reads length [150] \n')
    func_genome_parser.add_argument('-g', '--gap_size', type=int, required=False, default=None,
                                    help='Minimum number of consecutive N/n used to split sequence '
                                         '(default: any length)')
    func_genome_parser.add_argument('-i', '--inputGenome', help='The genome fasta file')
    func_genome_parser.add_argument('--telomere-script', type=str, required=False,
                                    default='02_find_telomere.py',
                                    help='path to 02_find_telomere.py for genome telomere [02_find_telomere.py]')
    func_genome_parser.add_argument('--genomescope1-script', type=str, required=False,
                                    default='genomescope.R',
                                    help='GenomeScope v1 R script path for genome survey [genomescope.R]')
    func_genome_parser.add_argument('--genomescope2-script', type=str, required=False,
                                    default='genomescope.R',
                                    help='GenomeScope v2 R script path for genome survey [genomescope.R]')

    return parser


def _require(parser, value, option, model):
    """Raise a concise argparse error when a model-specific option is missing."""
    if value is None:
        parser.error(f'{model} requires {option}')


def validate_args(parser, args):
    """Validate options whose requirement depends on the selected model."""
    if args.function == 'fasta':
        if args.model in ('extractk', 'extractc', 'count'):
            _require(parser, args.keywords, '-k/--keywords', f'fasta {args.model}')
        elif args.model == 'IDsimp':
            _require(parser, args.separator, '-s/--separator', 'fasta IDsimp')
        elif args.model == 'compare':
            _require(parser, args.newFASTA, '-x/--newFASTA', 'fasta compare')
        elif args.model in ('changeid', 'reverseq'):
            _require(parser, args.chrlist, '-c/--chrlist', f'fasta {args.model}')
    elif args.function == 'read':
        if args.model == 'extract':
            _require(parser, args.readid, '-k/--readid', 'read extract')
        elif args.model == 'filterONT':
            _require(parser, args.ONTquality, '-q/--ONTquality', 'read filterONT')
            _require(parser, args.ONTlength, '-l/--ONTlength', 'read filterONT')
    elif args.function == 'genome':
        if args.model == 'survey':
            _require(parser, args.path, '-p/--path', 'genome survey')
        elif args.model in ('telomere', 'split_by_gap'):
            _require(parser, args.inputGenome, '-i/--inputGenome', f'genome {args.model}')


def dispatch_fasta(args):
    from fasta_fun import (
        extract_based_coords,
        extract_based_keywords,
        extract_based_length,
        func_fasta_IDsimplified,
        func_fasta_compare,
        func_fasta_count,
        func_fasta_toPhy,
        func_get_fasta_length,
        genome_change_chr_name,
        genome_karyotype,
        genome_reverse_some_chr,
    )
    from util import split_fasta_based_bp, split_fasta_based_number

    if args.model == 'length':
        split_fasta_based_bp(args.input, args.length)
    elif args.model == 'number':
        split_fasta_based_number(args.input, args.number)
    elif args.model == 'extractk':
        extract_based_keywords(args.input, args.keywords, args.greedy)
    elif args.model == 'extractl':
        extract_based_length(args.input, args.length)
    elif args.model == 'extractc':
        extract_based_coords(args.input, args.keywords)
    elif args.model == 'falen':
        func_get_fasta_length(args.input)
    elif args.model == 'IDsimp':
        func_fasta_IDsimplified(args.input, args.separator)
    elif args.model == 'compare':
        func_fasta_compare(args.input, args.newFASTA)
    elif args.model == 'toPhy':
        func_fasta_toPhy(args.input)
    elif args.model == 'count':
        func_fasta_count(args.input, args.keywords)
    elif args.model == 'changeid':
        genome_change_chr_name(args.input, args.chrlist)
    elif args.model == 'reverseq':
        genome_reverse_some_chr(args.input, args.chrlist)
    elif args.model == 'karyotype':
        genome_karyotype(args.input)


def dispatch_args(args):
    """Run the selected subcommand after CLI validation has succeeded."""
    if args.function == 'fasta':
        dispatch_fasta(args)
    elif args.function == 'stat':
        from fasta_fun import func_stat
        func_stat(args.input, args.n50)
    elif args.function == 'gff':
        from gff_fun import func_get_ExonIntronInfo, func_get_LongTranscript, func_gff_NewGFF
        if not args.input.lower().endswith(('.gff', '.gff3')):
            raise ValueError('Error: Input file must have a .gff or .gff3 extension.')
        if args.model == 'extract':
            func_get_ExonIntronInfo(args.input, args.split)
        elif args.model == 'GetLongTrans':
            func_get_LongTranscript(args.input)
        elif args.model == 'NewGff':
            func_gff_NewGFF(args.input)
    elif args.function == 'gap':
        from gap_fun import func_get_gap_location
        if not args.input.lower().endswith(('.fasta', '.fa')):
            raise ValueError('Error: Input file must have a .fasta or .fa extension.')
        func_get_gap_location(args.input)
    elif args.function == 'read':
        from reads_fun import (
            func_get_reads_phreads,
            func_reads_ID_simplified,
            func_reads_deduplication_based_ID,
            func_reads_extract_based_ID,
            func_reads_filter_ONT_reads,
            func_reads_stat_from_seqkit_stat,
        )
        print('Input FASTQ reads: ', args.input)
        print('model: ', args.model)
        if args.model == 'phreads':
            print('Parse the quality value of reads.')
            func_get_reads_phreads(args.input)
        elif args.model == 'deduplication':
            print('Deduplication based on reads ID, if reads ID not uniq')
            func_reads_deduplication_based_ID(args.input)
        elif args.model == 'extract':
            print('Extract reads based on read ID')
            func_reads_extract_based_ID(args.input, args.readid)
        elif args.model == 'IDsimp':
            print('Simplified reads ID')
            func_reads_ID_simplified(args.input)
        elif args.model == 'filterONT':
            print('Filter ONT reads')
            func_reads_filter_ONT_reads(args.input, args.ONTquality, args.ONTlength, args.chopper)
        elif args.model == 'stat':
            print('Statistical reads based on seqkit stat')
            func_reads_stat_from_seqkit_stat(args.input)
    elif args.function == 'telomere':
        from telomere_fun import func_telomere_info
        func_telomere_info(args.input, args.length, args.minrptnum, args.output, args.telomere_script)
    elif args.function == 'genome':
        from genome_fun import genome_split_chr_by_gap, genome_survey, genome_telomere
        if args.model == 'survey':
            genome_survey(
                args.kmersize,
                args.length,
                args.count,
                args.threads,
                args.path,
                args.poly,
                args.genomescope1_script,
                args.genomescope2_script,
            )
        elif args.model == 'telomere':
            genome_telomere(args.inputGenome, args.telomere_script)
        elif args.model == 'split_by_gap':
            genome_split_chr_by_gap(args.inputGenome, args.gap_size)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.function is None:
        parser.print_help()
        return 0
    validate_args(parser, args)
    dispatch_args(args)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
