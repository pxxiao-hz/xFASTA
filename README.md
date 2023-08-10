# xFASTA
```shell
usage: xFASTA [-h]  ...

pxx 的分析小工具

positional arguments:

    fasta     处理 FASTA 文件
    stat      统计基因组的相关信息
    gff       Perform operations on GFF file
    gap       Genome gap statistics
    read      Reads information

optional arguments:
  -h, --help  show this help message and exit
```
## Function1: fasta
```
usage: xFASTA fasta [-h] -i INPUT [-m] [-l LENGTH] [-n NUMBER] [-k KEYWORDS] [-g GREEDY]
                    [-s SEPARATOR] [-x NEWFASTA] [-p TOPHY]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        file input: .fasta or .fa
  -m , --model          Process on .fasta or .fa
                        length: Split File Based on Length
                        number: Split File Based on Number
                        extractk: Sequence Extraction Based on Keywords,
                        	 -k required
                        extractl: Sequence Extraction Based on seq length,
                        	 -l required
                        extractc: Sequence Extraction Based on coords,
                        	 -k required, e.g.: chr1-10-20; chr1-10-, means 10 - end of chr1
                        	 Note: 1-based in this script
                        falen: Length of each sequence
                        IDsimp: Simplified ID Name,
                        	 -s required
                        compare: Compare the contents of two FASTA files for equality
                        toPhy: Convert FASTA format to phy format
                        count: Count the number of strings
  -l LENGTH, --length LENGTH
                        拆分成多少bp的小文件 or 根据序列的长度从基因组文件提取序列 [20000000]
  -n NUMBER, --number NUMBER
                        拆分成多个文件，每个文件包含number个序列 [50]
  -k KEYWORDS, --keywords KEYWORDS
                        Keywords in the Sequence to be Extracted or Count.
                        For Extracted, Chr1 or Chr1,Chr2,Chr3
  -g GREEDY, --greedy GREEDY
                        Greedy match or not, True or False
  -s SEPARATOR, --separator SEPARATOR
                        Separator of sequence ID; "None" == "\TAB"
  -x NEWFASTA, --newFASTA NEWFASTA
                        File 2 of two FASTA files for compare
  -p TOPHY, --toPhy TOPHY
                        FASTA file to Phy
```
## Function2: stat
```
usage: xFASTA stat [-h] -i INPUT [-n N50]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        输入文件，待处理的 FASTA 文件
  -n N50, --n50 N50     基因组组装结果的统计
```
## Function3: gff
```
usage: xFASTA gff [-h] -i INPUT -m  [-s SPLIT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        file input: GFF file
  -m , --model          GFF sub command, including sort, extract, GetLongTrans
                        sort: based on chr and start
                        extract: extract exon, intron information, output format: bed
                        GetLongTrans: only retain longest information, output format: gff
                        NewGff: if gff not have "gene" mark, please use NewGff
  -s SPLIT, --split SPLIT
                        the 9th column of gff file, maybe startswith "ID=", or "Parent=
```
## Function4: gap
```
usage: xFASTA gap [-h] -i INPUT -m

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        finle input: Genome fasta file
  -m , --model          gap sub coomand, including coord
                        coord: Obtaining the location information of gaps, NNNN... or nnnn... was defined as gap
```
## Function5: read
```
usage: xFASTA read [-h] -i INPUT -m

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        The input FASTQ file, can be in .fastq or .fastq.gz format
  -m , --model          read sub command, including phreads
                        phreads: Parse the quality value of reads
                        deduplication: remove reads based on reads ID, if reads ID not uniq
```
