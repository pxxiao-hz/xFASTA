# xFASTA
## Usage
```shell
$ xFASTA -h
usage: xFASTA [-h]  ...

pxx 的分析小工具

positional arguments:

    fasta     处理 FASTA 文件
    stat      统计基因组的相关信息
    gff       Perform operations on GFF file
    gap       Genome gap statistics
    read      Reads information
    telomere  Telomere information
    genome    Genome operation

optional arguments:
  -h, --help  show this help message and exit
```
You can view the help documentation for each subcommand by running:
```shell
xFASTA <subcommand> -h
```
## Release notes
### v1.0.2
* Add some features, such as telomere identification, genome survey, etc.
### v1.0.1
* Add extract read sequence based read ID.
### v1.0.0
* The first version.

