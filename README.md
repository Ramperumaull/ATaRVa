# ATaRVa - a tandem repeat genotyper
![Badge-PyPI](https://img.shields.io/badge/PyPI-v0.5.0-brightgreen)
![Badge-License](https://img.shields.io/badge/License-MIT-blue)

<p align=center>
  <img src="lib/atrv_logo.png" alt="Logo of ATaRVa" width="200"/>
</p>

ATaRVa (pronounced uh-thur-va, IPA: /əθərvə/, Sanskrit: अथर्व) is a technology-agnostic tandem repeat genotyper, specially designed for long read data. The name expands to **A**nalysis of **Ta**ndem **R**epeat **Va**riation, and is derived from the the Sanskrit word _Atharva_ meaning knowledge.


## Motivation
Long-read sequencing propelled comprehensive analysis of tandem repeats (TRs) in genomes. Current long-read TR genotypers are either platform specific or computationally inefficient. ATaRva outperforms existing tools while running an order of magnitude faster. ATaRVa also supports multi-threading, haplotyping, and motif decomposition, making it an invaluable tool for population scale TR analyses.

## Table of contents:

* [Installation](#installation)
  * [PyPI installation](#pypi-installation)
  * [Docker installation](#docker-installation)
* [Usage](#usage)
  * [`genotype` command](#genotype-command)
    * [Reference genome](#reference-genome)
    * [Alignment file](#alignment-file)
    * [Region file](#region-file)
  * [`merge` command](#merge-command)
* [Changelog](#changelog)
* [Analysis script](#analysis-script)
* [Citation](#citation)
* [Contact](#contact)

## Installation

### PyPI installation

ATaRVa can be directly installed using pip with the package name `ATaRVa`.

```bash
$ pip install ATaRVa
```
Alternatively, it can be installed from the source code:<br>
It is recommended to install this inside a Python virtual environment.

```bash
# Create a python env
$ python -m venv atarva_env

# Activate the env
$ source atarva_env/bin/activate
$ pip install build

# Download the git repo
$ git clone https://github.com/SowpatiLab/ATaRVa.git

# Install
$ cd ATaRVa
$ python -m build
$ pip install .

# Deactivate the env
$ deactivate
```
Both of the methods add a console command `atarva`, which can be executed from any directory

<!-- **NOTE: This tool has been tested and is recommended to be used with Python versions between 3.9 and 3.12 (inclusive).** -->

### Docker installation
ATaRVa can also be installed using the provided **Docker** image with the following steps:
```bash
$ cd ATaRVa
$ docker build --network host -t atarva
```
## Usage
The help message and available subcommands can be accessed using

```bash
$ atarva genotype -h
#  or
$ atarva genotype --help
```
which gives the following output

```
ATaRVa - Analysis of Tandem Repeat Variants
Sowpati Lab

Usage:
    atarva [OPTIONS] <COMMAND>

Commands:
  genotype  Tandem Repeat Genotyper
  merge     Merge ATaRVa VCF files

Options:
  -h, --help     Print help
  -v, --version  Print version
```

## `genotype` command
`atarva genotype` accepts read alignments and a set of TR regions of interest and outputs TR genotypes, including the consensus sequence, allele length, and decomposed motifs.

Overview of the ATaRVa worflow:
1. ATaRVa processes the input BAM file read-wise, assuming that most long reads span multiple TR loci.
2. After flank realignment and adjustment of read-wise allele lengths, ATaRVa clusters reads into haplotypes using nearby informative *SNV*s, or applies a *k-means* clustering approach when SNV information is unavailable.
3. It derives consensus allele sequences using partial order alignment, decomposes each TR allele into motif-level representations, and outputs the results in VCF format.

The help message and available options can be accessed using

```bash
$ atarva genotype -h
#  or
$ atarva genotype --help
```
which gives the following output

```
usage: atarva genotype  [-h] -f <FILE> -b <FILE> [<FILE> ...] -r <FILE> [--format <STR>] [-q <INT>]
                        [--contigs CONTIGS [CONTIGS ...]] [--min-reads <INT>] [--max-reads <INT>]
                        [--snp-dist <INT>] [--snp-count <INT>] [--snp-qual <INT>] [--flank <INT>]
                        [--snp-read <FLOAT>] [--meth-prob <FLOAT>] [--phasing-read <FLOAT>] [-o <FILE>]
                        [--karyotype KARYOTYPE [KARYOTYPE ...]] [-t <INT>] [--haplotag <STR>]
                        [--decompose] [--amplicon] [--read-wise] [--loci-wise] [-log] [-v]

Required arguments:
  -f <FILE>, --fasta <FILE>
                        input reference fasta file
  -b <FILE> [<FILE> ...], --bam <FILE> [<FILE> ...]
                        samples alignment files. allowed formats: SAM, BAM, CRAM
  -r <FILE>, --regions <FILE>
                        input regions file. the regions file should be strictly in bgzipped
                        tabix format. If the regions input file is in bed format. First sort it
                        using bedtools. Compress it using bgzip. Index the bgzipped file with
                        tabix command from samtools package.

Optional arguments:
  --format <STR>        format of input alignment file. allowed options: [cram, bam, sam]. default: [bam]
  -q <INT>, --map-qual <INT>
                        minimum mapping quality of the reads to be considered. [default: 5]
  --contigs CONTIGS [CONTIGS ...]
                        contigs to get genotyped [chr1 chr12 chr22 ..]. If not mentioned every
                        contigs in the region file will be genotyped.
  --min-reads <INT>     minimum read coverage after quality cutoff at a locus to be genotyped. [default: 10]
  --max-reads <INT>     maximum number of reads to be used for genotyping a locus. [default: 100]
  --snp-dist <INT>      maximum distance of the SNP from repeat region to be considered for
                        phasing. [default: 3000]
  --snp-count <INT>     number of SNPs to be considered for phasing (minimum value = 1).
                        [default: 3]
  --snp-qual <INT>      minimum basecall quality at the SNP position to be considered for
                        phasing. [default: 13]
  --flank <INT>         length of the flanking region (in base pairs) to search for insertion
                        with a repeat in it. [default: 10]
  --snp-read <FLOAT>    a positive float as the minimum fraction of snp's read contribution to
                        be used for phasing. [default: 0.25]
  --meth-prob <FLOAT>   a minimum probability cutoff for methylation. [default: 0.5]
  --phasing-read <FLOAT>
                        a positive float as the minimum fraction of total read contribution from
                        the phased read clusters. [default: 0.4]
  -o <FILE>, --vcf <FILE>
                        name of the output file, output is in vcf format. [default: sys.stdout]
  --karyotype KARYOTYPE [KARYOTYPE ...]
                        karyotype of the samples [XY XX]
  -t <INT>, --threads <INT>
                        number of threads. [default: 1]
  --haplotag <STR>      use haplotagged information for phasing. eg: [HP]. [default: None]
  --decompose           write the motif-decomposed sequence to the vcf. [default: False]
  --amplicon            genotype mode for targeted-sequenced samples.
                        In this mode, the default values for `max-reads` and `flank` values are 1000 and 20 respectively [default: False]
  --read-wise           Read-wise genotyping mode for BED file with dense regions. [default: False]
  --loci-wise           Loci-wise genotyping mode instead of Read-wise for BED file with sparse regions. [default: False]
  -log, --debug_mode    write the debug messages to log file. [default: False]
  -v, --version         show program's version number and exit
```

The details of required input files are given below:

### Reference genome
#### `-f or --fasta`
**Expects**: *FILE*<br>
**Default**: *None*<br>
The `-f` or `--fasta` option is used to specify the input FASTA file. The corresponding index file (`.fai`) should be in the same directory. ATaRVa uses [pysam](https://github.com/pysam-developers/pysam)'s `FastaFile` parser to read the input FASTA file.

### Alignment file
#### `-b or --bam`
**Expects**: *FILE*<br>
**Default**: *None*<br>
The `-b` or `--bam` option is used to specify one or more input alignment files in the same format. ATaRVa accepts any of the three alignment formats: SAM, BAM, or CRAM. The alignment file should be sorted by coordinates. The format should be specified using the `--format` option. The corresponding index file (`.bai` or `.csi`) should be located in the same directory. An alignment file can be sorted and indexed using the following commands:

```bash
# to sort the alignment file
$ samtools sort -o sorted_output.bam input.bam

# to generate .bai index file
$ samtools index -b sorted_output.bam
```

An alignment file containing at least one of the following tags is preferred for faster processing: `MD` tag, `CS` tag, or a `CIGAR` string with `=/X` operations.

- The CS tag is generated using the --cs option when aligning reads with the [minimap2](https://github.com/lh3/minimap2) aligner. (`--cs=short` is prefered over `--cs=long`)
- The MD tag can be generated using the --MD option in minimap2.

If the alignment files were generated without any of these tags, you can generate the `MD` tag by running the following command to 

```bash
# input: reference genome fasta file & alignment file
# output: an alignment file with MD tag in it

# for generating MD tag
$ samtools calmd -b aln.bam ref.fa > aln_md.bam
```
### Region file
#### `-r or --regions`
**Expects**: *FILE*<br>
**Default**: *None*<br>
The `-r` or `--regions` option is used to specify the input TR regions file. ATaRVa requires a sorted, bgzipped BED file of TR repeat regions, along with its corresponding tabix-indexed file. The BED file should contain the following columns:

1. Chromosome name where TR is located
2. Start position of the TR
3. End position of the TR
4. Repeat motif
5. Motif length

Below is an example of a repeat region BED file. **NOTE: The BED file should either have no header or a header that starts with `#` symbol. The .gz and .tbi files should be in same directory**

| #CHROM | START | END | MOTIF | MOTIF_LEN |
|--------|-------|-----|-------|-----------|
| chr1   | 10000 | 10467 | TAACCC | 6    |
| chr1   | 10481 | 10497 | GCCC | 4      |
| chr2   | 10005 | 10173 | CCCACACACCACA | 13 |
| chr2   | 10174 | 10604 | ACCCTA | 6    |
| chr17  | 60483 | 60491 | AGA    | 3    |

To sort, bgzip, and index the BED file, use the following commands:

#### Sort
```bash
# input: Unsorted bed file
# output: Sorted bed file

# Sorting the BED file using sort
$ sort -k1,1 -k2,2n input.bed > sorted_output.bed
# or using bedtools
$ bedtools sort -i input.bed > sorted_output.bed
```
#### Bgzip
```bash
# input: Sorted bed file
# output: bgzipped bed file

# To keep the original file unchanged and generate separate gz file
$ bgzip -c sorted_output.bed > sorted_output.bed.gz
# or to compress the original file; converts sorted_output.bed to sorted_output.bed.gz
$ bgzip sorted_output.bed
```
#### Index
```bash
# input: bgzipped bed file
# output: tabix indexed file (.tbi)

# install samtools to use tabix
$ tabix -p bed sorted_output.bed.gz
```
For detailed information on advanced genotyping options, refer to the [Advanced Commands and Usage](./docs/genotype_usage.md) documentation.

## `merge` command
`atarva merge` merges VCF files generated by ATaRVa. The tool is optimized to handle large datasets efficiently by reading and processing multiple files in small chunks, thereby avoiding excessive memory usage and ensuring fast, memory-efficient execution.

The tool requires the following inputs:
- BGZipped and tabix-indexed VCF files
- A BGZipped and tabix-indexed BED file specifying the regions of interest

The help message and available options can be accessed using

```bash
$ atarva merge -h
#  or
$ atarva merge --help
```
which gives the following output

```
usage: atarva merge [-h] -r <FILE> -i <FILE> [<FILE> ...] -f <FILE> [--contigs CONTIGS [CONTIGS ...]] [-o <STR>] [-p <INT>]

Required arguments:
  -r <FILE>, --regions <FILE>
                        input regions file. the regions file should be strictly in bgzipped tabix format. If the regions input file is in bed format. First sort it using bedtools. Compress it using
                        bgzip. Index the bgzipped file with tabix command from samtools package.
  -i <FILE> [<FILE> ...], --vcfs <FILE> [<FILE> ...]
                        text file containing paths to input vcf files to be merged. The text file should list each path on a separate line. The vcf files should be strictly in bgzipped tabix format. If
                        the vcfs input file is in vcf format. First sort it using bcftools. Compress it using bgzip. Index the bgzipped file with tabix command from samtools package.
  -f <FILE>, --fasta <FILE>
                        input reference fasta file. The file should be indexed.

Optional arguments:
  --contigs CONTIGS [CONTIGS ...]
                        contigs to get merged [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be merged.
  -o <STR>, --outname <STR>
                        name of the output file, output is in vcf format.
  -t <INT>, --threads <INT>
                        number of threads. [default: 1]
```
**NOTE: This will merge only those loci that are present in the input BED file** <br>
For detailed information on advanced merging options, refer to the [Tamatr](./docs/merge_usage.md) documentation.

## Changelog
### v0.5.0
* Changed the VCF-START column into 1-based coordinate system
* Included `START` tag in VCF-INFO column with 0-based coordinate system
* Added `MR` tags in the VCF_SAMPLE column to report the supporting read count for mean methylation level
* Added a confirmation step to check `MM` extraction from the reverse strand
* Added a `loci-wise` flag to perform region-wise genotyping (instead of the default read-wise mode) for BED files with sparse regions
* Improved Motif-decomposition script to maintain consistent representation of a motif (cyclic variation check)

### v0.4.0
* Added `MM` tag in VCF_SAMPLE column for mean methylation level
* Modified `AR` tag in VCF-SAMPLE column with central 95% allele range
* Implemented DBSCAN clustering in `amplicon` mode to check for multiple clusters
* Fixed bugs in decomposition function [#8](https://github.com/SowpatiLab/ATaRVa/issues/8)

### v0.3.1
* Added checkpoint in amplicon mode for non-repeatedness in ALT sequence
* Refined Motif-decomposition sequence for motif breaks
* Added `AR` tag in VCF-SAMPLE column for allele range

### v0.3.0
* Added `--amplicon` mode for targeted sequencing data
* Added function to convert eqx read sequence
* Improved Outlier cleaning in K-Means clustering
* Implemented De-novo motif identification in motif-decomposition
* Added optional tag `ID` in INFO field if BED input has additional column

### v0.2.0
* Added `--haplotag` argument to enable the use of haplotag information for genotyping.
* Fixed bugs in SNP-based clustering.
* Replaced the use of the mode function with a consensus-based approach for final allele derivation.
* Removed `PC` tag from the FORMAT field of the output VCF.

### v0.1.2
* Modified input arguments.

### v0.1.1
* Added a Mac OS compatible <code>.so</code> file.

### v0.1
* First release.

## Analysis script
All scripts used for analysis are provided in [ATaRVa_Manuscript](https://github.com/SowpatiLab/ATaRVa_Manuscript)

## Citation
If you find ATaRVa useful for your research, please cite it as follows:

ATaRVa: Analysis of Tandem Repeat Variation from Long Read Sequencing data  
_Abishek Kumar Sivakumar, Sriram Sudarsanam, Anukrati Sharma, Akshay Kumar Avvaru, Divya Tej Sowpati_ <br>
_BioRxiv_, **doi:** https://doi.org/10.1101/2025.05.13.653434

## Contact
For queries or suggestions, please contact:

Divya Tej Sowpati - tej at csirccmb dot org

Abishek Kumar S - abishekks at csirccmb dot org

Akshay Kumar Avvaru - avvaruakshay at gmail dot com
