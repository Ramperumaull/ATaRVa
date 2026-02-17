# Advanced Commands and Output format

## Table of contents

* [Advanced options](#advanced-options)
* [Examples](#examples)
* [Output format](#output-format)

## Advanced options
Optional arguments:
```
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
                        phasing. [default: 20]
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
  --methviz             write the methylation encoded sequence to the vcf for visualization
                        purpose. [default: False]
  --amplicon            genotype mode for targeted-sequenced samples.
                        In this mode, the default values for `max-reads` and `flank` values are 1000 and 20 respectively [default: False]
  --read-wise           Read-wise genotyping mode for BED file with dense regions. [default: False]
  --loci-wise           Loci-wise genotyping mode instead of Read-wise for BED file with sparse regions. [default: False]
  -log, --debug_mode    write the debug messages to log file. [default: False]
  -v, --version         show program's version number and exit
```
The details of each option are given below:

#### `--format`
**Expects**: *STRING*<br>
**Default**: *bam*<br>
This option sets the format of the alignment file. The default format is BAM. Specify the input format as `sam` for SAM files, `cram` for CRAM files, or `bam` for BAM files.  

#### `-q or --map-qual`
**Expects**: *INTEGER*<br>
**Default**: *5*<br>
Minimum mapping quality for the reads to be considered. All reads with a mapping quality below the specified value will be excluded during genotyping.

#### `--contigs`
**Expects**: *STRING*<br>
**Default**: *None*<br>
Specify the chromosome(s) for genotyping; repeat loci on all other chromosomes will be skipped. If no chromosomes are mentioned, repeats on all chromosomes in the BED file will be genotyped. eg: `--contigs chr1 chr12 chr22` will genotype only the repeat loci in these mentioned chromosomes in the BED file.

#### `--min-reads`
**Expects**: *INTEGER*<br>
**Default**: *10*<br>
Minimum number of the supporting reads required to genotype a locus. If the number of reads is less than this value, the locus will be skipped.

#### `--max-reads`
**Expects**: *INTEGER*<br>
**Default**: *100*<br>
Maximum number of supporting reads allowed for a locus to be genotyped. If the number of reads exceeds this limit, only this specified number of reads will be used for genotyping the locus.

#### `--snp-dist`
**Expects**: *INTEGER*<br>
**Default**: *3000*<br>
Maximum base pair (bp) distance from the flanks of the repeat locus to fetch SNPs from each read considered for phasing.

#### `--snp-count`
**Expects**: *INTEGER*<br>
**Default**: *3*<br>
Maximum number of SNPs to be used for read clustering and phasing.

#### `--snp-qual`
**Expects**: *INTEGER*<br>
**Default**: *20*<br>
Minimum Q value of the SNPs to be used for phasing.

#### `--flank`
**Expects**: *INTEGER*<br>
**Default**: *10*<br>
The number of base pairs in the flanking regions to be used for realignment.

#### `--snp-read`
**Expects**: *FLOAT*<br>
**Default**: *0.2*<br>
Minimum fraction of SNPs in the supporting reads of the repeat locus allowed for phasing.

#### `--meth-prob`
**Expects**: *FLOAT*<br>
**Default**: *0.5*<br>
Minimum probability value of methylation call to be considered for calculation, in the supporting reads of the repeat locus.

#### `--phasing-read`
**Expects**: *FLOAT*<br>
**Default**: *0.4*<br>
Minimum fraction of reads required in both clusters relative to the total supporting reads for the repeat locus after phasing.

#### `--karyotype`
**Expects**: *STRING*<br>
**Default**: *XX*<br>
Karyotype of the samples eg. XX or XY.

#### `-t or --threads`
**Expects**: *INTEGER*<br>
**Default**: *1*<br>
Number of threads to use for the process.

#### `--haplotag`
**Expects**: *STRING*<br>
**Default**: *None*<br>
Specify the haplotype tag to utilize phased information for genotyping. eg `HP`

#### `--decompose`
Performs motif-decomposition on ALT sequences.<br>
**NOTE: Only applicable for motif length <= 10**

#### `--methviz`
Calculates site-level methylation levels for CG bases(5mC) within the repeat region and writes them as a base64-encoded MV tag in the VCF SAMPLE column.

#### `--amplicon`
Genotyping mode for targeted sequencing data. In this mode, the default values for `max-reads` and `flank` values are 1000 and 20 respectively.

#### `--read-wise`
Classical ATaRVa genotyping mode, where loci are genotyped read-wise, utilizing the length advantage of the long reads to genotype multiple loci simultaneously (default : True)

#### `--loci-wise`
Genotyping mode for BED files with sparse regions across chromosomes or for BED files containing a small number of loci(<500). In this mode, loci are genotyped independently rather than using a read-wise approach.

#### `-v or --version`
Prints the version info of ATaRVa.

## Examples
The following examples assume the input reference genome is in `FASTA` format and is named ref.fa, the alignment file is in `BAM` format and is named input.bam, and the TR regions file is in `BED` format and is named regions.bed.gz.

### Basic usage
To run ATaRVa with default parameters, use the following command:
```bash
$ atarva -f ref.fa --bam input.bam -r regions.bed.gz
```
### With karyotype
To run ATaRVa with sex chromosome karyotype, use the following command:
```bash
$ atarva -f ref.fa --bam input.bam -r regions.bed.gz --karyotype XY
```
With multiple bams:
```bash
$ atarva -f ref.fa --bam input1.bam input2.bam -r regions.bed.gz --karyotype XY XX
```
### With haplotag
To run ATaRVa on haplotagged alignment file, use the following command:
```bash
$ atarva -f ref.fa --bam input.bam -r regions.bed.gz --haplotag HP
```
### With amplicon
To run ATaRVa on targeted sequencing file, use the following command:
```bash
$ atarva -f ref.fa --bam input.bam -r regions.bed.gz --amplicon
```
### Stringent parameter usage
To run ATaRVa with stringent parameters, use the following command:
```bash
$ atarva -q 30 --snp-count 5 --snp-qual 25 --min-reads 20 -t 32 -fi ref.fa --bam input.bam -r regions.bed.gz
# The above command with --snp-count 5 will use a maximum of five heterozygous SNPs to provide accurate genotypes, but only when phasing is based on SNPs and not on length.
```
### Genotyping TRs from specific chromosome/s
To genotype TRs from specific chromosomes only, run ATaRVa with the following command:
```bash
$ atarva --contigs chr9 chr15 chr17 chrX -t 32 -f ref.fa --bam input.bam -r regions.bed.gz
```
### For input alignment file other than `bam`
```bash
# input cram file
$ atarva --format cram -f ref.fa --bam input.cram -r regions.bed.gz

# input sam file
$ atarva --format sam -f ref.fa --bam input.sam -r regions.bed.gz
```
### Usage in docker
To run ATaRVa in docker container, use the following command:
```bash
$ docker run -i -t --rm -v /path_of_necessary_files/:/folder_name atarva:latest -f /folder_name/ref.fa --bam /folder_name/input.bam -r /folder_name/regions.bed.gz
``` 

In all the above examples, the output of ATaRVa is saved to input.vcf unless -o is specified.

## Output format

### `-o or --vcf`
**Expects**: *STRING (to be used as filename)*<br>
**Default**: *Input Alignment Filename + .vcf*<br>
If this option is not provided, the default output filename will be the same as the input alignment filename, with its extension replaced with '.vcf'. For example, if the input filename is `input.bam`, the default output filename will be `input.vcf`. If the input filename does not have any extension, .vcf will be appended to the filename.
Each entry includes the fields specified in the [Variant Calling Format (VCF)](https://samtools.github.io/hts-specs/VCFv4.3.pdf), as described in the table below.
|     FIELD     |           DESCRIPTION           | 
|---------------|---------------------------------|
| CHROM | Chromosome that contains the repeat region |
| POS | Start position of the repeat region |
| ID |  Region identifier (set to '.') |
| REF | Reference sequence of the repeat region |
| ALT | Sequence of the repeat alleles in the sample |
| QUAL | Quality score of the genotype (set to '0') |
| FILTER | Filter status (PASS, LESS_READS) |
| INFO | Information about the TR region |
| FORMAT | Data type of the genotype information |
| SAMPLE | Values of the genotype information for the TR region |

#### INFO fields
The `INFO` field describes the general structure of the repeat region and includes the following details:
|     INFO      |           DESCRIPTION           | 
|---------------|---------------------------------|
| AC | Total number of respective ALT alleles in called genotypes |
| AN | Total number of alleles in called genotypes |
| MOTIF | Motif of the repeat region |
| START | Start position of the repeat region |
| END | End position of the repeat region |
| ID  | Tag fetched form the extra column in BED file |
| REFCN | Reference allele copy number |

#### FORMAT fields
The `FORMAT` fields and their values are provided in the last two columns of the VCF file, containing information about each genotype call. These columns include the following fields:  
|     FORMAT      |           DESCRIPTION           | 
|-----------------|---------------------------------|
| GT | Genotype of the sample |
| AL | Length of the alleles in base pairs |
| CN | Motif copy number for each allele |
| AR | Central 95% range of allele lengths in each cluster |
| SD | Number of supporting reads for each alleles |
| DP | Number of the supporting reads for the repeat locus |
| SN | Number of SNPs used for phasing |
| SQ | Phred-scale qualities of the SNPs used for phasing |  
| MA | Mean methylation level for each allele | 
| MR | Number of reads supporting methylation info for each allele | 
| DS | Motif decomposed sequence of the alternate alleles |
| MV | Visual methylation encodings for the alleles |

**NOTE: Loci missing in the VCF either have no reads mapped to them, contain reads that do not fully enclose the repeat region, or have reads with low mapping quality (mapQ).**