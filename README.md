[![Build Status](https://travis-ci.org/mskilab/RSeqLib.svg?branch=master)](https://travis-ci.org/mskilab/RSeqLib)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/RSeqLib.svg)](https://codecov.io/github/mskilab/RSeqLib?branch=master)



# RSeqLib
R interface to SeqLib C++ package 

## Table of contents
* [Installation](#installation)
* [Examples](#examples)

## Installation

```R
devtools::install_github('mskilab/RSeqLib')
```

## Examples 

### BWA on character vector reference

```R

## instantiate the BWA object around a character vector custom reference sequence 
bwa = BWA(seq = c("CACTAGCTAGCTACGCGGGGGCGCGCGCGCGCGAAAAACACTTTCACAG"))

## align two sequences to this reference using "[" operator
## returned result is a granges on the above sequence
bwa[c("CACTAGCTAGCTACGCGGGGGCGCG", "CACTAGCTAGCTACGCGCGAAAAACACTTTCACAG")]
```

```
GRanges object with 2 ranges and 12 metadata columns:                                                                         
      seqnames    ranges strand |       qname        flag        mapq                                              
         <Rle> <IRanges>  <Rle> | <character> <character> <character>
  [1]        1  [ 0, 24]      + |     myquery           0          36
  [2]        1  [27, 48]      + |     myquery           0          18
            cigar       rnext       pnext        tlen
      <character> <character> <character> <character>
  [1]         25M           0          -1           0
  [2]      13S22M           0          -1           0
                                      seq        qual        AS        DD
                              <character> <character> <integer> <integer>
  [1]           CACTAGCTAGCTACGCGGGGGCGCG           *        25         0
  [2] CACTAGCTAGCTACGCGCGAAAAACACTTTCACAG           *        22         0
         qwidth
      <integer>
  [1]        25
  [2]        35
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```R

## you can make a reference genome with a multiple named contigs
## in this case we make them identical 
bwa = BWA(seq = c(
allele1 = "CACTAGCTAGCTACGCGGGGGCGCGCGCGCGCGAAAAACACTTTCACAG",
allele2 = "CACTAGCTAGCTACGCGGGGGCGCGCGCGCGCGAAAAACACTTTCACAG"))

## notice that now the seqnames of alignments are "allele1" and "allele2"
## and we get multiple alignments per query, each with mapq0, as expected
bwa[c("CACTAGCTAGCTACGCGGGGGCGCG", "CACTAGCTAGCTACGCGCGAAAAACACTTTCACAG")]
```

```
GRanges object with 4 ranges and 12 metadata columns:
      seqnames    ranges strand |       qname        flag        mapq
         <Rle> <IRanges>  <Rle> | <character> <character> <character>
  [1]  allele1      0-24      + |     myquery           0           0
  [2]  allele2      0-24      + |     myquery         256           0
  [3]  allele2     27-48      + |     myquery           0           0
  [4]  allele1     27-48      + |     myquery         256           0
            cigar       rnext       pnext        tlen
      <character> <character> <character> <character>
  [1]         25M           0          -1           0
  [2]         25M           0          -1           0
  [3]      13S22M           0          -1           0
  [4]      13S22M           0          -1           0
                                      seq        qual        AS        DD
                              <character> <character> <integer> <integer>
  [1]           CACTAGCTAGCTACGCGGGGGCGCG           *        25         0
  [2]           CACTAGCTAGCTACGCGGGGGCGCG           *        25         0
  [3] CACTAGCTAGCTACGCGCGAAAAACACTTTCACAG           *        22         0
  [4] CACTAGCTAGCTACGCGCGAAAAACACTTTCACAG           *        22         0
         qwidth
      <integer>
  [1]        25
  [2]        25
  [3]        35
  [4]        35
  -------
  seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

### BWA on pre-existing fasta index 

```R
## we can also use a pre existing fasta index
## this mini-genome was built using a chunk of hg19 chromosome 2
bwa = BWA(system.file('ext/mini.genome.fasta', package = 'RSeqLib'))

## inspect the bwa object
bwa

```
```
RSeqLib BWA object with params mc.cores = 1, hardclip = 0, keep_sec_with_frac_of_primary_score = 0.9, max_secondary = 10
```

```R
## can retrieve the fasta path or character vector
## corresponding of the reference of the  bwa object using genome(bwa)
genome(bwa)
```

```
[1] "/gpfs/commons/groups/imielinski_lab/git/mskilab/rSeqLib/inst/ext/mini.genome.fasta"

```

``` R

## now let's query a character vector against this fasta reference
qstring = c("TGCTCTGGCACAAAGATAGGCATGCTCAGCCCACCTCTGCCAGCACTGGGGCTGAAGACTGGCCCACGTGTCCTCTCAATCCCCAGGACAACTTCACTATGGCTTTCATTAATAA",
             "GCACAGATATCAACATAAGGACACAGAAAGCAGGAAAAGCAAGGAAATATGACTCCTTCAAAGGAACACAATAATTTTTCAGCATTAGATCTTAATCAGAAAGAACTTACTCCCAGATAAATAATTCAAAATAA")

## again we use the "[" operator to retrieve results
bwa[qstring]
```

```
GRanges object with 2 ranges and 12 metadata columns:
      seqnames        ranges strand |       qname        flag        mapq
         <Rle>     <IRanges>  <Rle> | <character> <character> <character>
  [1]        1 108294-108408      - |     myquery          16          60
  [2]        1 107951-108090      - |     myquery          16          60
            cigar       rnext       pnext        tlen
      <character> <character> <character> <character>
  [1]        115M           0          -1           0
  [2]   24M6D110M           0          -1           0
                                                                                                                                         seq
                                                                                                                                 <character>
  [1]                    TTATTAATGAAAGCCATAGTGAAGTTGTCCTGGGGATTGAGAGGACACGTGGGCCAGTCTTCAGCCCCAGTGCTGGCAGAGGTGGGCTGAGCATGCCTATCTTTGTGCCAGAGCA
  [2] TTATTTTGAATTATTTATCTGGGAGTAAGTTCTTTCTGATTAAGATCTAATGCTGAAAAATTATTGTGTTCCTTTGAAGGAGTCATATTTCCTTGCTTTTCCTGCTTTCTGTGTCCTTATGTTGATATCTGTGC
             qual        AS        DD    qwidth
      <character> <integer> <integer> <integer>
  [1]           *       115         0       115
  [2]           *       122         0       134
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths

```


### Fermi on character vectors

```R

## character vector of "reads"
reads =
c('CTGCTAAGGAGATGGTGTGGCATCTGTGGAAGTCTGCTTGCACACAGCTGGAAGCCTTCAGGGGAAATAACAAACTGCTTCTCTTGCTGGGATTCCGACAT',
'CTAAGGAGATGGTGTGGCATCTGTGGAAGTCTGCTTGCACACAGCTGGAAGCCTTCAGGGGAAATAACAAACTGCTTCTCTTGCTGGGATTCCGACATGTC',
'TAAGGAGATGGTGTGGCATCTGTGGAAGTCTGCTTGCACACAGCTGGAAGCCTTCAGGGGAAATAACAAACTGCTTCTCTTGCTGGGATTCCGACATGTCC',
'GATGGTGTGGCATCTGTGGAAGTCTGCTTGCACACAGCTGGAAGCCTTCAGGGGAAATAACAAACTGCTTCTCTTGCTGGGATTCCGACATGTCCAAATAT',
'GGTGTGGCATCTGTGGAAGTCTGCTTGCACACAGCTGGAAGCCTTCAGGGGAAATAACAAACTGCTTCTCTTGCTGGGATTCCGACATGTCCAAATATGTC')
fermi = Fermi(reads, assemble = TRUE)
contigs(fermi)
```
```
[1] "GACATATTTGGACATGTCGGAATCCCAGCAAGAGAAGCAGTTTGTTATTTCCCCTGAAGGCTTCCAGCTGTGTGCAAGCAGACTTCCACAGATGCCACACCATCTCCTTAGCAG"
```

### Fermi on objects coercible to data.tables

```R
## can also assemble a GRanges of reads (or any object
## coercible to a data.frame, or data.table)
## input just needs $seq, (optional) $qual, and (optional) $qname fields 
## (so fermi takes base qualities into account)
reads = readRDS(system.file('ext/reads.rds', package = 'RSeqLib'))

names(as.data.table(reads))
```
```
 [1] "seqnames" "start"    "end"      "width"    "strand"   "qname"
 [7] "flag"     "qwidth"   "mapq"     "cigar"    "mrnm"     "mpos"
[13] "isize"    "seq"      "qual"     "MD"       "MQ" 
```

```R
## make fermi object and assemble with error correction
fermi = Fermi(reads, assemble = TRUE)
fermi
```

```
RSeqLib Fermi object with 6365 reads and 43 contigs:

```

```R

## retrieve the contigs and align back to the reference
gr = bwa[contigs(fermi)]

## this result is a granges aligned to the fasta corresponding to the bwa object above
gr[1]

```
```
GRanges object with 1 range and 12 metadata columns:
      seqnames        ranges strand |       qname        flag        mapq
         <Rle>     <IRanges>  <Rle> | <character> <character> <character>
  [1]        1 107677-108408      - |     myquery          16          60
            cigar       rnext       pnext        tlen
      <character> <character> <character> <character>
  [1]        732M           0          -1           0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               seq
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       <character>
  [1] TTTTTTTTTTAATTGTCTGATTGGAGTATTTAAAAAGATCTGTCTTCAAGTTCTGAGATTCTTTATTCTACCTGATCTAATCTATTGTTGATGCTTTCAAATGTGTTTTTTTTATTTCCTTCAATGAATTCTTCAGTTCCAGAATTTCTATTTGGTTCATTTAAAAATGTCTATTTCTGTGTAAATTTCTCATTTATATTCTGAATTGATTTTCTAATTTCTTTCTATTGTTTTTCAGAATTTTCTTGTATTTTACTGAGTTTCTTTAAAATCATTATTTTGAATTATTTATCTGGGATGTTATGTAAGTTCTTTCTGATTAAGATCTAATGCTGAAAAATTATTGTGTTCCTTTGAAGGAGTCATATTTCCTTGCTTTTCCTGCTTTCTGTGTCCTTATGTTGATATCTGTGCATCTGGCATAATAGTCACCTTCTATTTTTGAATTTGCTTTTATAGGGGAGGACTTTTTCCTGAAGATGTATCTCTGGTATCAGTTGGGTAGAGCACTTTGGCTTTGATTCTGGGTGCATCCAGTAGTGTAGTCTCCATATGATTTCTTTGGCTGTAAACAGTGTTAGTAGCATCTGTGATTTCCTTCATGGATTAGAGTGTGGTTATTAATGAAAGCCATAGTGAAGTTGTCCTGGGGATTGAGAGGACACGTGGGCCAGTCTTCAGCCCCAGTGCTGGCAGAGGTGGGCTGAGCATGCCTATCTTTGTGCCAGAGCA
             qual        AS        DD    qwidth
      <character> <integer> <integer> <integer>
  [1]           *       732         0       732

```

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill Cornell Medicine
> Core Member, New York Genome Center

> Jeremiah Wala - MD PhD Student, Harvard / MIT / HST

> Khagay Nagdimov - Intern, Imielinski Lab, New York Genome Center
