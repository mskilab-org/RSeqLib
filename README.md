
# RSeqLib
R interface to SeqLib C++ package 

## Table of contents
* [Installation](#installation)
* [Examples](#examples)

## Installation


```bash 
# Clone repository in a directory 
cd ~/temp/git/
git clone --recursive git@github.com:mskilab/RSeqLib.git
```



```R
## install package via devtools::load_all(). Future plans: get install_local() and install_github() to work. 
library(devtools)
setwd('~/temp/git/RSeqLib')
install_local('~/temp/git/RSeqLib')

```

    Installing RSeqLib
    '/nfs/sw/R/R-3.3.0/lib64/R/bin/R' --no-site-file --no-environ --no-save  \
      --no-restore --quiet CMD INSTALL '/tmp/RtmpNtm7wg/file39e07396510f/RSeqLib'  \
      --library='/gpfs/commons/groups/imielinski_lab/lib/R-3.3.0' --install-tests 
    
    Reloading installed RSeqLib
    Warning message:
    “replacing previous import ‘data.table::melt’ by ‘reshape2::melt’ when loading ‘RSeqLib’”Warning message:
    “replacing previous import ‘data.table::dcast’ by ‘reshape2::dcast’ when loading ‘RSeqLib’”

Examples 
--------




```R
## via functions available in package.
library(RSeqLib)
ls.str('package:RSeqLib')
```


    BWA : function (fasta = NULL, seq = NULL, seqname = "myseq", mc.cores = 1, hardclip = FALSE, 
        keep_sec_with_frac_of_primary_score = 0.9, max_secondary = 10)  
    initialize : Formal class 'nonstandardGenericFunction' [package "methods"] with 8 slots
    query : Formal class 'standardGeneric' [package "methods"] with 8 slots

```R
## create a bwa object by using the BWA() function and supply a random reference fasta to it. 
bwa <- BWA(seq = "CACTAGCTAGCTACGCGGGGGCGCGCGCGCGCGAAAAACACTTTCACAG")
## align two sequences with the reference fasta. 
## notice the result is in a GRanges format and there is a cigar field. 
query(bwa, c("CACTAGCTAGCTACGCGGGGGCGCG", "CACTAGCTAGCTACGCGCGAAAAACACTTTCACAG"))
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
## now, the best part! let's query to hg19 (replace the path with your own)
bwa <- BWA(fasta = "human_g1k_v37_decoy.fasta") # takes a few minutes.
query(bwa, c("TGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCAT"))
```

```
7                                                                                                        
GRanges object with 7 ranges and 12 metadata columns:      
    seqnames                 ranges strand |       qname        flag                                                              <Rle>              <IRanges>  <Rle> | <character> <character>                                                         
[1]        1 [    12060,     12119]      + |     myquery           0                                                         
[2]       23 [155257601, 155257660]      - |     myquery         272                                                         
[3]       15 [102519049, 102519108]      - |     myquery         272                                                         
[4]       16 [    61741,     61800]      + |     myquery         256                                                         
[5]        9 [    12173,     12232]      + |     myquery         256                                                         
[6]        2 [114358895, 114358954]      - |     myquery         272                                                         
[7]       12 [    93488,     93547]      - |     myquery         272                                                         
    mapq       cigar       rnext       pnext        tlen                                                                 
    <character> <character> <character> <character> <character>                                                             
[1]           0         60M           0          -1           0                                                             
[2]           0         60M           0          -1           0                                                             
[3]           0         60M           0          -1           0                                                             
[4]           0         60M           0          -1           0                                                             
[5]           0         60M           0          -1           0                                                             
[6]           0         60M           0          -1           0                                                             
[7]           0         60M           0          -1           0                                                             
    seq        qual                                                                                                      
    <character> <character>                                                 
[1] TGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCAT           *                                                 
[2] ATGAACAATGGCCCATCCCACTCTAGGCATGGCTCCTCTCCACAGGAAAACTCCACTCCA           *                                                 
[3] ATGAACAATGGCCCATCCCACTCTAGGCATGGCTCCTCTCCACAGGAAAACTCCACTCCA           *                                                 
[4] TGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCAT           *                                                 
[5] TGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCAT           *                                                 
[6] ATGAACAATGGCCCATCCCACTCTAGGCATGGCTCCTCTCCACAGGAAAACTCCACTCCA           *                                                 
[7] ATGAACAATGGCCCATCCCACTCTAGGCATGGCTCCTCTCCACAGGAAAACTCCACTCCA           *                                                 
    AS        DD    qwidth                                                                                               
    <integer> <integer> <integer>                                                                                           
[1]        60         0        60                                                                                           
[2]        60         0        60                                                                                           
[3]        60         0        60                                                                                           
[4]        60         0        60                                                                                           
[5]        60         0        60                                                                                           
[6]        60         0        60                                                                                           
[7]        55         0        60                                                                                           
-------                                                                                                                     
seqinfo: 7 sequences from an unspecified genome; no seqlengths         
```
