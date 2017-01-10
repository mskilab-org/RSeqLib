
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
## query two sequences with the reference fasta. 
## notice the result is in a GRanges format and there is a cigar field. 
query(bwa, c("CACTAGCTAGCTACGCGGGGGCGCG", "CACTAGCTAGCTACGCGCGAAAAACACTTTCACAG"))
```




```R

```
