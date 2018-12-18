#' @useDynLib RSeqLib
#' @importFrom Rcpp sourceCpp
#' @import data.table reshape2
NULL



#' @name BWA-claas
#' @title BWA-claas
#' @description
#' BWA class
#' @export
#' @author Marcin Imielinski
setClass( "BWA", representation( bwa = "externalptr",
                                params = "list"
                                ), contains = "externalptr" )

#' @name show
#' @title show
#' @description show
#' @author Marcin Imielinski
setMethod('show', 'BWA', function(object)
{
  message('RSeqLib BWA object with params ', paste(names(object@params), '=', unlist(object@params), collapse = ', '))
})


#' @name load_index
#' @title load_index
#' @description load_index
#' @export 
setMethod('initialize', 'BWA', function(.Object, fasta = NULL, seq = NULL, seqname = 'seq',
                                        mc.cores = 1,
                                        hardclip = FALSE,
                                        keep_sec_with_frac_of_primary_score = 0.9,
                                        max_secondary = 10
                                        )
{
    .Object@bwa =.Call(BWA_method( "new" ))
    .Object@params = list(
        mc.cores = mc.cores,
        hardclip = hardclip,
        keep_sec_with_frac_of_primary_score = keep_sec_with_frac_of_primary_score,
        max_secondary = max_secondary)
    if (!is.null(fasta)) {
        BWA__from_fasta(.Object@bwa, fasta)
    }
    else if (!is.null(seq)) {
        BWA__from_string(.Object@bwa, seq, seqname)
    }
    else {
        stop('Either fasta or sequence must be provided')
    }
    return(.Object)
})

BWA_method <- function(name){
     paste( "_RSeqLib_BWA", name, sep = "__" )

}

#' @name BWA
#' @title  BWA
#' @description BWA
#' @export
BWA = function( fasta = NULL, seq = NULL, seqname = 'myseq', mc.cores = 1,
                                        hardclip = FALSE,
                                        keep_sec_with_frac_of_primary_score = 0.9,
                                        max_secondary = 10)
{
    new('BWA', fasta = fasta, seq = seq, seqname = seqname,
        mc.cores,
        hardclip,
        keep_sec_with_frac_of_primary_score,
        max_secondary)
}

#' @name getparam
#' @title getparam
#' @description setparam
#' @export
 setMethod('$', 'BWA', function(x, name)
{
    return(x@params[[name]])
})


#' @name $<-
#' @title $<-
#' @description
#'
#' Setting query params of BWA object 
#'
#' Usage:
#' bwa$mc.cores = 10
#' bwa$max_secondary = 5
#'
#' @param x \code{BWA} object to alter \code{params} field of
#' @param name \code{params} field to alter
#' @param value New value
#' @docType methods
#' @rdname cash-set-methods
#' @aliases $<-,BWA-method
#' @export
#' @author Marcin Imielinski
setMethod('$<-', 'BWA', function(x, name, value)
{
  x@params[[name]] = value
  return(x)
})


#' @name query
#' @title query
#' @description  query
#' @export
setGeneric('query', function(.Object,...) standardGeneric('query'))
setMethod("query", "BWA", function(.Object, qstring, qname = 'myquery',
                                   mc.cores = .Object$mc.cores,
                  hardclip = FALSE,
                  keep_sec_with_frac_of_primary_score = 0.9,
                  max_secondary = 10)
{
    if (!is.null(names(qstring))) {
        qname = names(qstring)
    }
    else if (length(qname)==1) {
        qname = rep(qname, length(qstring))
    }
    tmp = unlist(mcmapply(function(q, qn)
        strsplit(BWA__query(.Object@bwa, q, qn, hardclip,
                   keep_sec_with_frac_of_primary_score,
                   max_secondary), '\n'), qstring, qname, mc.cores = mc.cores))
    
    if (sum(nchar(tmp))==0) {
        return(GRanges())
    }
    
    tmp = .parse_bam(tmp)
    return(tmp)
})


#' @name [
#' @title [
#' @description  query
#' @export
setMethod("[", "BWA", function(x, i) query(x, i))                                      

.parse_bam = function(lines, tags = NULL)
{
    fields = c('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual')
#    lines = gsub('\n', '', lines)   
    linesp = strsplit(lines, '\t')
    chunk = data.table::as.data.table(do.call(rbind, lapply(linesp, function(x) x[1:11])))
    chunk$line = 1:length(chunk$V1)
    ## figure out the optional columns of the BAM read and annotate them. 
    m = munlist(lapply(linesp, function(x) x[-c(1:11)]))    
    tagchunk = fread(paste(m[,3], collapse = '\n'), sep = ':')
    tagchunk$line = as.numeric(m[,1])
    if (is.null(tags))
        tags = unique(tagchunk$V1)
    utags = union(tags, unique(tagchunk$V1))
    tagchunk$V1 = factor(tagchunk$V1, utags)
                                        #    tagchunk = dcast.data.table(tagchunk[tagchunk$V1 %in% tags, ], line ~ V1, value.var = 'V3', fill = NA, drop = FALSE)
    tagchunk = dcast.data.table(tagchunk, line ~ V1, value.var = 'V2', fill = NA, drop = FALSE)
    out = merge(chunk, tagchunk, by = 'line', all.x = TRUE)[, -1, with = FALSE]
    setnames(out, 1:11, fields)
    out = out[, c(fields, tags), with = FALSE]

    cigs <- countCigar(out$cigar)

    out$pos = as.numeric(out$pos)
    
    out$pos2 <- out$pos + rowSums(cigs[, c("D", "M"), drop = FALSE], na.rm=T) - 1
        
    out$qwidth = nchar(out$seq)
    out$strand = bamflag(out$flag)[, "isMinusStrand"] == 1
    out$strand = ifelse(out$strand, "-", "+")
    
    unmapped = bamflag(out$flag)[,"isUnmappedQuery"] == 1
    if (any(unmapped))
    {
        out$pos[unmapped] = 1
        out$pos2[unmapped] = 0
        out$strand[unmapped] = "*"
    }
    
    bf = out$flag
    
    newdt <- data.table(pos = out$pos, pos2 = out$pos2, strand = out$strand, rname = out$rname)   # create data.table of start, end, strand, seqnames
    
    rr <- IRanges(newdt$pos, newdt$pos2)
    sf <- factor(newdt$strand, levels=c('+', '-', '*'))
    ff <- factor(newdt$rname, levels=unique(newdt$rname))
    gr.fields <- c("rname", "strand", "pos", "pos2")
    grobj <- GRanges(seqnames=ff, ranges=rr, strand=sf)
    
    vals = out[, setdiff(names(out), gr.fields), with=FALSE]
    values(grobj) <- vals    
    return(grobj)
}


munlist = function(x, force.rbind = F, force.cbind = F, force.list = F)
  {
    if (!any(c(force.list, force.cbind, force.rbind)))
      {
        if (any(sapply(x, function(y) is.null(dim(y))))) {
            force.list = T
        }
        if (length(unique(sapply(x, function(y) dim(y)[2]))) == 1) {
            force.rbind = T
        }
        if ((length(unique(sapply(x, function(y) dim(y)[1]))) == 1)) {
            force.cbind = T
        }
      }
    else {
          force.list = T
    }

    if (force.list) {
      return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, length(x[[y]])))),
                   iix = unlist(lapply(1:length(x), function(y) if (length(x[[y]])>0) 1:length(x[[y]]) else NULL)),
                   unlist(x)))
    }
    else if (force.rbind) {
      return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))),
                   iix = unlist(lapply(1:length(x), function(y) if (nrow(x[[y]])>0) 1:nrow(x[[y]]) else NULL)),
                   do.call('rbind', x)))
    }
    else if (force.cbind) {
      return(t(rbind(ix = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                     iix = unlist(lapply(1:length(x), function(y) if (ncol(x[[y]])>0) 1:ncol(x[[y]]) else NULL)),
                     do.call('cbind', x))))
    }
  }


#' @name
#' countCigar
#' @title countCigar
#' @description
#'
#' Count bases in cigar string
#'
#' Counts the total number of bases, per cigar, that fall into D, I, M, S categories.
#' countCigar makes no distinction between, for instance 1S2M2S, 2S2M1S, or 3S2M
#' @param cigar character vector of cigar strings
#' @return a 4-column, length(cigar)-row matrix with the total counts for each type
countCigar <- function(cigar) {
    
    cigar.vals <- unlist(strsplit(cigar, "\\d+"))
    cigar.lens <- strsplit(cigar, "[A-Z]")
    lens <- nchar(gsub('\\d+', '', cigar))
    lens[is.na(cigar)] <- 1
    
    cigar.lens <- as.numeric(unlist(cigar.lens))
    cigar.vals <- cigar.vals[cigar.vals != ""]
    repr       <- rep(seq_along(cigar), lens)
    dt         <- data.table(val=cigar.vals, lens=cigar.lens, group=repr, key="val")
    
    smr.d      <- dt["D",][, sum(lens), by=group]
    smr.i      <- dt["I",][, sum(lens), by=group]
    smr.m      <- dt["M",][, sum(lens), by=group]
    smr.s      <- dt["S",][, sum(lens), by=group]
    
    out <- matrix(nrow=length(cigar), ncol=4, 0)
    out[smr.d$group,1] <- smr.d$V1
    out[smr.i$group,2] <- smr.i$V1
    out[smr.m$group,3] <- smr.m$V1
    out[smr.s$group,4] <- smr.s$V1
    colnames(out) <- c('D','I','M','S')
    
    return(out)
}


bamflag = function(reads)
{
    if (inherits(reads, 'GappedAlignments') | inherits(reads, 'data.frame') | inherits(reads, 'GRanges')) {
        bf = reads$flag
    }
    else {
        bf = reads
    }

    out = matrix(as.numeric(intToBits(bf)), byrow = T, ncol = 32)[, 1:12, drop = FALSE]
    colnames(out) = c('isPaired', 'isProperPair', 'isUnmappedQuery', 'hasUnmappedMate', 'isMinusStrand', 'isMateMinusStrand', 'isFirstMateRead', 'isSecondMateRead', 'isNotPrimaryRead', 'isNotPassingQualityControls', 'isDuplicate', 'isSupplementary')

    return(out)
}

