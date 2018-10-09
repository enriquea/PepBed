#' calculate percent overlap between two GRanges objects.
#'
#' @param query GRanges object
#' @param subject GRanges object
#'
#' @return percent overlap between \code{query} and \code{subject}, as defined by the
#'   ratio of the intersection of \code{query} and \code{subject} to the union of
#'   \code{query} and \code{subject}.
#'
#' @details In the this context, \code{query} and \code{subject} are two
#'   GRanges objects. The percent overlap is the number of
#'   nucleotides falling within both query (i.e. genomic features) and subject (i.e. transcripts),
#'   divided by the number of nucleotides falling within either object.
#'   It is essentially Jaccard distance and can be used as a measure of % coverage of qry on sub.
#'
#' @author Function modified from ballgown r-package (Alyssa Frazee)
#'
#' @export
#'
pctOverlap = function(query, subject){
  stopifnot(class(query) == 'GRanges' & class(subject) == 'GRanges')
  ch1 = as.character(runValue(seqnames(query)))
  ch2 = as.character(runValue(seqnames(subject)))
  if(ch1 != ch2){
    return(0)
  }
  tmp1 = reduce(query)
  tmp2 = reduce(subject)
  if(!identical(ranges(tmp1), ranges(query))){
    warning('query contained overlapping ranges and was reduced.')
    query = tmp1
  }
  if(!identical(ranges(tmp2), ranges(tmp2))){
    warning('subject contained overlapping ranges and was reduced.')
    subject = tmp2
  }

  # remove meta-columns
  mcols(query) = NULL
  mcols(subject) = NULL

  ntcov = coverage(c(query, subject))
  ind = which(names(ntcov)==ch1)
  covrle = ntcov[[ind]]
  covrle_val = runValue(covrle)
  covrle_len = runLength(covrle)
    return(sum(covrle_len[covrle_val==2]) / sum(covrle_len[covrle_val==1 | covrle_val == 2]))
}



#' Count number of features by chromosome.
#'
#' @param gr GRanges object
#'
#' @return A dataframe summarizing counts by chromosome.
#'
#' @author Enrique Audain
#'
#' @export
#'
countsByChromosome <- function(gr, colName = 'Peptides') {

        if(length(gr) == 0) {stop('Input GRanges has length zero...')}

         # getting number of elements by chromosome as list...
         l <- lapply(split(gr, seqnames(gr)), function(x) length(x)) # split by chrom and get lenghts
         chrs <- names(l)
         chrs <- gsub(pattern = 'chr', replacement = '', x = chrs, fixed = TRUE)
         counts <- as.vector(unlist(l))
         ## build dataframe with summary
         df <- data.frame(x = as.character(chrs), y = as.numeric(counts))
         names(df) <- c('Chromosome', colName) # rename df
         # sort by Chromosome
             df_ordered <- orderByChromosome(df, colName = 'Chromosome', ref.chr = c(1:22, 'X', 'Y', 'M'))
         # remove row index
             row.names(df_ordered) <- NULL
         return(df_ordered)
}



#' Get the unique features from the original GRanges object.
#'
#' @param gr GRanges object
#' @param colFeatures column name of the features
#'
#' @return A GRanges object with unique features.
#'
#' @author Enrique Audain
#'
#' @export
#'
getUniqueFeatures <- function(gr, colFeatures) {

     # valid GRanges?
      if(class(gr) != 'GRanges' | length(gr) == 0){
        stop('Invalid input GRanges object...')
      }
     # column names exist?
       if(!colFeatures %in% names(mcols(gr))){
          stop('the specified column do not exist...')
       }
      # getting unique features
      all_features <- mcols(gr)[,colFeatures]
      gr_unique <- GenomicRanges::subset(gr, !duplicated(all_features))

     return(gr_unique)
}


#' compute % coverage by chromosome.
#'
#' @param query GRanges object
#' @param subject GRanges object
#' @param colName Set column name with coverage values
#'
#' @return percent coverage of \code{query} on \code{subject} by chromosome.
#'
#' @details This function uses \code{pctOverlap} to compute relative percent overlap.
#'          In this context, interpreted as percent coverage of \code{query} on \code{subject}.
#'
#' @author Enrique Audain
#'
#' @export
#'
computeCoverageByChromosome <- function(query, subject, colName) {

  # getting sequences (chromosomes) info
  seq_qry <- levels(runValue(seqnames(query)))
  seq_sub <- levels(runValue(seqnames(subject)))

  # Are the sequences comparable?
  if (!identical(seq_qry, seq_sub)){
    # stop('Entries contain different chromosome sizes/names...')
  }

  # compute coverage by chromosome
  coverage <- vector()
  for (i in 1:length(seq_qry)){
    chrom <- seq_qry[i]
    coverage[i] <- pctOverlap( query =  GenomicRanges::subset(query, seqnames(query) == chrom),
                               subject = GenomicRanges::subset(subject, seqnames(subject) == chrom))
  }

  # build dataframe chr-coverage
  df <- data.frame(Chromosome = seq_qry, Coverage = round(coverage*100, 3))
  names(df) <- c('Chromosome', colName) # rename df
  df$Chromosome <- gsub(pattern = 'chr', replacement = '', x = df$Chromosome, fixed = TRUE)
  df_ordered <- orderByChromosome(df, colName = 'Chromosome', ref.chr = c(1:22, 'X', 'Y', 'M'))
  row.names(df_ordered) <- NULL # remove row index/names
  return(df_ordered)
}



#' ordering dataframe by crhomosome
#'
#' @param df dataframe
#' @param colName column name with chromosome description
#' @param ref.chr chromosome order
#'
#' @return Return ordered \code{df} by chromosome.
#'
#' @author Enrique Audain
#'
#' @export
#'
orderByChromosome <- function(df, colName, ref.chr = c(1:22, 'X','Y','M')){
  return(df[order(match(as.character(df[,colName]), ref.chr)),])
}


#### count overlapped features ####

## Compute overlap frequency of query on sugject

computeOverlapFrequencyOnSubject <- function(query, subject, sizePopulation = NULL, overlap = 0.50, normalize = TRUE, cCol = 'counts_overlap', fCol = 'frequency' ) {

  #add column of counts. Initially it will be zero...
  mcols(query)[,cCol] <- 0
  mcols(query)[,fCol] <- 0

  # Find self hits
  hits <- findOverlaps(query, subject, type = 'any', select = 'all', ignore.strand = TRUE)

  # If there aren't overlap between features, return the original set
  if(length(hits)==0){
    message('There is not overlapping features...')
    return(query)
  }
  # Keep the hits that achieve the (or passing) pre-defined percent overlap:

  x <- query[queryHits(hits)]
  y <- subject[subjectHits(hits)]
  relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y)) # relative overlap to width union query-subject
  #relative_overlap <- width(pintersect(x, y)) / width(x) # relative overlap to width query
  hits_set <- hits[relative_overlap >= overlap] #keep hits with at least x% overlap...

  # If there are not overlapping up defined threshold, return original set
  if(length(hits_set)==0){
    message('There is not overlapping features at this cutoff...')
    return(query)
  }

  #compute frecuency table
  hits_frecuencies <- as.data.frame(table(queryHits(hits_set)))
  names(hits_frecuencies) <- c("index_query", "counts_overlap")

  #append frecuencies to the original subset
  granges_overlapping <- query[as.integer(hits_frecuencies$index_query)]
  mcols(granges_overlapping)[,cCol] <- hits_frecuencies$counts_overlap

  #retrive non-overlapping regions fron original subset
  granges_non_overlapping <- query[-as.integer(hits_frecuencies$index_query)]

  #build original GRanges object with frequencies
  merged <- c(granges_overlapping, granges_non_overlapping)

  if (!is.null(sizePopulation) & normalize == TRUE ) {
    mcols(merged)[,fCol] <- (mcols(merged)[,cCol]/sizePopulation)
    mcols(merged)[,fCol] <- mcols(merged)[,fCol]/max(mcols(merged)[,fCol]) # normalize based on max F value
    return(merged)
  } else if (!is.null(sizePopulation) & normalize == FALSE ) {
    mcols(merged)[,fCol] <- mcols(merged)[,cCol]/sizePopulation
    return(merged)
  } else {
    mcols(merged)[,fCol] <- NULL
    return(merged)
  }
}



## Compute simple overlap counts of query on sugject
getOverlapCountOnSubject <- function(query, subject, cCol = 'counts'){
  hits <- findOverlaps(query, subject, type = "any", select = "all", ignore.strand = TRUE)
  q <- queryHits(hits)
  freqTable <- as.data.frame(table(q), stringsAsFactors = FALSE) # overlap frequency of query on subject (any overlap)
  overlapping_index <- as.integer(freqTable$q) # retriving index of overlapping features
  query_overlapping <- query[overlapping_index] # getting overlapping features from original query set
  query_non_overlapping <- query[-overlapping_index] # getting non-overlapping features from original query set
  mcols(query_overlapping)[,cCol] <- freqTable$Freq # adding 'number of overlaps on subject' for each query overlapping instance
  mcols(query_non_overlapping)[,cCol] <- 0 # adding 'zero' for each query non-overlapping instance
  query_merged <- c(query_overlapping, query_non_overlapping) # merging
  return(query_merged) # return original query with extra 'counts' column
}

## Convert Granges to DataFrame
granges2dataframe <- function(grange, keep.metacols = TRUE){
  # valid GRanges?
  if(class(grange) != 'GRanges' | length(grange) == 0){
    stop('Invalid input GRanges object...')
  }
     chrs <- as.character(GenomicRanges::seqnames(grange)) # getting chromosome info
     starts <- as.integer(GenomicRanges::start(grange)) # getting starts position
     ends <- as.integer(GenomicRanges::end(grange)) # getting ends position

     # crate basic df
     df <- data.frame(chr=chrs, start=starts, end=ends)

     # checking for meta columns
       if(keep.metacols){
          metacols <- as.data.frame(mcols(grange))
            if(ncol(metacols) != 0){
              df <- cbind(df, metacols) # append meta info
            }
       }
  return(df)
}

