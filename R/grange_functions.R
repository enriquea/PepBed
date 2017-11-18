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
countsByChromosome <- function(gr) {
      ## getting seq info
       rle <- seqnames(gr)
       if(length(rle) != 0){
         ## build dataframe with summary
         df <- data.frame(chromosome = levels(runValue(rle)), countFeatures = runLength(rle))
         return(df)
       } else {
         return(NULL) # print msg here!
       }
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
computeCoverageByChromosome <- function(query, subject) {
  
  # getting sequences (chromosomes) info
  seq_qry <- levels(runValue(seqnames(query)))
  seq_sub <- levels(runValue(seqnames(subject)))
  
  # Are the sequences comparable?
  if (!identical(seq_qry, seq_sub)){
    stop('Entries contain different chromosome sizes/names...')
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
  return(df)
}