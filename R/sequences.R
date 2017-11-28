
#' Giving a dataframe (in BED format) with modified peptide sequences return the list of modifications.
#' Modification ID expected to be in brackets.
#'
#' @param bedmod Dataframe in BED format with modified peptide sequences
#' @param colName Name of the column with sequences
#' 
#' @return A character vector with modification description(ID).
#' 
#' @author Enrique Audain
#'
#' @export
#' 
getModifications <- function(bedmod, colName) {
       seq <- as.vector(bedmod[,colName])
       mods <- unique(unlist(lapply(seq, function(x) gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]))))
       return(mods)
}



#' Giving a dataframe (in BED format) with modified peptide sequences,
#' return a dataframe containing sequences with (at least) the specified
#' modification (expected to be in brackets).
#'
#' @param bedmod Dataframe in BED format with modified peptide sequences
#' @param colName Name of the column with sequences
#' @param modPattern The modification name (as in the sequence)
#' 
#' @return A dataframe 
#' 
#' @author Enrique Audain
#'
#' @export
#' 
getModifiedSeq <- function(bedmod, colName, modPattern = 'phospho'){
       temp <- subset(bedmod, grepl(x = bedmod[,colName], pattern = modPattern, fixed = TRUE))
}



#' Giving a dataframe (in BED format) with modified peptide sequences,
#' return the number of modified features by chromosome.
#'
#' @param bedmod Dataframe in BED format with modified peptide sequences
#' @param colName Name of the column with sequences
#' @param modPattern The modification name (as in the sequence)
#' 
#' @return A dataframe summarizing the number of modified sequences by chromosome
#' 
#' @author Enrique Audain
#'
#' @export
#' 
summariseModByChromosome <- function(bedmod, colName, modPattern = 'phospho'){
       subset_df <- getModifiedSeq(bedmod, colName, modPattern = modPattern)
       l <- split.data.frame(subset_df, f = subset_df[,1])
       counts <- lapply(l, nrow)
       output <- data.frame(chromosome = names(counts), mod = as.vector(unlist(counts)))
       names(output) <- c('Chromosome', modPattern)
       output$Chromosome <- gsub(pattern = 'chr', output$Chromosome, replacement = '', fixed = TRUE)
       return(output)
}


#' Giving a dataframe (in BED format) with modified peptide sequences,
#' return the number of modified features by chromosome for every modification
#' present in the input data.
#'
#' @param bedmod Dataframe in BED format with modified peptide sequences
#' @param colName Name of the column with sequences
#' 
#' @return A dataframe with all modifications (counts) by chromosome
#' 
#' @author Enrique Audain
#'
#' @export
#' 
summariseBedMod <- function(bedmod, colName) {
  # gettting all present modification
  modPatterns <- getModifications(bedmod, colName)
  # split data by modifications
  l <- lapply(modPatterns, function(x) summariseModByChromosome(bedmod, colName=colName, modPattern = x))
  # merged data frames by Chromosome
  merged <- Reduce(function(df1, df2) merge(df1, df2, by = "Chromosome", all = TRUE), l)
  # replace NA by zeros
  merged[is.na(merged)] <- 0
  return(merged)
}
