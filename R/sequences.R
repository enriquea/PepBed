
#' Giving a dataframe (in BED format) with modified peptide sequences return the list of modifications.
#' Modification ID expected to be in brackets.
#'
#' @param bedmod Dataframe in BED format with modified peptide sequences
#' @param colName Name of the column with sequences
#' 
#'
#' @export
#' 
getModifications <- function(bedmod, colName) {
       seq <- as.vector(bedmod[,colName])
       mods <- unique(unlist(lapply(seq, function(x) gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]]))))
       return(mods)
}



#' Giving a dataframe (in BED format) with modified peptide sequences,
#' return the list of sequences with a specific modification (expected to be in brackets).
#'
#' @param bedmod Dataframe in BED format with modified peptide sequences
#' @param colName Name of the column with sequences
#' @param modPattern The modification name (as in the sequence)
#' 
#'
#' @export
#' 
getModifiedSeq <- function(bedmod, colName, modPattern = 'phospho'){
       temp <- subset(bedmod, grepl(x = bedmod[,colName], pattern = modPattern, fixed = TRUE))
}


# Count the number of modifications by chromosome
summariseModByChromosome <- function(df, colName, modPattern = 'phospho'){
       subset_df <- getModifiedSeq(df, colName, modPattern = modPattern)
       l <- split.data.frame(subset_df, f = subset_df[,1])
       counts <- lapply(l, nrow)
       output <- data.frame(chromosome = names(counts), mod = as.vector(unlist(counts)))
       names(output) <- c('Chromosome', modPattern)
       output$Chromosome <- gsub(pattern = 'chr', output$Chromosome, replacement = '', fixed = TRUE)
       return(output)
}


# summarising bed file by modifications
summariseModBed <- function(modbed, colName) {
  
  # gettting all present modification
  modPatterns <- getModifications(modbed, colName)
  
  # split data by modifications
  l <- lapply(modPatterns, function(x) summariseModByChromosome(modbed, colName=colName, modPattern = x))
    
  # merged data frames by Chromosome
  merged <- Reduce(function(df1, df2) merge(df1, df2, by = "Chromosome", all = TRUE), l)
  
  # replace NA by zeros
  merged[is.na(merged)] <- 0
  
  return(merged)
}
