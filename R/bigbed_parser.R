## Simple wrapper functions to parse bigbed/bed files


#' Wrapper to convert bigBed to bed files, using UCSC bedToBigBed tool
#'
#' @param inputFile Path to the bigBed file
#' @param compress Must the output bed file be compressed?
#' @param keep.header If available, Must the header be included in the output bed file?
#' 
#' @details You must have the UCSC bedToBigBed tool installed. 
#' @details UCSC provides pre-built binaries here: http://hgdownload.soe.ucsc.edu/admin/exe/.
#'
#' @export

bigbed2bed <- function(inputFile, compress = TRUE, keep.header = TRUE){
  
  # set output bedfile name path
  output_bed <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", inputFile),'.bed')
  
  # invoke cmd line tool bigbedtobed 
  if (!file.exists(output_bed)) {
    system(paste('bigBedToBed', inputFile, output_bed, sep = '\t'))
       # insert header into the bed output file?
       # if(keep.header){
       #    fields <- getBigBedFieldNames(inputFile = inputFile)
       #      if(!is.null(fields)){
       #         names <- as.character(fields)
       #         # parsing header
       #           temp <- paste(names, collapse = ' ')
       #         # insert header in-place using GNU sed tool
       #         system(paste('sed', '-i', paste("'1 ", "i ", temp, "'", sep = ""), output_bed, sep = '\t'))
       #      }
               # must the output to be compressed
              if(compress){
                system(paste('bgzip', output_bed, '-f', sep = '\t'))
              } 
  } else {
    message('The output bed file already exists...')
  }
}


#' Wrapper to get basic info from bigBed files, using UCSC bigBedInfo tool
#'
#' @param inputFile Path to the bigBed file
#' 
#' @details You must have the UCSC bigBedInfo tool installed. 
#' @details UCSC provides pre-built binaries here: http://hgdownload.soe.ucsc.edu/admin/exe/.
#' 
#' @return This function return a text file summarizing the basic metrics from the input bigBed file.
#'
#' @export

getBigBedInfo <- function(inputFile) {
  
  # set output file path for basic metrics
  output_file <- paste0(sub(pattern = "(.*)\\..*$", replacement = "\\1", inputFile),'_metrics.txt')
  
  # invoke cmd line tool bigBedInfo
  if (!file.exists(output_file)) {
    system(paste('bigBedInfo', inputFile, '-chroms', '>', output_file, sep = '\t'))
  } else {
     message('The output file already exists...')
  }
}


#' Wrapper to get field names and its description from bigBed files, using UCSC bigBedSummary tool
#'
#' @param inputFile Path to the bigBed file
#' @param only.names If TRUE, a character vector with the column names will be returned.
#' 
#' @details You must have the UCSC bigBedSummary tool installed. 
#' @details UCSC provides pre-built binaries here: http://hgdownload.soe.ucsc.edu/admin/exe/.
#' 
#' @return This function return a dataframe with the field names and its description from the input bigBed file.
#'
#' @export

getBigBedFieldNames <- function(inputFile, only.names = TRUE) {

    # invoke cmd line tool bigBedInfo
    object <- system(paste('bigBedSummary', inputFile, '-fields', sep = '\t'), intern = TRUE)
    
    if (is.character(object) & !is.null(object) & length(object) > 1){
        # parse output from bigBedSummary tool
         object <- object[-1] # remove first element (general descrition)
         list.fields <- sapply(object, FUN = strsplit, split='\t', fixed=TRUE) # getting list of fields
         df <- as.data.frame(list.fields, check.names = FALSE)
         df <- as.data.frame(t(df))
         rownames(df) <- NULL
         df[,1] <- NULL
         names(df) <- c('field_name', 'field_description')
           if(only.names){
               return(as.vector(df$field_name))
           } else {
               return(df)
           }
  
    } else {
         message('Not field names available...')
         return(NULL)
    }
}


#' Simple function to read bed files
#'
#' @param inputFile Path to the Bed file
#' 
#' @return This function return a dataframe.
#'
#' @export

readBedFile <- function(inputFile, header = FALSE){
          bed_df <- data.table::fread(input = inputFile, 
                                      header = header, 
                                      data.table = FALSE,
                                      fill = TRUE,
                                      sep = '\t',
                                      quote = '')
          return(bed_df)
}


importBEDasGRange <- function(inputFile){
         grange <- rtracklayer::import.bed(rtracklayer::BEDFile(inputFile))
  return(grange)
}


#' buildGRangesFromData
#'
#'
#' This function builds GRanges Object from data.frame (or file, expected, at least: chr, start, and end columns separeted by tab), 
#' remains columns in the data.frame will be retained as metacolumns.
#'
#' @param data Input dataframe.
#' @param chrColName Column with sequence information.
#' @param startColName Column with start position information.
#' @param endColNmae Column with end position information.
#'
#' @return This function return a GRanges object with extra column like meta-columns.
#'
buildGRangesFromData <- function(data, chrColName = "chrom", startColName = "chromStart", endColName = "chromEnd") {
  
  # rename columns
  if (is.data.frame(data)) {
    
    names(data)[names(data) == chrColName]   <- "seqid"
    names(data)[names(data) == startColName] <- "start"
    names(data)[names(data) == endColName]   <- "end"
    
    # replace sequence identifiers
    id_vector <- c(`1` = "chr1", `2` = "chr2", `3` = "chr3", `4` = "chr4", `5` = "chr5", `6` = "chr6",
                   `7` = "chr7", `8` = "chr8", `9` = "chr9", `10` = "chr10", `11` = "chr11", `12` = "chr12",
                   `13` = "chr13", `14` = "chr14", `15` = "chr15", `16` = "chr16", `17` = "chr17", `18` = "chr18",
                   `19` = "chr19", `20` = "chr20", `21` = "chr21", `22` = "chr22", x = "chrX", y = "chrY", X = "chrX",
                    Y = "chrY", M = "chrM")
    
    # renaming column sequence id values...
    data <- plyr::mutate(data, seqid = plyr::revalue(x = as.character(data$seqid), replace = id_vector))
    
    # create GenomicsRanges
    data_granges <- GenomicRanges::makeGRangesFromDataFrame(df = data,
                                                            seqnames.field = "seqid", 
                                                            start.field = "start", 
                                                            end.field = "end",
                                                            strand.field = "strand",
                                                            keep.extra.columns = TRUE, 
                                                            ignore.strand = FALSE)
    
    # save file with GRanges objects save(data_granges, file = 'granges.rda')
    
    return(data_granges)
  }
  return(NULL)
}
