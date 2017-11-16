## read bigbed file

# path to bigbed file
bigbedpath <- paste(getwd(), '/inst/extdata/pride_cluster_peptides_9606_Human.pogo.bb', sep = '')

# convert bigbed to bed file (same directory for output)
bigbed2bed(inputFile = bigbedpath, compress = FALSE)

# getting basic information
getBigBedInfo(inputFile = bigbedpath)

# getting field names if available
fieldNames <- getBigBedFieldNames(inputFile = bigbedpath, only.names = TRUE)

# import bed file as  dataframe
bedpath <- paste(getwd(), '/inst/extdata/pride_cluster_peptides_9606_Human.bed', sep = '')
df <- readBedFile(inputFile = bedpath)

# rename columns
names(df) <- fieldNames

# convert dataframe to GRanges
granges_peptide <- buildGRangesFromData(data = df, chrColName = "chrom", startColName = "chromStart", endColName = "chromEnd") 


## getting number of features(peptides) by chromosome
stats <- countsByChromosome(gr = granges_peptide)
print(stats)

## compute coverage of query (peptide evidences) on subject (transcripts) by crhomosome
load(file = 'data/protein_coding_transcript.rda') # load protein coding transcript as GRanges object
coverage <- computeCoverageByChromosome(query = granges_peptide, subject = transcript)
print(coverage)
