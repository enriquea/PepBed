#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-infile", "--infile"), type="character", default=NULL, help="Bed file to be analized", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="pedbedr-out.pdf",
                help="pdf file ouput [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$infile)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

bedFilePath <- opt$infile
pdfFile <- opt$out

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### Parsing bigbed files

```{r , message=FALSE, warning=FALSE}

library(PepBed)

### Parsing Bed file

```{r , message=FALSE, warning=FALSE}


# import bed file as  dataframe
bed_df <- readBedFile(inputFile = bedFilePath)
bed_df <- setBed12Columns(bed_df)

# convert dataframe to GRanges
# all non-modified peptides
granges_peptide <- buildGRangesFromData(data = bed_df,
                                        chrColName = "chrom",
                                        startColName = "chromStart",
                                        endColName = "chromEnd")

```

\newpage

### Computing some basic stats from the data

```{r , message=FALSE, warning=FALSE}
# getting number of features(peptides) by chromosome
counts <- countsByChromosome(gr = granges_peptide, colName = 'Peptides')
counts_mod <- countsByChromosome(gr = granges_mod_peptide, colName = 'Peptides_mod')

# merging dfs
merged_counts <- merge.data.frame(counts, counts_mod, by = 'Chromosome')
# ordering by chromosome
merged_counts <- orderByChromosome(df = merged_counts, colName = 'Chromosome')

print(merged_counts)
```

\newpage

### Getting stats for unique peptides

```{r , message=FALSE, warning=FALSE}
# removing duplicated entries from original granges_peptide
unique_pep <- getUniqueFeatures(granges_peptide, colFeatures = 'name')
unique_pep_mod <- getUniqueFeatures(granges_mod_peptide, colFeatures = 'name')

# getting unique number of features(peptides) by chromosome
counts_unique <- countsByChromosome(gr = unique_pep, colName = 'Peptides')
counts_mod_unique <- countsByChromosome(gr = unique_pep_mod, colName = 'Peptides_mod')

# merging dfs
merged_counts_unique <- merge.data.frame(counts_unique,
                                         counts_mod_unique,
                                         by = 'Chromosome')

# ordering by chromosome
merged_counts_unique <- orderByChromosome(df = merged_counts_unique,
                                          colName = 'Chromosome')

print(merged_counts_unique)
```

\newpage

### Computing % coverage

```{r , message=FALSE, warning=FALSE}
## compute coverage of query (peptide evidences) on subject (transcripts) by crhomosome
data("protein_coding_transcript_hg38") # load protein coding transcript as GRanges object

coverage <- computeCoverageByChromosome(query = granges_peptide,
                                        subject = transcripts_hg38,
                                        colName = 'Coverage')

coverage_mod <- computeCoverageByChromosome(query = granges_mod_peptide,
                                            subject = transcripts_hg38,
                                            colName = 'Coverage_mod')

# merging dfs
merged_coverage <- merge.data.frame(coverage, coverage_mod, by = 'Chromosome')

# ordering by chromosome
merged_coverage <- orderByChromosome(df = merged_coverage, colName = 'Chromosome')

print(merged_coverage)
```

\newpage

### Visualizing the data

* The distribution of peptides by chromosome. (blue_track: modified peptide;  red_track: non-modified)

```{r density, echo=FALSE, fig.width=12, fig.height=12}
library(circlize)
circos.initializeWithIdeogram(species = 'hg19')
bed <- bed_df
bed_mod <- bed_mod_df
circos.genomicDensity(bed, col = c("#FF000080"), track.height = 0.1, baseline = 0)
circos.genomicDensity(bed_mod, col = c("#0000FF80"), track.height = 0.1, baseline = 0)
circos.clear()
```

\newpage

* barplot with number of peptides (modified and non-modified) by chromosome

```{r counts, echo=FALSE, message=FALSE, warning=FALSE, fig.width=12, fig.height=8}

dat <- reshape2::melt(merged_counts)

plot1 <- ggplot(dat, aes(x = factor(Chromosome, levels = unique(dat$Chromosome)), y = value, fill=variable)) +
    geom_col(position = 'dodge', alpha = 0.7) +
    labs(x = 'Chromosome', y = 'Number of Peptides', fill = '') +
    scale_fill_discrete("Peptides", labels=c("non-modified", "modified")) +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
plot1
```

\newpage

* barplot with number of unique peptides (modified and non-modified) by chromosome

```{r unique_counts, echo=FALSE, message=FALSE, warning=FALSE, fig.width=12, fig.height=8}

dat <- reshape2::melt(merged_counts_unique)

plot2 <- ggplot(dat, aes(x = factor(Chromosome, levels = unique(dat$Chromosome)), y = value, fill=variable)) +
    geom_col(position = 'dodge', alpha = 0.7) +
    labs(x = 'Chromosome', y = 'Number of Peptides', fill = '') +
    scale_fill_discrete("Unique peptides", labels=c("non-modified", "modified")) +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
plot2
```

\newpage

* barplot with coverage (all peptides) by chromosome

```{r coverage, echo=FALSE, message=FALSE, warning=FALSE, fig.width=12, fig.height=8}

dat <- coverage

plot3 <- ggplot(dat, aes(x = factor(Chromosome, levels = unique(dat$Chromosome)), y = Coverage)) +
    geom_col(fill='darkgreen', alpha=0.4) +
    labs(x = 'Chromosome', y = '% coverage', fill = '') +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
plot3
```

