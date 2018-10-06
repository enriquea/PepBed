library(PepBed)
library(VennDiagram)
library(ggplot2)
library(gplots)

bedFileDataPepGenome <- readBedFile("/Users/yperez/IdeaProjects/github-repo/BDP/benchmark/PepGenome-PeptideAtlas.bed");
bedFileDataPepGenome <- setBed12Columns(bedFileDataPepGenome);

bedFileDataPGx <- readBedFile("/Users/yperez/IdeaProjects/github-repo/BDP/benchmark/PGx-PeptideAtlas-Human.bed")
bedFileDataPGx <- setBed12Columns(bedFileDataPGx);

#Generate Venn Diagram Overlap the overlap of peptides.
peptideLists <-list(PepGenome = bedFileDataPepGenome$name, PGx = bedFileDataPGx$name)

venn.plot <- venn.diagram(peptideLists , filename = "../benchmark/PepGenome.png", fill=c("red", "green"), alpha=c(0.2,0.2), cex = 1, cat.fontface=1, cat.cex=1, category.names=c("PepGenome", "PGx"), main="Number of Peptides map to the Human ENSEMBL Genome", scaled = FALSE, cat.default.pos='text')
