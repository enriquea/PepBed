## Simple routines to extract features from GTF file

## import GTF file
## valid ftp gencode url (Version 27 (January 2017 freeze, GRCh38) - Ensembl 90)
gtf <- rtracklayer::import(con = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz")

## getting protein coding transcripts as GRanges object
transcript <- subset(gtf, GenomicRanges::mcols(gtf)$type == 'transcript' & GenomicRanges::mcols(gtf)$transcript_type == 'protein_coding')

## save protein coding transcripts as GRanges
save(transcript, file = 'data/protein_coding_transcript.rda')

## getting proteing coding genes
genes <- subset(gtf, GenomicRanges::mcols(gtf)$type == 'gene' & GenomicRanges::mcols(gtf)$gene_type == 'protein_coding')

## save protein coding genes as GRanges
save(genes, file = 'data/protein_coding_genes.rda')
