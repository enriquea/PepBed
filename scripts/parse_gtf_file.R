## Simple routines to extract features from GTF file(s)

## import GTF file(s)

## valid ftp gencode url (Version 27 (January 2017 freeze, GRCh38) - Ensembl 90)
gtf <- rtracklayer::import(con = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz")

## getting protein coding transcripts as GRanges object
transcripts_hg38 <- subset(gtf, GenomicRanges::mcols(gtf)$type == 'transcript' & GenomicRanges::mcols(gtf)$transcript_type == 'protein_coding')

## save protein coding transcripts as GRanges
save(transcripts_hg38, file = 'data/protein_coding_transcript_hg38.rda')

## getting proteing coding genes
genes_hg38 <- subset(gtf, GenomicRanges::mcols(gtf)$type == 'gene' & GenomicRanges::mcols(gtf)$gene_type == 'protein_coding')

## save protein coding genes as GRanges
save(genes_hg38, file = 'data/protein_coding_genes_hg38.rda')

## getting "all exons" from protein-coding transcripts
exons_hg38 <- subset(gtf, GenomicRanges::mcols(gtf)$type == 'exon' & GenomicRanges::mcols(gtf)$transcript_type == 'protein_coding')

## save protein coding exons (and isoforms) as GRanges (~700,000 exons)
save(exons_hg38, file = 'data/protein_coding_exons_hg38.rda')

################################################ downloading/parsing GTF file for mouse (mm10) ##############################################################

## valid ftp gencode url (Mouse Genome Version M10 (mm10))
gtf <- rtracklayer::import(con = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M10/gencode.vM10.basic.annotation.gtf.gz")

## getting protein coding transcripts as GRanges object
transcripts_mm10 <- subset(gtf, GenomicRanges::mcols(gtf)$type == 'transcript' & GenomicRanges::mcols(gtf)$transcript_type == 'protein_coding')

## save protein coding transcripts (mm10) as GRanges
save(transcripts_mm10, file = 'data/protein_coding_transcript_mm10.rda')

## getting proteing coding genes (mm10)
genes_mm10 <- subset(gtf, GenomicRanges::mcols(gtf)$type == 'gene' & GenomicRanges::mcols(gtf)$gene_type == 'protein_coding')

## save protein coding genes as GRanges
save(genes_mm10, file = 'data/protein_coding_genes_mm10.rda')
