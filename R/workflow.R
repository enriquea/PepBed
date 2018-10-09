## Build stats report from peptide BED file

generate_pepbed_report <- function(bedfile = NULL, ref_genome = NULL, output_file = 'report.pdf'){
       
  # import bed file as GRanges
    granges_peptide <- importBEDasGRange(inputFile = bedfile)
       
      # print number of features imported
      n_features <- length(granges_peptide)
      message(c('Imported ', n_features, ' peptides...'))
    
  # getting number of features(peptides) by chromosome
    counts <- countsByChromosome(gr = granges_peptide, colName = 'Number of peptides')
    
  # removing duplicated entries from original granges_peptide
    unique_pep <- getUniqueFeatures(granges_peptide, colFeatures = 'name')
    
  # getting number of features(peptides) by chromosome
    counts_unique <- countsByChromosome(gr = unique_pep, colName = 'Number of peptides')
    
  genome_grange <- NULL
  
  # compute coverage of query (peptide evidences) on subject (transcripts) by crhomosome
    
  message(c('Computing coverage on genome ensembl ', ref_genome, '...'))
  
    if(ref_genome == 'hg38'){
      data("protein_coding_transcript_hg38") # load protein coding transcript as GRanges object
      coverage <- computeCoverageByChromosome(query = granges_peptide, 
                                              subject = transcripts_hg38, 
                                              colName = 'Coverage')
      
    } else {
      data("protein_coding_transcript_hg19")
      coverage <- computeCoverageByChromosome(query = granges_peptide, 
                                              subject = transcripts_hg19, 
                                              colName = 'Coverage')
      
    }

    # plotting covergae per cromosome
    plot <- ggplot(coverage, aes(x = factor(Chromosome, levels = unique(coverage$Chromosome)), y = Coverage)) + 
      geom_col(fill='darkgreen', alpha=0.4) +
      labs(x = 'Chromosome', y = '% coverage', fill = '') +
      theme_bw() +
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=14),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"))
    
    message('Rendering circular genome density plot, It could take several minutes...')
    
    # building circular plotin external device
      temporal_file <- paste0(getwd(),'/','circos_plot.png')
      png(filename = temporal_file, width = 1600, height = 1600, res = 300)
      
      if(ref_genome == 'hg19'){
        circos.initializeWithIdeogram(species='hg19') # initialize circos plot
      } else {
        circos.initializeWithIdeogram(species='hg38', chromosome.index = paste0("chr", c(1:22, "X", "Y")))
      }
        # convert to basic data frame
        pepbed_df <- granges2dataframe(grange = granges_peptide, keep.metacols = FALSE)
        # plot density
        circos.genomicDensity(pepbed_df, col = c("#FF000080"), track.height = 0.1, baseline = 0)
        circos.clear()
        
      dev.off()
      
    message('Building final report...')
      
    ## rendering PDF report
      # getting the default RMD template to render document
      report_template <- system.file('extdata', 'Template/report_template.Rmd', package = 'PepBed')
      
      rmarkdown::render(input = report_template,
                        output_dir = getwd(),
                        output_file = output_file,
                        output_format = 'pdf_document')
      
    ## cleaning up temporal file(s)
      if (file.exists(temporal_file)) file.remove(temporal_file)
    
    message('Done!')
}