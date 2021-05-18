
phyloseq$muscleDextran$deseq_tax_muscleDextran <- purrr::map(
  .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'), 
  .f = function(taxon_level){
    glom_ <- phyloseq::tax_glom(physeq = phyloseq[["datasets"]][["dada2"]][["phylo_filtered"]], taxrank = taxon_level)
    
    matrix_ <- phyloseq::phyloseq_to_deseq2(
      physeq = glom_,
      design = ~mother + group)
    
    deseq_ <- DESeq2::DESeq(matrix_, test = "Wald")
    
    lrt_ <- DESeq2::results(deseq_, contrast = c('group', 'fe_dextran_muscle', 'anemia'))
    
    if (length(lrt_[[1]]) == 0) {
      return('result' = NA)
    }
    
    result_ = cbind(as(lrt_, "data.frame"), as(phyloseq::tax_table(glom_)[rownames(lrt_), ], "matrix"))
    
    openxlsx::addWorksheet(
      wb = phyloseq$muscleDextran$deseq_wb, 
      sheetName = paste0(taxon_level, "_", factor))
    
    openxlsx::writeDataTable(
      wb = phyloseq$muscleDextran$deseq_wb, 
      sheet = paste0(taxon_level, "_", factor),
      x = subset(result_, subset = result_$pvalue <0.05))
    
    return('result' = result_)
  })