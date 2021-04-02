phyloseq$datasets$group_effect_rm_dada2$deseq <- purrr::map(
  .x = list('mother' = 'mother', 'cluster' = 'cluster'), 
  .f = function(factor){
    
    matrix <- phyloseq::phyloseq_to_deseq2(
      physeq = phyloseq$datasets$group_effect_rm_dada2$phylo_filtered,
      design = as.formula(paste0('~ ', factor)))
    
    deseq <- DESeq2::DESeq(matrix, test = "LRT", reduced = ~ 1)
    
    lrt <- DESeq2::results(deseq)
    
    results <- get_nice_results_from_lrt_deseq2(results_lrt = lrt, phyloseq_object = phyloseq[["datasets"]][["group_effect_rm_dada2"]][["phylo_filtered"]])
    
    return(list('matrix' = matrix, 'deseq' = deseq, 'lrt' = lrt, 'results' = results))
  })