corn_analysis$analysis_taxa <- purrr::map2(
  .x = corn_analysis$formulas,
  .y = names(corn_analysis$formulas),
  .f = function(formula_, formula_name) {
    
    
    purrr::map2(
      .x = corn_analysis$datasets_taxa,
      .y = names(corn_analysis$datasets_taxa),
      .f = function(dataset, dataset_name) {
        
        
        purrr::map2(
          .x = dataset,
          .y = names(dataset),
          .f = function(tax_lev, tax_lev_name) {
            
            formula_proper_format <- as.formula(formula_name)
            
            
            
            diff_obj_abundance <- corncob::differentialTest(formula = formula_proper_format,
                                                            phi.formula = formula_proper_format,
                                                            formula_null = ~ 1,
                                                            phi.formula_null = formula_proper_format,
                                                            test = "Wald", 
                                                            boot = FALSE,
                                                            data = tax_lev,
                                                            fdr_cutoff = 0.05)
            
            if (length(diff_obj_abundance$significant_taxa) != 0) {
              diff_tax_abundance <- corncob::otu_to_taxonomy(OTU = diff_obj_abundance$significant_taxa, data = tax_lev)
            } else {diff_tax_abundance <- NA}
            
            
            
            diff_obj_dispersion <- corncob::differentialTest(formula = formula_proper_format,
                                                             phi.formula = formula_proper_format,
                                                             formula_null = formula_proper_format,
                                                             phi.formula_null = ~ 1,
                                                             test = "LRT", 
                                                             boot = FALSE,
                                                             data = tax_lev,
                                                             fdr_cutoff = 0.05)
            
            if (length(diff_obj_dispersion$significant_taxa) != 0) {
              diff_tax_dispersion <- corncob::otu_to_taxonomy(OTU = diff_obj_dispersion$significant_taxa, data = tax_lev)
            } else {diff_tax_dispersion <- NA}
            
            
            
            contrasts_da <- contrastsTest(formula = formula_proper_format,
                                                   phi.formula = formula_proper_format,
                                                   contrasts_DA = formula_,
                                                   data = tax_lev,
                                                   fdr_cutoff = 0.99)
            
            contrasts_dv <- corncob::contrastsTest(formula = formula_proper_format,
                                                   phi.formula = formula_proper_format,
                                                   contrasts_DV = formula_,
                                                   data = tax_lev,
                                                   fdr_cutoff = 0.99)
            
            return(list('diff_obj_abundance' = diff_obj_abundance, 'diff_tax_abundance' = diff_tax_abundance, 'diff_obj_dispersion' = diff_obj_dispersion, 'diff_tax_dispersion' = diff_tax_dispersion, 'contrasts_da' = contrasts_da, 'contrasts_dv' = contrasts_dv))
          })
      })
  })
