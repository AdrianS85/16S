# phyloseq::psmelt - phyloseq to df
# BiocManager::install("limma")

library(corncob)
library(phyloseq)
library(magrittr)
data(soil_phylo)

corn_analysis <- list(
  'import_dir' = list(
    'feat' = '/home/adrian/Desktop/nx_temp/nx_data_q/', 
    'tree' = '/home/adrian/Desktop/nx_temp/nx_phylogeny/',
    'tax' = '/home/adrian/Desktop/nx_temp/nx_taxonomy/',
    'metadata' = '/home/adrian/Desktop/qiime/qiime_metadata_nf.tsv'),
  'formulas' = list('~ group' = list("groupanemia-groupfe_dextran_muscle", "groupanemia-groupsiderai_fe_liposomal_oral", 'groupfe_dextran_muscle-groupsiderai_fe_liposomal_oral')),
  'corncob_dir' = 'corncob',
  'vis_dir' = 'corncob/vis'
)
  
dir.create(corn_analysis$corncob_dir)
dir.create(corn_analysis$vis_dir)

### !!! add sklearn!, create and use non-filtered taxonomies
corn_analysis$datasets <- list(
  'dada2' = qiime2R::qza_to_phyloseq(
    features = paste0(corn_analysis$import_dir$feat, 'filtered_sepp_dada2.FeatureTableFrequency.qza'), 
    tree = paste0(corn_analysis$import_dir$tree, 'sepp_dada2.PhylogenyRooted.qza'),
    taxonomy = paste0(corn_analysis$import_dir$tax, 'vsearch_filtered_sepp_dada2.FeatureDataTaxonomy.qza'),
    metadata = corn_analysis$import_dir$metadata),
  'deblur' = qiime2R::qza_to_phyloseq(
    features = paste0(corn_analysis$import_dir$feat, 'filtered_sepp_deblur.FeatureTableFrequency.qza'), 
    tree = paste0(corn_analysis$import_dir$tree, 'sepp_deblur.PhylogenyRooted.qza'),
    taxonomy = paste0(corn_analysis$import_dir$tax, 'vsearch_filtered_sepp_deblur.FeatureDataTaxonomy.qza'),
    metadata = corn_analysis$import_dir$metadata),
  'deblur_s' = qiime2R::qza_to_phyloseq(
    features = paste0(corn_analysis$import_dir$feat, 'filtered_sepp_deblur_singletons.FeatureTableFrequency.qza'), 
    tree = paste0(corn_analysis$import_dir$tree, 'sepp_deblur_singletons.PhylogenyRooted.qza'),
    taxonomy = paste0(corn_analysis$import_dir$tax, 'vsearch_filtered_sepp_deblur_singletons.FeatureDataTaxonomy.qza'),
    metadata = corn_analysis$import_dir$metadata)
)






corn_analysis$datasets_taxa <- purrr::map(
  .x = corn_analysis$datasets,
  .f = function(dataset) {
    
    
    purrr::map(
      .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'),
      .f = function(tax_level) {
        
        phyloseq::tax_glom(physeq = dataset, taxrank = tax_level) 
      })
  })





### !!! does not catch formulas as it should
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
                
                
                
                contrasts_da <- corncob::contrastsTest(formula = formula_proper_format,
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












