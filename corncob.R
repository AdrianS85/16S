# phyloseq::psmelt - phyloseq to df
# BiocManager::install("limma")
# BiocManager::install("MMUPHin")
furrr

### !!! To avoid testing rare ASVs with our smaller sample size, we removed ASVs that did not have a mean count at or above the 25th percentile in at least 10% of samples. [Timing of complementary feeding is associated with gut microbiota diversity and composition and short chain fatty acid concentrations over the first year of life]

# try batch correct with this: https://bioconductor.org/packages/release/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html#performing-batch-study-effect-adjustment-with-adjust_batch

library(corncob)
library(phyloseq)
library(magrittr)


load('corn_analysis.save')




corn_analysis <- list(
  'import_dir' = list(
    'feat' = '/home/adrian/Desktop/qiime/full_data_server/nx_data_q/', 
    'tree' = '/home/adrian/Desktop/qiime/full_data_server/nx_phylogeny/',
    'tax' = '/home/adrian/Desktop/qiime/full_data_server/nx_taxonomy/',
    'metadata' = '/home/adrian/Desktop/qiime/metadata_enhanced_phylo.tsv'),
  'formulas' = list(
    '~ group' = list("groupfe_dextran_muscle", "groupsiderai_fe_liposomal_oral", "groupnanofe_dextran_50nm_oral", "groupnanofe_phospholipid_70nm_oral", "groupnanofe_synomag_30nm_oral", "groupfeso4_oral"),
    '~ anemia_v_treatment' = list("anemia_v_treatmenttreatment"),
    '~ dex_muscle_v_treatment' = list("dex_muscle_v_treatmentnon_dextran_muscle_treatment")),
  'sample_to_remove_for_formulas' = 
    list(
      '~ group' = NA,
      '~ anemia_v_treatment' = NA,
      '~ dex_muscle_v_treatment' = c('A1', 'A2', 'A3', 'A4', 'A5')),
  'corncob_dir' = 'corncob',
  'vis_dir' = 'corncob/vis'
)
# dex_muscle_v_treatment does not converge








dir.create(corn_analysis$corncob_dir)
dir.create(corn_analysis$vis_dir)


corn_analysis$datasets_base <- list(
  'dada2' = qiime2R::qza_to_phyloseq(
    features = paste0(corn_analysis$import_dir$feat, 'filtered_sepp_dada2.FeatureTableFrequency.qza'), 
    tree = paste0(corn_analysis$import_dir$tree, 'sepp_dada2.PhylogenyRooted.qza'),
    taxonomy = paste0(corn_analysis$import_dir$tax, 'vsearch_filtered_sepp_dada2.FeatureDataTaxonomy.qza'),
    metadata = corn_analysis$import_dir$metadata),
  'deblur' = qiime2R::qza_to_phyloseq(
    features = paste0(corn_analysis$import_dir$feat, 'filtered_sepp_deblur.FeatureTableFrequency.qza'), 
    tree = paste0(corn_analysis$import_dir$tree, 'sepp_deblur.PhylogenyRooted.qza'),
    taxonomy = paste0(corn_analysis$import_dir$tax, 'vsearch_filtered_sepp_deblur.FeatureDataTaxonomy.qza'),
    metadata = corn_analysis$import_dir$metadata)
)

### !!! is metadata imported as factors?




temp_batch_dada2 <- MMUPHin::adjust_batch(
  feature_abd = corn_analysis$datasets_base$dada2@otu_table, 
  batch = 'mother', 
  data = data.frame(corn_analysis$datasets_base$dada2@sam_data),
  control = list('zero_inflation' = T, 'pseudo_count' = NULL, 'diagnostic_plot' = paste0(corn_analysis$vis_dir, '/dada2_batch_adj.pdf')))[['feature_abd_adj']]

temp_batch_deblur <- MMUPHin::adjust_batch(
  feature_abd = corn_analysis$datasets_base$deblur@otu_table, 
  batch = 'mother', 
  data = data.frame(corn_analysis$datasets_base$deblur@sam_data),
  control = list('zero_inflation' = T, 'pseudo_count' = NULL, 'diagnostic_plot' = paste0(corn_analysis$vis_dir, '/deblur_batch_adj.pdf')))[['feature_abd_adj']]

corn_analysis$datasets_base$batch_corr_dada2 <- phyloseq::phyloseq(
  phyloseq::otu_table(temp_batch_dada2, taxa_are_rows = T), 
  phyloseq::phy_tree(corn_analysis$datasets_base$dada2@phy_tree), 
  phyloseq::tax_table(corn_analysis$datasets_base$dada2@tax_table), 
  phyloseq::sample_data(corn_analysis$datasets_base$dada2@sam_data))

corn_analysis$datasets_base$batch_corr_deblur <- phyloseq::phyloseq(
  phyloseq::otu_table(temp_batch_deblur, taxa_are_rows = T), 
  phyloseq::phy_tree(corn_analysis$datasets_base$deblur@phy_tree), 
  phyloseq::tax_table(corn_analysis$datasets_base$deblur@tax_table), 
  phyloseq::sample_data(corn_analysis$datasets_base$deblur@sam_data))














corn_analysis$datasets <- purrr::map(
  .x = corn_analysis$datasets_base,
  .f = function(dataset) {
    
    list_taxons <- purrr::map(
      .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'),
      .f = function(tax_level) {
        
        phyloseq::tax_glom(physeq = dataset, taxrank = tax_level) 
      })
    
    list_taxons[['asv']] <- dataset
    
    return(list_taxons)
  })




corn_analysis$datasets_only_most_interesting$batch_corr_dada2 <- corn_analysis$datasets$batch_corr_dada2

save(corn_analysis, file= 'corn_analysis.save')



load('corn_analysis.save')

corn_analysis$analysis <- purrr::map2(
  .x = corn_analysis$formulas,
  .y = names(corn_analysis$formulas),
  .f = function(formula_, formula_name) {
    
    library(corncob)
    library(phyloseq)
    
    
    
    test <- purrr::map2(
      .x = corn_analysis$datasets_only_most_interesting,
      .y = names(corn_analysis$datasets_only_most_interesting),
      .f = function(dataset, dataset_name) {
        
        
        
        furrr::furrr_options(seed = 1)
        
        return_data <- furrr::future_map2(
          .x = dataset,
          .y = dataset_name,
          .f = function(tax_lev, tax_lev_name) {
            
            
            formula_proper_format <- as.formula(formula_name)
            
            
            
            if (formula_name == '~ dex_muscle_v_treatment') {
              tax_lev = phyloseq::subset_samples(tax_lev, !is.na(dex_muscle_v_treatment))
            }
            
            tryCatch(
              {
                diff_obj_abundance <- corncob::differentialTest(formula = formula_proper_format,
                                                                phi.formula = formula_proper_format,
                                                                formula_null = ~ 1,
                                                                phi.formula_null = formula_proper_format,
                                                                test = "Wald", 
                                                                boot = T, 
                                                                B = 999,
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
                                                                 boot = T,
                                                                 B = 999,
                                                                 data = tax_lev,
                                                                 fdr_cutoff = 0.05)
                
                if (length(diff_obj_dispersion$significant_taxa) != 0) {
                  diff_tax_dispersion <- corncob::otu_to_taxonomy(OTU = diff_obj_dispersion$significant_taxa, data = tax_lev)
                } else {diff_tax_dispersion <- NA}
                
                
                contrasts_da <- corncob::contrastsTest(formula = formula_proper_format,
                                                       phi.formula = formula_proper_format,
                                                       contrasts_DA = formula_,
                                                       data = tax_lev,
                                                       fdr_cutoff = 0.05)
                
                contrasts_dv <- corncob::contrastsTest(formula = formula_proper_format,
                                                       phi.formula = formula_proper_format,
                                                       contrasts_DV = formula_,
                                                       data = tax_lev,
                                                       fdr_cutoff = 0.05)
                
                diff_obj_dispersion <- NA
                diff_tax_dispersion <- NA
                contrasts_dv <- NA
                
                return(list('diff_obj_abundance' = diff_obj_abundance, 'diff_tax_abundance' = diff_tax_abundance, 'diff_obj_dispersion' = diff_obj_dispersion, 'diff_tax_dispersion' = diff_tax_dispersion, 'contrasts_da' = contrasts_da, 'contrasts_dv' = contrasts_dv))
                
              },
              error=function(cond) {
                message(paste0(cond, '\t for \t', tax_lev_name, ' ', formula_name))
                
                return(NA)
              })
          })
        return(return_data)
      })
  })

save(corn_analysis, file= 'corn_analysis_from_server.save')




# B = 1: ~40min
# B = 10: ~180min
# B = 1 (1 group): 

# corn_analysis$json$group_batch_corr_dada2 <- jsonlite::toJSON(corn_analysis$analysis_taxa$`~ group`$batch_corr_dada2, force = T)
# corn_analysis$json$anemia_v_treatment_batch_corr_dada2 <- jsonlite::toJSON(corn_analysis$analysis_taxa$`~ anemia_v_treatment`$batch_corr_dada2, force = T)
# corn_analysis$json$dex_muscle_v_treatment_batch_corr_dada2 <- jsonlite::toJSON(corn_analysis$analysis_taxa$`~ dex_muscle_v_treatment`$batch_corr_dada2, force = T)
# 
# write(corn_analysis$json$group_batch_corr_dada2, "corncob_group_batch_corr_dada2.json")
# write(corn_analysis$json$anemia_v_treatment_batch_corr_dada2, "corncob_anemia_v_treatment_batch_corr_dada2.json")
# write(corn_analysis$json$dex_muscle_v_treatment_batch_corr_dada2, "corncob_dex_muscle_v_treatment_batch_corr_dada2.json")

plot(corncob)



meta_color <- subset(x = metadata_analysis$metadata, subset = metadata_analysis$metadata$SampleID != 'F5')

cornc_df <- as.data.frame(corn_analysis[["datasets"]][["batch_corr_dada2"]]@otu_table)

dir.create('corncob/diff_abund')
dir.create('corncob/diff_var')

##############################
### DIFFERENTIAL ABUNDANCE ###

cornc_anals <- list(
  'dex_muscle_v_treatment _ Species' = corn_analysis[["analysis_taxa"]][["~ dex_muscle_v_treatment"]][["batch_corr_dada2"]][["Species"]],
  'anemia_v_treatment _ Family' = corn_analysis[["analysis_taxa"]][["~ anemia_v_treatment"]][["batch_corr_dada2"]][["Family"]],
  'anemia_v_treatment _ Species' = corn_analysis[["analysis_taxa"]][["~ anemia_v_treatment"]][["batch_corr_dada2"]][["Species"]],
  'group _ Phylum' = corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Phylum"]],
  'group _ Family' = corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Family"]],
  'group _ Genus' = corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Genus"]],
  'group _ Species' = corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Species"]]
  )

### !!! check this out

plot(corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Species"]])
plot(corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Species"]])

corn_analysis$diff_dfs <- purrr::map2(
  .x = cornc_anals, 
  .y = names(cornc_anals), 
  .f = function(cornc_anal, cornc_comp){
    
    dir.create(paste0('corncob/diff_abund/', cornc_comp))
    
    
    
    purrr::map2_df(
      .x = cornc_anal[["diff_tax_abundance"]],
      .y = names(cornc_anal[["diff_tax_abundance"]]),
      .f = function(taxons, taxon_seqs){
        
        taxon_counts <- subset(x = cornc_df, subset = rownames(cornc_df) == taxon_seqs)
        
        cornc_ggplot <- as.data.frame(t(taxon_counts))
        cornc_ggplot$names <- rownames(cornc_ggplot)
        cornc_ggplot$counts <- cornc_ggplot[[1]]
        cornc_ggplot[[1]] <- NULL
        
        ggplot(cornc_ggplot, aes(y = counts, x = names, colour = meta_color$group)) + 
          geom_point() +
          xlab(taxons) +
          labs(title = cornc_comp)
        
        ggsave(filename = paste0('corncob/diff_abund/', cornc_comp, '/', taxons, '.png'), width = 9.6, height = 5.7)
        
        taxon_counts$tax_name <- taxons
        taxon_counts$comp <- cornc_comp
        
        return(taxon_counts)
      })
  })

### DIFFERENTIAL ABUNDANCE ###
##############################








##########################
### VARIABLE ABUNDANCE ###

cornc_var_anals <- list(
  'dex_muscle_v_treatment _ Order' = corn_analysis[["analysis_taxa"]][["~ dex_muscle_v_treatment"]][["batch_corr_dada2"]][["Order"]],
  'dex_muscle_v_treatment _ Species' = corn_analysis[["analysis_taxa"]][["~ dex_muscle_v_treatment"]][["batch_corr_dada2"]][["Species"]],
  'anemia_v_treatment _ Species' = corn_analysis[["analysis_taxa"]][["~ anemia_v_treatment"]][["batch_corr_dada2"]][["Species"]],
  'group _ Phylum' = corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Phylum"]],
  'group _ Order' = corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Order"]],
  'group _ Family' = corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Family"]],
  'group _ Genus' = corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Genus"]],
  'group _ Species' = corn_analysis[["analysis_taxa"]][["~ group"]][["batch_corr_dada2"]][["Species"]]
)

corn_analysis$diff_var_dfs <- purrr::map2(
  .x = cornc_var_anals, 
  .y = names(cornc_var_anals), 
  .f = function(cornc_anal, cornc_comp){
    
    dir.create(paste0('corncob/diff_var/', cornc_comp))
    
    
    
    purrr::map2_df(
      .x = cornc_anal[["diff_tax_dispersion"]],
      .y = names(cornc_anal[["diff_tax_dispersion"]]),
      .f = function(taxons, taxon_seqs){
        
        taxon_counts <- subset(x = cornc_df, subset = rownames(cornc_df) == taxon_seqs)
        
        cornc_ggplot <- as.data.frame(t(taxon_counts))
        cornc_ggplot$names <- rownames(cornc_ggplot)
        cornc_ggplot$counts <- cornc_ggplot[[1]]
        cornc_ggplot[[1]] <- NULL
        
        ggplot(cornc_ggplot, aes(y = counts, x = names, colour = meta_color$group)) + 
          geom_point() +
          xlab(taxons) +
          labs(title = cornc_comp)
        
        ggsave(filename = paste0('corncob/diff_var/', cornc_comp, '/', taxons, '.png'), width = 9.6, height = 5.7)
        
        taxon_counts$tax_name <- taxons
        taxon_counts$comp <- cornc_comp
        
        return(taxon_counts)
      })
  })

### VARIABLE ABUNDANCE ###
##########################



























# 
# 
# corn_analysis$analysis_taxa <- purrr::map2(
#   .x = corn_analysis$formulas,
#   .y = names(corn_analysis$formulas),
#   .f = function(formula_, formula_name) {
#     
#     
#     purrr::map2(
#       .x = corn_analysis$datasets_taxa,
#       .y = names(corn_analysis$datasets_taxa),
#       .f = function(dataset, dataset_name) {
#         
#         
#         purrr::map2(
#           .x = dataset,
#           .y = names(dataset),
#           .f = function(tax_lev, tax_lev_name) {
#             
#             
#             formula_proper_format <- as.formula(formula_name)
#             
#             
#             
#             if (tax_lev_name == 'dex_muscle_v_treatment') {
#               tax_lev = subset_samples(tax_lev, !is.na(dex_muscle_v_treatment))
#             }
#             
#             
#             ### !!! is this good?
#             diff_obj_abundance <- corncob::differentialTest(formula = formula_proper_format,
#                                                             phi.formula = formula_proper_format,
#                                                             formula_null = ~ 1,
#                                                             phi.formula_null = formula_proper_format,
#                                                             test = "Wald", 
#                                                             boot = FALSE,
#                                                             data = tax_lev,
#                                                             fdr_cutoff = 0.05)
#             
#             if (length(diff_obj_abundance$significant_taxa) != 0) {
#               diff_tax_abundance <- corncob::otu_to_taxonomy(OTU = diff_obj_abundance$significant_taxa, data = tax_lev)
#             } else {diff_tax_abundance <- NA}
#             
#             
#             ### !!! is this good?
#             diff_obj_dispersion <- corncob::differentialTest(formula = formula_proper_format,
#                                                              phi.formula = formula_proper_format,
#                                                              formula_null = formula_proper_format,
#                                                              phi.formula_null = ~ 1,
#                                                              test = "LRT", 
#                                                              boot = FALSE,
#                                                              data = tax_lev,
#                                                              fdr_cutoff = 0.05)
#             
#             if (length(diff_obj_dispersion$significant_taxa) != 0) {
#               diff_tax_dispersion <- corncob::otu_to_taxonomy(OTU = diff_obj_dispersion$significant_taxa, data = tax_lev)
#             } else {diff_tax_dispersion <- NA}
#             
#             
#             ### !!! ("groupfe_dextran_muscle", "groupsiderai_fe_liposomal_oral", "groupnanofe_dextran_50nm_oral", "groupnanofe_phospholipid_70nm_oral", "groupnanofe_synomag_30nm_oral", "groupfeso4_oral") is a formula for contrasts, is this good?
#             contrasts_da <- corncob::contrastsTest(formula = formula_proper_format,
#                                                    phi.formula = formula_proper_format,
#                                                    contrasts_DA = formula_,
#                                                    data = tax_lev,
#                                                    fdr_cutoff = 0.99)
#             
#             contrasts_dv <- corncob::contrastsTest(formula = formula_proper_format,
#                                                    phi.formula = formula_proper_format,
#                                                    contrasts_DV = formula_,
#                                                    data = tax_lev,
#                                                    fdr_cutoff = 0.99)
#             
#             return(list('diff_obj_abundance' = diff_obj_abundance, 'diff_tax_abundance' = diff_tax_abundance, 'diff_obj_dispersion' = diff_obj_dispersion, 'diff_tax_dispersion' = diff_tax_dispersion, 'contrasts_da' = contrasts_da, 'contrasts_dv' = contrasts_dv))
#           })
#       })
#   })












# corn_analysis$analysis <- purrr::map2(
#   .x = corn_analysis$formulas,
#   .y = names(corn_analysis$formulas),
#   .f = function(formula_, formula_name) {
#     
#     library(corncob)
#     library(phyloseq)
#     
#     
#     purrr::map2(
#       .x = corn_analysis$datasets,
#       .y = names(corn_analysis$datasets),
#       .f = function(dataset, dataset_name) {
#         
#         
#         purrr::map2(
#           .x = dataset,
#           .y = names(dataset),
#           .f = function(tax_lev, tax_lev_name) {
#             
#             
#             formula_proper_format <- as.formula(formula_name)
#             
#             
#             
#             if (tax_lev_name == 'dex_muscle_v_treatment') {
#               tax_lev = phyloseq::subset_samples(tax_lev, !is.na(dex_muscle_v_treatment))
#             }
#             
#             
#             
#             diff_obj_abundance <- corncob::differentialTest(formula = formula_proper_format,
#                                                             phi.formula = formula_proper_format,
#                                                             formula_null = ~ 1,
#                                                             phi.formula_null = formula_proper_format,
#                                                             test = "Wald", 
#                                                             boot = T, 
#                                                             B = 1,
#                                                             data = tax_lev,
#                                                             fdr_cutoff = 0.05)
#             
#             if (length(diff_obj_abundance$significant_taxa) != 0) {
#               diff_tax_abundance <- corncob::otu_to_taxonomy(OTU = diff_obj_abundance$significant_taxa, data = tax_lev)
#             } else {diff_tax_abundance <- NA}
#             
#             
#             
#             diff_obj_dispersion <- corncob::differentialTest(formula = formula_proper_format,
#                                                              phi.formula = formula_proper_format,
#                                                              formula_null = formula_proper_format,
#                                                              phi.formula_null = ~ 1,
#                                                              test = "LRT", 
#                                                              boot = T,
#                                                              B = 1,
#                                                              data = tax_lev,
#                                                              fdr_cutoff = 0.05)
#             
#             if (length(diff_obj_dispersion$significant_taxa) != 0) {
#               diff_tax_dispersion <- corncob::otu_to_taxonomy(OTU = diff_obj_dispersion$significant_taxa, data = tax_lev)
#             } else {diff_tax_dispersion <- NA}
#             
#             
#             ### !!! ("groupfe_dextran_muscle", "groupsiderai_fe_liposomal_oral", "groupnanofe_dextran_50nm_oral", "groupnanofe_phospholipid_70nm_oral", "groupnanofe_synomag_30nm_oral", "groupfeso4_oral") is a formula for contrasts, is this good?
#             contrasts_da <- corncob::contrastsTest(formula = formula_proper_format,
#                                                    phi.formula = formula_proper_format,
#                                                    contrasts_DA = formula_,
#                                                    data = tax_lev,
#                                                    fdr_cutoff = 0.99)
#             
#             contrasts_dv <- corncob::contrastsTest(formula = formula_proper_format,
#                                                    phi.formula = formula_proper_format,
#                                                    contrasts_DV = formula_,
#                                                    data = tax_lev,
#                                                    fdr_cutoff = 0.99)
#             
#             return(list('diff_obj_abundance' = diff_obj_abundance, 'diff_tax_abundance' = diff_tax_abundance, 'diff_obj_dispersion' = diff_obj_dispersion, 'diff_tax_dispersion' = diff_tax_dispersion, 'contrasts_da' = contrasts_da, 'contrasts_dv' = contrasts_dv))
#           })
#       })
#   })