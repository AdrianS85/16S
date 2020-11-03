devtools::install_github("jbisanz/qiime2R")
devtools::install_github("bryandmartin/corncob")
BiocManager::install("DESeq2")
library(Hmisc)
# read_qza() - Function for reading artifacts (.qza).
# qza_to_phyloseq() - Imports multiple artifacts to produce a phyloseq object.
# read_q2metadata() - Reads qiime2 metadata file (containing q2-types definition line)
# write_q2manifest() - Writes a read manifest file to import data into qiime2
# theme_q2r() - A ggplot2 theme for publication-type figures.
# print_provenance() - A function to display provenance information.
# is_q2metadata() - A function to check if a file is a qiime2 metadata file.
# parse_taxonomy() - A function to parse taxonomy strings and return a table where each column is a taxonomic class.
# parse_ordination() - A function to parse the internal ordination format.
# read_q2biom - A function for reading QIIME2 biom files in format v2.1
# make_clr - Transform feature table using centered log2 ratio.
# make_proportion - Transform feature table to proportion (sum to 1).
# make_percent - Transform feature to percent (sum to 100).
# interactive_table - Create an interactive table in Rstudio viewer or rmarkdown html.
# summarize_taxa- Create a list of tables with abundances sumed to each taxonomic level.
# taxa_barplot - Create a stacked barplot using ggplot2.
# taxa_heatmap - Create a heatmap of taxonomic abundances using gplot2.

### !!! we should probably also analyze taxonomic tables?
DESeq_analysis <- list('freq_dir' = '/home/adrian/Desktop/nx_temp/nx_data_q/',
                       'freq_patt' = '^filtered.*FeatureTable.*',
                       'tax_dir' = '/home/adrian/Desktop/nx_temp/nx_taxonomy/',
                       'tax_patt' = '^vsearch.*Taxonomy\\.qza$',
                       'metadata' = qiime2R::read_q2metadata(file = '../qiime_metadata_nf.tsv'),
                       'formulas' = 
                         list('~ group' = list(
                           'f_d_m' = c('group', 'fe_dextran_muscle', 'anemia'),
                           's_f_l_o' = c('group', 'siderai_fe_liposomal_oral', 'anemia'))),
                       'seq_col_name' = 'seqs')

DESeq_analysis$freq_files = list.files(path = DESeq_analysis$freq_dir, pattern = DESeq_analysis$freq_patt)
DESeq_analysis$tax_files = list.files(path = DESeq_analysis$tax_dir, pattern = DESeq_analysis$tax_patt)

DESeq_analysis$tax_merged <- prepare_merged_taxonomies(list_tax_files = DESeq_analysis$tax_files, str_tax_dir = DESeq_analysis$tax_dir, seq_col_name_ = DESeq_analysis$seq_col_name)

rownames(DESeq_analysis$metadata) <- as.character(DESeq_analysis$metadata[[1]])





DESeq_analysis$freq <- purrr::map(
  .x = DESeq_analysis$freq_files, 
  .f = function(x){
    qiime2R::read_qza(paste0(DESeq_analysis$freq_dir, x))
  })
names(DESeq_analysis$freq) <- stringr::str_remove(string = DESeq_analysis$freq_files, pattern = '_FeatureTable.*')



DESeq_analysis$freq_reorder <- purrr::map(
  .x = DESeq_analysis$freq, 
  .f = function(freq_data){
    new_order <- as.data.frame(freq_data$data[, rownames(DESeq_analysis$metadata)])
  })
names(DESeq_analysis$freq_reorder) <- names(DESeq_analysis$freq)





DESeq_analysis$deseq_data <- purrr::map(
  .x = DESeq_analysis$freq_reorder,
  .f = function(freq_data){

    temp <- purrr::map(
      .x = names(DESeq_analysis$formulas),
      .f = function(formula_){
        
        if (all(colnames(freq_data) == rownames(DESeq_analysis$metadata))) {
          formula_dataset <- DESeq2::DESeqDataSetFromMatrix(
            countData = freq_data,
            colData = DESeq_analysis$metadata,
            design = as.formula(formula_) )

        }
      })
    names(temp) <- as.character(names(DESeq_analysis$formulas))

    return(temp)
  })
names(DESeq_analysis$deseq_data) <- names(DESeq_analysis$freq)









DESeq_analysis$analysis <- purrr::map2(
  .x = DESeq_analysis$deseq_data,
  .y = names(DESeq_analysis$deseq_data),
  .f = function(dataset, dataset_name) {
    
    dataset_name <- stringr::str_remove(string = dataset_name, pattern = '.FeatureTableFrequency.qza$')
    
    
    
    formulas_ <- purrr::pmap(
      .l = list(dataset, names(dataset), DESeq_analysis$formulas, names(DESeq_analysis$formulas)),
      .f = function(formula_dataset, formula_dataset_name, formula_expected, formula_expected_name) {
        DESeq_ <- DESeq2::DESeq(formula_dataset)
        assertthat::assert_that(formula_dataset_name == formula_expected_name)
        
        
        
        contrasts <- purrr::map2(
          .x = formula_expected,
          .y = names(formula_expected),
          .f = function(contrast, contrast_name){
            results <- DESeq2::results(DESeq_, contrast = contrast)
            
            results_df <- as.data.frame(results)
            results_df$seqs <- rownames(results_df)
            results_df <- merge(results_df, DESeq_analysis$tax_merged, by = DESeq_analysis$seq_col_name, all.x = T)
            results_df$formula <- formula_expected_name
            results_df$dataset <- dataset_name
            results_df$contrast <- contrast_name
            
            return(list('results' = results, 'results_df' = results_df))
          })
        
        return(list('DESeq' = DESeq_, 'contrasts' = contrasts))
      })
    
    names(formulas_) <- as.character(names(DESeq_analysis$formulas))
    
    return(formulas_)
  })










DESeq_analysis$merged_group_data <- purrr::map2(
  .x = DESeq_analysis$formulas,
  .y = names(DESeq_analysis$formulas),
  .f = function(formula, formula_name){
    
    
    
    formulas_<- purrr::map2(
      .x = formula,
      .y = names(formula),
      .f = function(contrast, contrast_names){
        
        
        
        contrasts_ <- purrr::map(
          .x = DESeq_analysis$analysis,
          .f = function(analysis){
            
            contrast_ <- analysis[[formula_name]]$contrasts[[contrast_names]]$results_df
            
            contrast_ <- subset(x = contrast_, subset = contrast_$pvalue < 0.05 )
          })
        
        contrasts_ <- rlist::list.rbind(contrasts_)
    })
})    








DESeq_analysis$tax_merged_names <- tidyr::unite(data = DESeq_analysis$tax_merged, 1:7, col = name, sep = '_X_')

DESeq_analysis$taxonomic_data_for_visualization <- purrr::map2(
  .x = DESeq_analysis$freq_reorder,
  .y = names(DESeq_analysis$freq_reorder),
  .f = function(dataset, dataset_name){
    

})    










