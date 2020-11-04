# devtools::install_github("jbisanz/qiime2R")
# devtools::install_github("bryandmartin/corncob")
# BiocManager::install("DESeq2")

##############################
### PREPARE DATA FOR DESEQ ###

library(Hmisc)
source('functions.R')

### 'formulas' element has actual formula in the name and contrasts as list elements. Variable of interest needs to be the last element of the formula
### !!! we should probably also analyze taxonomic tables?
DESeq_analysis <- list('freq_dir' = '/home/adrian/Desktop/nx_temp/nx_data_q/',
                       'freq_asv_patt' = '^filtered.*FeatureTable.*',
                       'freq_tax_patt' = '^lev_.*FeatureTable.*',
                       'tax_dir' = '/home/adrian/Desktop/nx_temp/nx_taxonomy/',
                       'tax_patt' = '^vsearch.*Taxonomy\\.qza$',
                       'metadata' = qiime2R::read_q2metadata(file = '../qiime_metadata_nf.tsv'),
                       'formulas' = 
                         list('~ group' = list(
                           'f_d_m' = c('group', 'fe_dextran_muscle', 'anemia'),
                           's_f_l_o' = c('group', 'siderai_fe_liposomal_oral', 'anemia'))),
                       'seq_col_name' = 'seqs')



### PREPARE MERGED TAXONOMY
DESeq_analysis$tax_files = list.files(path = DESeq_analysis$tax_dir, pattern = DESeq_analysis$tax_patt)

DESeq_analysis$tax_merged <- prepare_merged_taxonomies(list_tax_files = DESeq_analysis$tax_files, str_tax_dir = DESeq_analysis$tax_dir, seq_col_name_ = DESeq_analysis$seq_col_name)



### PREPARE METADATA
rownames(DESeq_analysis$metadata) <- as.character(DESeq_analysis$metadata[[1]])

### All categorical variables should be factors I think! Control group needs to be the first factor. This can be set using relevel ###
DESeq_analysis$metadata$group <- relevel(DESeq_analysis$metadata$group, ref = "anemia")



### READ DATA FILES
DESeq_analysis$freq_files = list.files(path = DESeq_analysis$freq_dir, pattern = paste0(DESeq_analysis$freq_asv_patt, '|', DESeq_analysis$freq_tax_patt))

DESeq_analysis$freq <- purrr::map(
  .x = DESeq_analysis$freq_files, 
  .f = function(x){
    
    qiime2R::read_qza(paste0(DESeq_analysis$freq_dir, x))
  })
names(DESeq_analysis$freq) <- stringr::str_remove(string = DESeq_analysis$freq_files, pattern = '_FeatureTable.*')






DESeq_analysis$freq_reorder <- purrr::map(
  .x = DESeq_analysis$freq, 
  .f = function(freq_data){
    
    as.data.frame(freq_data$data[, rownames(DESeq_analysis$metadata)])
  })

### PREPARE DATA FOR DESEQ ###
##############################





#############
### DESEQ ###

DESeq_analysis$deseq_data <- purrr::map(
  .x = DESeq_analysis$freq_reorder,
  .f = function(freq_data){

    
    
    temp <- purrr::map(
      .x = names(DESeq_analysis$formulas),
      .f = function(formula_){
        
        if (all(colnames(freq_data) == rownames(DESeq_analysis$metadata))) {
          DESeq2::DESeqDataSetFromMatrix(
            countData = freq_data,
            colData = DESeq_analysis$metadata,
            design = as.formula(formula_) )
        }
      })
    names(temp) <- as.character(names(DESeq_analysis$formulas))

    return(temp)
  })






DESeq_analysis$analysis <- purrr::map2(
  .x = DESeq_analysis$deseq_data,
  .y = names(DESeq_analysis$deseq_data),
  .f = function(dataset, dataset_name) {
    
    dataset_name <- stringr::str_remove(string = dataset_name, pattern = '.FeatureTableFrequency.qza$')
    
    
    
    purrr::pmap(
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
  })

### DESEQ ###
#############






#############################
### ASV CENTERED ANALYSIS ###

DESeq_analysis$asv_centered$analysis <- DESeq_analysis$analysis[stringr::str_detect(string = names(DESeq_analysis$analysis), pattern = DESeq_analysis$freq_asv_patt)]



DESeq_analysis$asv_centered$merged_group_data <- purrr::map2(
  .x = DESeq_analysis$formulas,
  .y = names(DESeq_analysis$formulas),
  .f = function(formula, formula_name){


    formulas_<- purrr::map2(
      .x = formula,
      .y = names(formula),
      .f = function(contrast, contrast_names){


        contrasts_ <- purrr::map(
          .x = DESeq_analysis$asv_centered$analysis,
          .f = function(analysis){

            contrast_ <- analysis[[formula_name]]$contrasts[[contrast_names]]$results_df

            contrast_ <- subset(x = contrast_, subset = contrast_$pvalue < 0.05 )
          })

        contrasts_ <- rlist::list.rbind(contrasts_)

        contrasts_ <- wrapper_get_all_values_for_seqences(
          list_of_original_datasets = DESeq_analysis$freq_reorder,
          df_with_sequences = contrasts_,
          df_with_sequences_sequence_colname = 'seqs')
        contrasts_ <- unique(dplyr::select(contrasts_, -dataset))

        return(contrasts_)
      })
  })



DESeq_analysis$asv_centered$taxonomic_data_for_visualization <- purrr::map2(
  .x = DESeq_analysis$asv_centered$merged_group_data,
  .y = names(DESeq_analysis$asv_centered$merged_group_data),
  .f = function(variable, variable_name){
    

    purrr::map2(
      .x = variable,
      .y = names(variable),
      .f = function(contrast_, contrast_name){
        
        contrast_$full_name <- NA
        
        for (row in seq_along(contrast_[[1]])) {
          contrast_$full_name[[row]] <- paste(contrast_[row,8:14], collapse = '__')
        }
        
        return( dplyr::select(contrast_, 1, 18:27) )
      })    
  }) 

### ASV CENTERED ANALYSIS ###
#############################







#############################
### TAX CENTERED ANALYSIS ###

DESeq_analysis$tax_centered <- purrr::map(
  .x = seq(2, 7),
  .f = function(level){
    
    tax_centered_analysis <- DESeq_analysis$analysis[stringr::str_detect(string = names(DESeq_analysis$analysis), pattern = DESeq_analysis$freq_tax_patt)]
    
    tax_centered_analysis[stringr::str_detect(
      string = names(tax_centered_analysis), 
      pattern = paste0('^lev_', as.character(level)))]
    
  })
names(DESeq_analysis$tax_centered) <- as.character(seq(2, 7))






DESeq_analysis$tax_centered$merged_group_data <- purrr::map2(
  .x = DESeq_analysis$formulas,
  .y = names(DESeq_analysis$formulas),
  .f = function(formula, formula_name){

    
    purrr::map2(
      .x = formula,
      .y = names(formula),
      .f = function(contrast, contrast_names){
       
         
        purrr::map2(
          .x = DESeq_analysis$tax_centered,
          .y = names(DESeq_analysis$tax_centered),
          .f = function(tax_lev, tax_lev_name){ 
            
            
            
            contrasts_ <- purrr::map(
              .x = tax_lev,
              .f = function(analysis){
                
                contrast_ <- analysis[[formula_name]]$contrasts[[contrast_names]]$results_df
                
                contrast_ <- subset(x = contrast_, subset = contrast_$pvalue < 0.05 )
              })
            
            contrasts_ <- rlist::list.rbind(contrasts_)
            
            if (!is.null(contrasts_)) {
              contrasts_ <- contrasts_[,-c(8:15)]
            }
            
            return(contrasts_)
    })
  })
})

### TAX CENTERED ANALYSIS ###
#############################









# plotCounts - for extracting single genes for visualization (ggplot) i think?