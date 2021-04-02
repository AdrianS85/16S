source('https://raw.githubusercontent.com/AdrianS85/helper_R_functions/master/little_helpers.R')

prepare_merged_taxonomies <- function(list_tax_files, str_tax_dir, seq_col_name_)
{
  taxonomy <- purrr::map(
    .x = list_tax_files, 
    .f = function(file){
      file_tax <- qiime2R::read_qza(file = paste0(str_tax_dir, file))
      
      file_tax_df <- qiime2R::parse_taxonomy(file_tax$data)
      
      file_tax_df[[seq_col_name_]] <- rownames(file_tax_df)
      
      return(file_tax_df)
    })
  
  taxonomy <- unique(rlist::list.rbind(taxonomy))

  return(taxonomy)
}






wrapper_get_all_values_for_seqences <- function(list_of_original_datasets, df_with_sequences, df_with_sequences_sequence_colname, collapse_values_from_dataset_using_median = T)
{
  assertthat::assert_that(are_vectors_the_same(purrr::map(.x = list_of_original_datasets, .f = colnames))) ### !!! check this

  list_of_original_datasets_with_names <- purrr::map2(
    .x = list_of_original_datasets,
    .y = names(list_of_original_datasets),
    .f = function(dataset_, dataset_name) {
      
      if (collapse_values_from_dataset_using_median == F) {
        dataset_$dataset_of_origin <- stringr::str_remove(dataset_name, '\\..*')
      }
      
      dataset_[[df_with_sequences_sequence_colname]] <- stringr::str_remove_all(rownames(dataset_), '.*\\.') ### !!! ale dlaczego tutaj kropkę usuwamy. ja nie rozumie tego tutaj

      return(dataset_)
    })
  
  df_with_values_for_all_seqs_from_all_datasets <- rlist::list.rbind(list_of_original_datasets_with_names)
  
  if (collapse_values_from_dataset_using_median == T) {
    warning('BEWARE! Values from all datasets in list_of_original_datasets are medianed! This is ok for e.g. filtered data for deblur, deblur with singletons and dada2, but that may not be so every time')
    
    df_with_values_for_all_seqs_from_all_datasets <- group_and_then_summarize_all_cols_in_df(
      df_ = df_with_values_for_all_seqs_from_all_datasets, 
      group_by_col = df_with_sequences_sequence_colname, 
      summarize_method = 'median')
    }

  df_with_values_for_all_seqs_from_all_datasets <- merge(
    x = df_with_sequences, 
    y = df_with_values_for_all_seqs_from_all_datasets,
    by = df_with_sequences_sequence_colname,
    all.x = T)

  return(df_with_values_for_all_seqs_from_all_datasets)
}










calculations_of_metadata_plots <- function(prepared_metadata, dataset_name, full_metadata, full_metadata_groupid_col, full_metadata_sampleid_col = 'SampleID'){

  prepared_metadata <- prepared_metadata[ order(rownames(prepared_metadata)), ]
  
  
  
  habillage <- subset(
    x = full_metadata, 
    subset = full_metadata[[full_metadata_sampleid_col]] %in% rownames(prepared_metadata),
    select = c(full_metadata_sampleid_col, full_metadata_groupid_col))
  
  habillage[[full_metadata_sampleid_col]] <- as.character(habillage[[full_metadata_sampleid_col]])
  habillage[[full_metadata_groupid_col]] <- as.character(habillage[[full_metadata_groupid_col]])
  
  habillage <- habillage[order(habillage[[full_metadata_sampleid_col]]), ]
  
  assertthat::assert_that(are_vectors_the_same(list(habillage[[full_metadata_sampleid_col]], rownames(prepared_metadata))))
  
  habillage <- habillage[[full_metadata_groupid_col]]
  
  
  
  cor <- Hmisc::rcorr(as.matrix(prepared_metadata))
  pca_scaled <- prcomp(t(prepared_metadata), scale = F) ### !!! scale????
  pca_scaled_samples <- prcomp(prepared_metadata, scale = F) ### !!! scale????
  

  
  dir <- paste0(metadata_analysis$dir, '/', dataset_name, '/')
  
  dir.create(dir)
  
  return(list('data' = prepared_metadata, 'habillage' = habillage, 'cor' = cor, 'pca_scaled' = pca_scaled, 'pca_scaled_samples' = pca_scaled_samples, 'dir' = dir))
}













merge_measures <- function(prepared_metadata, cols_to_merge_char_vec, merged_col_name, scale_before_merging = F)
{
  to_merge <- subset(x = prepared_metadata, select = cols_to_merge_char_vec)
  
  if (scale_before_merging == T) {
    
    to_merge <- purrr::map_df(
      .x = to_merge, 
      .f = function(x){scale(to_merge, center = TRUE, scale = TRUE)})
  }

  prepared_metadata[[merged_col_name]] <- apply(X = to_merge, MARGIN = 1, FUN = median)
  
  prepared_metadata <- dplyr::select(prepared_metadata, -cols_to_merge_char_vec)

  return(prepared_metadata)
}













### !!! to do
draw_log_ratios <- function(taxon, fitted_taxon_, observed_taxon_) {
  data.frame(temperature = temp_continuous,
             log_ratio = fitted_taxon_[,taxon],
             log_ratio_observed = observed_taxon_[,taxon]) %>%
    ggplot(aes(x = temperature, y = log_ratio_observed)) +
    geom_point() +
    xlab("Temperature") +
    ylab("log-ratios") +
    geom_point(aes(y = log_ratio), col = "blue") +
    ggtitle("Observed (black) and fitted (blue) log-ratios against temperature")
}



















draw_alpha <- function(
  dir_ = divnet_analysis_2$vis_alpha_dir, 
  analysis_names_ = col_name,
  dataset_name_ = dataset_name,
  metric_,
  divnet_ = col,
  original_dataset_ = original_dataset)
{
  library(ggplot2)

  
  
  purrr::walk2(
    .x = divnet_, 
    .y = names(divnet_),
    .f = function(column, column_name){
      

      
      
      purrr::walk(
        .x = analysis_names_,
        .f = function(col_name){
          
          dir.create(paste0(dir_, '/', dataset_name_))
          
          plot_name <- paste0(dir_, '/', dataset_name_, '/', col_name, '__',   metric_, '__', column_name, '__', dataset_name_, '.png')
      
          temp <- summary(column[[metric_]])
          temp[[column_name]] <- recode_values_based_on_key(
            to_recode_chrvec = temp$sample_names, 
            replace_this_chrvec = rownames(original_dataset_@sam_data), 
            with_this_chrvec = original_dataset_@sam_data[[column_name]])
          
          temp <- temp[order(temp[[column_name]]),] 
          temp$sample_names_factor <- factor(x = temp$sample_names, levels = temp$sample_names)
    
          main_plot <- ggplot(data = temp, mapping = aes(x = sample_names_factor, y = estimate, color = !!as.symbol(column_name) )) +
            geom_point() +
            geom_errorbar(aes(ymin = lower, ymax = upper), width=.1)
          
          ggsave(
            filename = plot_name, 
            plot = main_plot,
            width = 8.5, 
            height = 5.4)
        
        })
  
  })
}




mean_sd_percent_by_rows_phyloseq <- function(phyloseq_object, taxon_of_the_object)
{
  library(dplyr)
  
  options(scipen = 999)
  
  msp_rowsum <- sum(rowSums(phyloseq::otu_table(phyloseq_object)))
  
  msp <- data.frame(phyloseq::otu_table(phyloseq_object))
  
  msp_rowsum_vector <- rep(msp_rowsum, times = length(msp[[1]]))

  
  
  msp$taxon <- recode_values_based_on_key(
    to_recode_chrvec = rownames(msp), 
    replace_this_chrvec = rownames(data.frame(phyloseq::tax_table(phyloseq_object))), 
    with_this_chrvec = data.frame(phyloseq::tax_table(phyloseq_object))[[taxon_of_the_object]]
  ) 
  
  last_col_ <- length(msp)
  
  msp <- select(.data = msp, taxon, everything()) %>%
    rowwise() %>%
    mutate(mean = mean(c_across(cols = 2:last_col_))) %>%
    mutate(sd = sd(c_across(cols = 2:last_col_))) %>%
    mutate(total = sum(c_across(cols = 2:last_col_))) %>%
    select(taxon, mean, sd, total, everything())
  
  msp$percentage <- msp$total/msp_rowsum_vector*100
  
  msp <- select(.data = msp, taxon, mean, sd, total, percentage, everything())
  
  return(msp)
}




get_nice_results_from_lrt_deseq2 <- function(results_lrt, phyloseq_object)
{
  lrt_df <- as.data.frame(results_lrt)
  lrt_df_signif <- subset(x = lrt_df, subset = lrt_df$pvalue < 0.05)
  
  lrt_df_signif_full <- merge(
    x = lrt_df_signif, 
    y = data.frame(phyloseq_object@tax_table@.Data), 
    by = 'row.names', 
    all.x = T)
  
  rownames(lrt_df_signif_full) <- lrt_df_signif_full$Row.names
  lrt_df_signif_full$Row.names <- NULL
  
  lrt_df_signif_full <- merge(
    x = lrt_df_signif_full, 
    y = data.frame(phyloseq_object@otu_table), 
    by = 'row.names', 
    all.x = T)
  
  lrt_df_signif_full_tidy <- tidyr::pivot_longer(
    data = lrt_df_signif_full,
    cols = c(16:47), 
    names_to = 'sample', 
    values_to = 'frequency')
  
  lrt_df_signif_full_tidy <- dplyr::select(
    .data = lrt_df_signif_full_tidy,
    Row.names, log2FoldChange, pvalue, padj, Kingdom:Species, sample, frequency)
  
  lrt_df_signif_full_tidy_metadata <- merge(
    x = lrt_df_signif_full_tidy, 
    y = data.frame(phyloseq_object@sam_data), 
    by.x = 'sample', 
    by.y = 0, 
    all.x = T)
  
  for (row_nb in seq_along(rownames(lrt_df_signif_full_tidy_metadata))) {
    
    lrt_df_signif_full_tidy_metadata$full_name[[row_nb]] <- paste(
      lrt_df_signif_full_tidy_metadata[row_nb,6:12],
      collapse = '__')
  }
  
  lrt_df_signif_full_tidy_metadata$full_name <- stringr::str_remove(string = lrt_df_signif_full_tidy_metadata$full_name, pattern = '^d__')
  
  # openxlsx::write.xlsx(
  #   x = subset(
  #     x = lrt_df_signif_full, 
  #     subset = lrt_df_signif_full$pvalue <0.05),
  #   file = 'diff_abund_significant_asvs.xlsx')
  
  return(list('lrt_df_signif_full' = lrt_df_signif_full, 'lrt_df_signif_full_tidy_metadata' = lrt_df_signif_full_tidy_metadata))
  
}