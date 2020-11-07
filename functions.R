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
      
      dataset_[[df_with_sequences_sequence_colname]] <- stringr::str_remove_all(rownames(dataset_), '.*\\.')

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
  
  prepared_metadata <- as.data.frame(scale(prepared_metadata, center = TRUE, scale = TRUE))
  
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
  pca_scaled <- prcomp(t(prepared_metadata), scale = TRUE)
  pca_scaled_samples <- prcomp(prepared_metadata, scale = TRUE)
  

  
  dir <- paste0(metadata_analysis$dir, '/', dataset_name, '/')
  
  dir.create(dir)
  
  return(list('data' = prepared_metadata, 'habillage' = habillage, 'cor' = cor, 'pca_scaled' = pca_scaled, 'pca_scaled_samples' = pca_scaled_samples, 'dir' = dir))
}













merge_measures <- function(prepared_metadata, cols_to_merge_char_vec, merged_col_name)
{
  to_merge <- subset(x = prepared_metadata, select = cols_to_merge_char_vec)
  
  prepared_metadata[[merged_col_name]] <- apply(X = to_merge, MARGIN = 1, FUN = median)
  
  prepared_metadata <- dplyr::select(prepared_metadata, -cols_to_merge_char_vec)

  return(prepared_metadata)
}













