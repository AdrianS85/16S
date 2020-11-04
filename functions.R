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


