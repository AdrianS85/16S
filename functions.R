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