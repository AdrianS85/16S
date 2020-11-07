library(magrittr)
library(ggplot2)

metadata_analysis <- list('metadata' = qiime2R::read_q2metadata(file = '../qiime_metadata.tsv'),
                          'dir' = 'metadata',
                          'distribution' = 'metadata/dist'
                          )

metadata_analysis$metadata_tidy <- tidyr::pivot_longer(
  data = metadata_analysis$metadata,
  cols = where(is.numeric),
  names_to = 'variable')
metadata_analysis$metadata_tidy$group <- as.character(metadata_analysis$metadata_tidy$group)

dir.create(metadata_analysis$dir)
dir.create(metadata_analysis$distribution)





##########################
### SEPARATE VARIABLES ###

purrr::walk(
  .x = unique(metadata_analysis$metadata_tidy$variable), 
  .f = function(col_name){
    
    plot_data <- metadata_analysis$metadata_tidy %>%
      dplyr::filter(variable == col_name)
    
    plot_ <- ggplot(plot_data, aes(x = value)) + geom_histogram() + geom_density(color = 'red')
    ggsave(filename = paste0(metadata_analysis$distribution, '/', col_name, '.png'), plot = plot_)
    
    boxplot_ <- ggplot(plot_data, aes(x = group, y = value, color = group)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    ggsave(filename = paste0(metadata_analysis$distribution, '/boxplot_', col_name, '.png'), plot = boxplot_)
    
    purrr::walk(
      .x = unique(plot_data$group), 
      .f = function(group_){
      
        dir_name <- paste0(metadata_analysis$distribution, '/', group_)
          
        dir.create(dir_name)
        
        plot_data_group <- plot_data %>%
          dplyr::filter(group == group_)
        
        plot_group <- ggplot(plot_data_group, aes(x = value)) + geom_histogram() + geom_density(color = 'red')
        ggsave(filename = paste0(dir_name, '/', col_name, group_, '.png'), plot = plot_group)
    })
  })

### SEPARATE VARIABLES ###
##########################






# Pearson will be good - group count is to small to see normality, but I think it should be normal, cause why not
########################
### PREPARE DATASETS ###

metadata_analysis$metadata_for_corr_all <- metadata_analysis$metadata
metadata_analysis$metadata_for_corr_all$group <- NULL
rownames(metadata_analysis$metadata_for_corr_all) <- metadata_analysis$metadata_for_corr_all$SampleID
metadata_analysis$metadata_for_corr_all$SampleID <- NULL



metadata_analysis$corr$no_nas$data <- metadata_analysis$metadata_for_corr_all[-1,-1]

metadata_analysis$corr$no_nas <- calculations_of_metadata_plots(metadata_analysis$corr$no_nas$data, dataset_name = 'no_nas', full_metadata = metadata_analysis$metadata, full_metadata_groupid_col = 'group')



metadata_analysis$corr$fecal$data <- metadata_analysis$metadata_for_corr_all[!is.na(metadata_analysis$metadata_for_corr_all$`non-heme_fe_faces_mg_fe_kg`),] %>%
  dplyr::select(-ret_k_uL, -`ret-he_pg`, -`rbc-he_pg`)

metadata_analysis$corr$fecal <- calculations_of_metadata_plots(metadata_analysis$corr$fecal$data, dataset_name = 'fecal', full_metadata = metadata_analysis$metadata, full_metadata_groupid_col = 'group')



metadata_analysis$corr$fecal_samples_but_no_fecal_data$data <- metadata_analysis$metadata_for_corr_all[!is.na(metadata_analysis$metadata_for_corr_all$`non-heme_fe_faces_mg_fe_kg`),] %>%
  dplyr::select(-ret_k_uL, -`ret-he_pg`, -`rbc-he_pg`, -`non-heme_fe_faces_mg_fe_kg`)

metadata_analysis$corr$fecal_samples_but_no_fecal_data <- calculations_of_metadata_plots(metadata_analysis$corr$fecal_samples_but_no_fecal_data$data, dataset_name = 'fecal_samples_but_no_fecal_data', full_metadata = metadata_analysis$metadata, full_metadata_groupid_col = 'group')



metadata_analysis$corr$fecal_samples_but_no_fecal_data_with_more_variables$data <- metadata_analysis$metadata_for_corr_all[!is.na(metadata_analysis$metadata_for_corr_all$`non-heme_fe_faces_mg_fe_kg`),] %>%
  dplyr::select(-`non-heme_fe_faces_mg_fe_kg`)
metadata_analysis$corr$fecal_samples_but_no_fecal_data_with_more_variables$data <- subset(x = metadata_analysis$corr$fecal_samples_but_no_fecal_data_with_more_variables$data, subset = !is.na(metadata_analysis$corr$fecal_samples_but_no_fecal_data_with_more_variables$data$ret_k_uL)) 

metadata_analysis$corr$fecal_samples_but_no_fecal_data_with_more_variables <- calculations_of_metadata_plots(metadata_analysis$corr$fecal_samples_but_no_fecal_data_with_more_variables$data, dataset_name = 'fecal_samples_but_no_fecal_data_with_more_variables', full_metadata = metadata_analysis$metadata, full_metadata_groupid_col = 'group')



metadata_analysis$corr$no_nas_rbc_morph$data <- merge_measures(
  prepared_metadata = metadata_analysis$corr$no_nas$data,
  cols_to_merge_char_vec = c('mcv_fL', 'mch_pg', 'rbc-he_pg'),
  merged_col_name = 'rbc_morphology')

metadata_analysis$corr$no_nas_rbc_morph <- calculations_of_metadata_plots(metadata_analysis$corr$no_nas_rbc_morph$data, dataset_name = 'no_nas_rbc_morph', full_metadata = metadata_analysis$metadata, full_metadata_groupid_col = 'group')



metadata_analysis$corr$fecal_rbc_morph$data <- metadata_analysis$metadata_for_corr_all[!is.na(metadata_analysis$metadata_for_corr_all$`non-heme_fe_faces_mg_fe_kg`),] %>%
  merge_measures(
    cols_to_merge_char_vec = c('mcv_fL', 'mch_pg', 'rbc-he_pg'),
    merged_col_name = 'rbc_morphology') %>%
  dplyr::select(-ret_k_uL, -`ret-he_pg`)

metadata_analysis$corr$fecal_rbc_morph$data <- subset(x = metadata_analysis$corr$fecal_rbc_morph$data, subset = !is.na(metadata_analysis$corr$fecal_rbc_morph$data[['rbc_morphology']])) 

metadata_analysis$corr$fecal_rbc_morph <- calculations_of_metadata_plots(metadata_analysis$corr$fecal_rbc_morph$data, dataset_name = 'fecal_rbc_morph', full_metadata = metadata_analysis$metadata, full_metadata_groupid_col = 'group')




metadata_analysis$corr$no_nas_rbc_all$data <- merge_measures(
  prepared_metadata = metadata_analysis$corr$no_nas$data,
  cols_to_merge_char_vec = c('mcv_fL', 'mch_pg', 'rbc-he_pg', 'hgb_g_dL', 'rbc_m_uL'),
  merged_col_name = 'rbc_all')

metadata_analysis$corr$no_nas_rbc_all <- calculations_of_metadata_plots(metadata_analysis$corr$no_nas_rbc_all$data, dataset_name = 'no_nas_rbc_all', full_metadata = metadata_analysis$metadata, full_metadata_groupid_col = 'group')



metadata_analysis$corr$fecal_rbc_all$data <- metadata_analysis$metadata_for_corr_all[!is.na(metadata_analysis$metadata_for_corr_all$`non-heme_fe_faces_mg_fe_kg`),] %>%
  merge_measures(
    cols_to_merge_char_vec = c('mcv_fL', 'mch_pg', 'rbc-he_pg', 'hgb_g_dL', 'rbc_m_uL'),
    merged_col_name = 'rbc_all') %>%
  dplyr::select(-ret_k_uL, -`ret-he_pg`)  
  
metadata_analysis$corr$fecal_rbc_all$data <- subset(x = metadata_analysis$corr$fecal_rbc_all$data, subset = !is.na(metadata_analysis$corr$fecal_rbc_all$data[['rbc_all']])) 

metadata_analysis$corr$fecal_rbc_all <- calculations_of_metadata_plots(metadata_analysis$corr$fecal_rbc_all$data, dataset_name = 'fecal_rbc_all', full_metadata = metadata_analysis$metadata, full_metadata_groupid_col = 'group')



metadata_analysis$corr$no_nas_rbc_2$data <- 
  merge_measures(
    prepared_metadata = metadata_analysis$corr$no_nas$data,
    cols_to_merge_char_vec = c('mcv_fL', 'mch_pg', 'rbc-he_pg'),
    merged_col_name = 'rbc_morphology') %>%
  merge_measures(
    cols_to_merge_char_vec = c('hgb_g_dL', 'rbc_m_uL'),
    merged_col_name = 'rbc_quant')

metadata_analysis$corr$no_nas_rbc_2 <- calculations_of_metadata_plots(metadata_analysis$corr$no_nas_rbc_2$data, dataset_name = 'no_nas_rbc_2', full_metadata = metadata_analysis$metadata, full_metadata_groupid_col = 'group')



metadata_analysis$corr$fecal_rbc_2$data <- 
  metadata_analysis$metadata_for_corr_all[!is.na(metadata_analysis$metadata_for_corr_all$`non-heme_fe_faces_mg_fe_kg`),] %>%
  merge_measures(
    cols_to_merge_char_vec = c('mcv_fL', 'mch_pg', 'rbc-he_pg'),
    merged_col_name = 'rbc_morphology') %>%
  merge_measures(
    cols_to_merge_char_vec = c('hgb_g_dL', 'rbc_m_uL'),
    merged_col_name = 'rbc_quant') %>%
  dplyr::select(-ret_k_uL, -`ret-he_pg`)  

metadata_analysis$corr$fecal_rbc_2$data <- subset(x = metadata_analysis$corr$fecal_rbc_2$data, subset = !is.na(metadata_analysis$corr$fecal_rbc_2$data[['rbc_morphology']])) 

metadata_analysis$corr$fecal_rbc_2$data <- subset(x = metadata_analysis$corr$fecal_rbc_2$data, subset = !is.na(metadata_analysis$corr$fecal_rbc_2$data[['rbc_quant']])) 

metadata_analysis$corr$fecal_rbc_2 <- calculations_of_metadata_plots(
  metadata_analysis$corr$fecal_rbc_2$data, 
  dataset_name = 'fecal_rbc_2', 
  full_metadata = metadata_analysis$metadata, 
  full_metadata_groupid_col = 'group')

### PREPARE DATASETS ###
########################






##############################
### VISUALIZE CORRELATIONS ###

purrr::walk2(
  .x = metadata_analysis$corr, 
  .y = names(metadata_analysis$corr),
  .f = function(dataset, dataset_name){
    
    cor_plot <- ggcorrplot::ggcorrplot(dataset$cor$r, 
                                       hc.order = TRUE, 
                                       type = "lower", 
                                       lab = TRUE, 
                                       p.mat = ggcorrplot::cor_pmat(dataset$cor$r))
    ggsave(filename = paste0(dataset$dir, dataset_name, '_correlation_metadata.png'), plot = cor_plot, dpi = 250)
    

    
    ggsave(filename = paste0(dataset$dir, dataset_name, '_correlation_metadata_2.png'), plot = psych::pairs.panels(dataset$data), dpi = 250)
  })



for (dataset_name in names(metadata_analysis$corr)) {
  
  z_standarized <- scale(metadata_analysis$corr[[dataset_name]]$data, center = TRUE, scale = TRUE)
  
  png(filename = paste0(metadata_analysis$corr[[dataset_name]]$dir, dataset_name, '_heatmap_metadata.png'), width = 1600, height = 900)
  gplots::heatmap.2(x = as.matrix(z_standarized), trace = 'none', margins = c(15, 5))
  dev.off()
}

### VISUALIZE CORRELATIONS ###
##############################






###########
### PCA ###

purrr::walk2(
  .x = metadata_analysis$corr, 
  .y = names(metadata_analysis$corr),
  .f = function(dataset, dataset_name){
    
    ggsave(
      filename = paste0(dataset$dir, dataset_name, '_pca_metadata.png'), 
      plot = factoextra::fviz_pca_ind(dataset$pca_scaled, repel = T), 
      dpi = 250)
    
    if (dataset_name != 'no_nas_pca_fil') {
      
      ggsave(
        filename = paste0(dataset$dir, dataset_name, '_pca_metadata_samples.png'), 
        plot = factoextra::fviz_pca_biplot(
          dataset$pca_scaled_samples, 
          repel = T, 
          habillage = dataset$habillage),
        dpi = 250)
    }
  })

### PCA ###
###########

















#######################
### ADD COMPARISONS ###
metadata_analysis$metadata_enhanced <- metadata_analysis$metadata

metadata_analysis$metadata_enhanced$anemia_v_rest <- as.character(metadata_analysis$metadata_enhanced$group)
metadata_analysis$metadata_enhanced$anemia_v_rest[metadata_analysis$metadata_enhanced$anemia_v_rest != 'anemia'] <- 'treatment'



metadata_analysis$metadata_enhanced$dex_muscle_v_rest <- as.character(metadata_analysis$metadata_enhanced$group)
metadata_analysis$metadata_enhanced$dex_muscle_v_rest[metadata_analysis$metadata_enhanced$dex_muscle_v_rest != 'fe_dextran_muscle
'] <- 'non_dextran_muscle'


# Fe dextran muscle
metadata_analysis$metadata_enhanced$ret_k_uL_cat <- ifelse(
  metadata_analysis$metadata_enhanced$ret_k_uL >= median(metadata_analysis$metadata_enhanced$ret_k_uL, na.rm=T), 
  'high', 
  'low')



metadata_analysis$metadata_enhanced$plasma_iron_ug_dl_cat <- ifelse(
  metadata_analysis$metadata_enhanced$plasma_iron_ug_dl >= median(metadata_analysis$metadata_enhanced$plasma_iron_ug_dl, na.rm=T), 
  'high', 
  'low')



temp_hgb_g_quant <- quantile(metadata_analysis$metadata_enhanced$hgb_g_dL)

metadata_analysis$metadata_enhanced$hgb_g_dL_cat <- ifelse(
  metadata_analysis$metadata_enhanced$hgb_g_dL <= temp_hgb_g_quant[names(temp_hgb_g_quant) == '25%'], 
  '1st_quantile', 
  
  ifelse(
    metadata_analysis$metadata_enhanced$hgb_g_dL <= temp_hgb_g_quant[names(temp_hgb_g_quant) == '50%'], 
    '2st_quantile', 
    
    ifelse(
      metadata_analysis$metadata_enhanced$hgb_g_dL <= temp_hgb_g_quant[names(temp_hgb_g_quant) == '75%'], 
      '3st_quantile', 
      
      ifelse(
        metadata_analysis$metadata_enhanced$hgb_g_dL <= temp_hgb_g_quant[names(temp_hgb_g_quant) == '100%'], 
        '4st_quantile', 
        NA
        )
      )
    )
  )



### ADD COMPARISONS ###
#######################






########################################
### FINAL TOUCHES FOR FINAL METADATA ###

colnames(metadata_analysis$metadata_enhanced) <- stringr::str_replace_all(
  string = colnames(metadata_analysis$metadata_enhanced), 
  pattern = '-', 
  replacement = '_')

colnames(metadata_analysis$metadata_enhanced) <- tolower(colnames(metadata_analysis$metadata_enhanced))

readr::write_tsv(x = metadata_analysis$metadata_enhanced, file = '../metadata_enhanced.tsv')

### FINAL TOUCHES FOR FINAL METADATA ###
########################################




















# metadata_analysis$corr$no_nas_pca_fil$data <- dplyr::select(.data = metadata_analysis$corr$no_nas$data, -plt_k_uL, -ret_k_uL)
# 
# metadata_analysis$corr$no_nas_pca_fil <- calculations_of_metadata_plots(metadata_analysis$corr$no_nas_pca_fil$data, dataset_name = 'no_nas_pca_fil')
# 
# metadata_analysis$corr$no_nas_pca_fil$habillage <- metadata_analysis$corr$no_nas$habillage