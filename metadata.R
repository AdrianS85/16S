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





# Pearson will be good - group count is to small to see normality, but I think it should be normal, cause why not
##################################
### REMOVE NAS FOR CORRELATION ###

metadata_analysis$metadata_for_corr$all <- metadata_analysis$metadata
metadata_analysis$metadata_for_corr$all$group <- NULL
rownames(metadata_analysis$metadata_for_corr$all) <- metadata_analysis$metadata_for_corr$all$SampleID
metadata_analysis$metadata_for_corr$all$SampleID <- NULL

metadata_analysis$metadata_for_corr$no_nas <- metadata_analysis$metadata_for_corr$all[-1,-1]

### REMOVE NAS FOR CORRELATION ###
##################################

metadata_analysis$corr$data <- Hmisc::rcorr(as.matrix(metadata_analysis$metadata_for_corr$no_nas))

##############################
### VISUALIZE CORRELATIONS ###

png(filename = paste0(metadata_analysis$dir, '/correlation_metadata.png'), width = 1600, height = 900)
ggcorrplot::ggcorrplot(metadata_analysis$corr$data$r, 
                       hc.order = TRUE, 
                       type = "lower", 
                       lab = TRUE, 
                       p.mat = ggcorrplot::cor_pmat(metadata_analysis$corr$data$r))
dev.off()

png(filename = paste0(metadata_analysis$dir, '/correlation_metadata_2.png'), width = 1600, height = 900)
psych::pairs.panels(metadata_analysis$metadata_for_corr$no_nas)
dev.off()

### VISUALIZE CORRELATIONS ###
##############################






###########################
### VARIABLES-BASED PCA ###

metadata_analysis$pca_scaled$data_variables <- prcomp(t(metadata_analysis$metadata_for_corr$no_nas), scale = TRUE)

png(filename = paste0(metadata_analysis$dir, '/pca_metadata_variables.png'), width = 1600, height = 900)
factoextra::fviz_pca_ind(metadata_analysis$pca_scaled$data_variables, repel = T)
dev.off()



metadata_analysis$metadata_for_corr$no_nas_pca_filt <- dplyr::select(.data = metadata_analysis$metadata_for_corr$no_nas, -plt_k_uL, -ret_k_uL)

metadata_analysis$pca_scaled$data_filt_variables <- prcomp(t(metadata_analysis$metadata_for_corr$no_nas_pca_filt), scale = TRUE)

png(filename = paste0(metadata_analysis$dir, '/pca_metadata_filt_variables.png'), width = 1600, height = 900)
factoextra::fviz_pca_ind(metadata_analysis$pca_scaled$data_filt_variables, repel = T)
dev.off()

### VARIABLES-BASED PCA ###
###########################



#########################
### SAMPLES-BASED PCA ###

metadata_analysis$pca_scaled$data_samples <- prcomp(metadata_analysis$metadata_for_corr$no_nas, scale = TRUE)

png(filename = paste0(metadata_analysis$dir, '/pca_metadata_samples.png'), width = 1600, height = 900)
factoextra::fviz_pca_biplot(metadata_analysis$pca_scaled$data_samples, repel = T, habillage = metadata_analysis$metadata[-1,]$group)
dev.off()

### SAMPLES-BASED PCA ###
#########################








# different variables: 1) plt_k_uL, 2) ret_k_uL. 3) plasma_iron_ug_dl, 4) non-heme_fe_liver_mg_fe_kg, mcv_fL, 5) rest



