### needs doParallel, foreach or doSNOW

library(ggplot2)
colnames(metadata_analysis$metadata)
divnet_analysis <- list(
  'import_dir' = list(
    'feat' = '/home/adrian/Desktop/nx_temp/nx_data_q/', 
    'tree' = '/home/adrian/Desktop/qiime/nx_phylogeny/',
    'tax' = '/home/adrian/Desktop/qiime/nx_taxonomy/',
    'metadata' = '/home/adrian/Desktop/qiime/qiime_metadata_nf.tsv'),
  'cols_to_analyze' = c('group', 'plasma_iron_ug_dl'),
  'cores' = 3,
  'tuning' = 'fast',
  'divnet_dir' = 'divnet',
  'vis_dir' = 'divnet/vis',
  'alpha_metrics' = c('shannon', 'simpson')
)

dir.create(divnet_analysis$divnet_dir)
dir.create(divnet_analysis$vis_dir)

### !!! add sklearn!
divnet_analysis$datasets <- list(
  'dada2' = qiime2R::qza_to_phyloseq(
    features = paste0(divnet_analysis$import_dir$feat, 'dada2.FeatureTableFrequency.qza'), 
    tree = paste0(divnet_analysis$import_dir$tree, 'sepp_dada2.PhylogenyRooted.qza'),
    taxonomy = paste0(divnet_analysis$import_dir$tax, 'for_divnet_vsearch_filtered_sepp_dada2.FeatureDataTaxonomy.qza'),
    metadata = divnet_analysis$import_dir$metadata),
  'deblur' = qiime2R::qza_to_phyloseq(
    features = paste0(divnet_analysis$import_dir$feat, 'deblur.FeatureTableFrequency.qza'), 
    tree = paste0(divnet_analysis$import_dir$tree, 'sepp_deblur.PhylogenyRooted.qza'),
    taxonomy = paste0(divnet_analysis$import_dir$tax, 'for_divnet_vsearch_filtered_sepp_deblur.FeatureDataTaxonomy.qza'),
    metadata = divnet_analysis$import_dir$metadata),
  'deblur_s' = qiime2R::qza_to_phyloseq(
    features = paste0(divnet_analysis$import_dir$feat, 'deblur_singletons.FeatureTableFrequency.qza'), 
    tree = paste0(divnet_analysis$import_dir$tree, 'sepp_deblur_singletons.PhylogenyRooted.qza'),
    taxonomy = paste0(divnet_analysis$import_dir$tax, 'for_divnet_vsearch_filtered_sepp_deblur_singletons.FeatureDataTaxonomy.qza'),
    metadata = divnet_analysis$import_dir$metadata)
)






divnet_analysis$divnet <- purrr::map(
  .x = divnet_analysis$datasets,
  .f = function(dataset){

    boop <- purrr::map(
      .x = divnet_analysis$cols_to_analyze,
      .f = function(col_name){

        DivNet::divnet(X = col_name, ncores = divnet_analysis$cores, tuning = divnet_analysis$tuning, W = dataset)
      })

    names(boop) <- divnet_analysis$cols_to_analyze

    return(boop)
  })
names(divnet_analysis$divnet) <- names(divnet_analysis$datasets)






purrr::walk(
  .x = divnet_analysis$alpha_metrics,
  .f = function(metric){
    
    purrr::pwalk(.l = list(
      divnet_analysis$divnet,
      divnet_analysis$datasets,
      names(divnet_analysis$datasets)
    ),
    .f = function(divnet_dataset, original_dataset, dataset_name){
      
      purrr::map2(
        .x = divnet_dataset,
        .y = names(divnet_dataset),
        .f = function(col, col_name){
          
          png_name <- paste0(divnet_analysis$vis_dir, '/', col_name, '_', dataset_name, '_', metric, '.png')
          
          x <- col[[metric]] %>%
            plot(original_dataset, color = col_name) +
            xlab(col_name) +
            ylab(metric)
          
          ggsave(filename = png_name, plot = x)
        })
    })
  })






divnet_analysis$beta <- purrr::pmap(.l = list(divnet_analysis$divnet, divnet_analysis$dataset, names(divnet_analysis$divnet)),
  .f = function(divnet_dataset, original_dataset, dataset_name){
    
    purrr::walk2(
      .x = divnet_dataset,
      .y = names(divnet_dataset),
      .f = function(cols, cols_names){
        
        png_name <- paste0(divnet_analysis$vis_dir, '/', cols_names, '_', dataset_name, '_bray_curtis.png')
        
        x <- ggplot(
          data = DivNet::simplifyBeta(dv = cols, physeq = original_dataset, measure = "bray-curtis", x = cols_names), 
          mapping = aes(x = interaction(Covar1, Covar2), y = beta_est, col = interaction(Covar1, Covar2))) +
          geom_point() +
          geom_linerange(aes(ymin = lower, ymax = upper)) + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          xlab("") + ylab("Estimates of Bray-Curtis distance")
        
        ggsave(filename = png_name, plot = x)
      })
    
    purrr::map2(
      .x = divnet_dataset,
      .y = names(divnet_dataset),
      .f = function(cols, cols_names){
        
        DivNet::simplifyBeta(dv = cols, physeq = original_dataset, measure = "bray-curtis", x = cols_names)
      })
})



  


divnet_analysis$hypothesis <- purrr::map(
  .x = divnet_analysis$alpha_metrics,
  .f = function(metric){

  purrr::pmap(
    .l = list(divnet_analysis$divnet, divnet_analysis$dataset, names(divnet_analysis$divnet)),
    .f = function(divnet_dataset, original_dataset, dataset_name){
      
      purrr::map2(
        .x = divnet_dataset,
        .y = names(divnet_dataset),
        .f = function(cols, cols_names){
          
          matric_var_name <- paste0(metric, '-variance')
          
          estimates <- cols[[metric]] %>% summary %$% estimate
          ses <- sqrt(cols[[matric_var_name]])
          X <- breakaway::make_design_matrix(original_dataset, cols_names)

          return(breakaway::betta(estimates, ses, X)$table)
      })
  })
})
names(divnet_analysis$hypothesis) <- divnet_analysis$alpha_metrics
















#####################################
### HYPOTHESIS ASSUMPTION TESTING ###

### This looks like it only pertians to continous data
# check model: samples (without group markers) to numeric and plot numerics
divnet_analysis$divnet_model_test <- purrr::map(
  .x = divnet_analysis$datasets, 
  .f = function(dataset){
    
    dataset$sam_data$group <- as.numeric(dataset$sam_data$group)
    
    return_ <- DivNet::divnet(X = dataset_continuous, ncores = 3, tuning = "careful", W = dataset, base = 1) # what does the base do really?
    
    return(return_)
  })

phyloseq::sample_data(object = divnet_analysis$datasets$dada2) %$% group %>% as.numeric()

divnet_analysis$datasets$dada2 %>% phyloseq::sample_data %$% group %>% as.character %>% as.numeric


shannon_table <- divnet_temperature %$%
  shannon %>%
  summary %>%
  cbind(temp_continuous)
shannon_table %>%
  ggplot(aes(x = temp_continuous, y = estimate)) + geom_point() +
  geom_linerange(aes(x = temp_continuous, ymin = lower, ymax = upper))

### HYPOTHESIS ASSUMPTION TESTING ###
#####################################

# 'formulas' = 
#   list('~ group' = list(
#     'f_d_m' = c('group', 'fe_dextran_muscle', 'anemia'),
#     's_f_l_o' = c('group', 'siderai_fe_liposomal_oral', 'anemia'))))

