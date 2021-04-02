### needs doParallel, foreach or doSNOW

# sampleid	group	mother	non_heme_fe_faces_mg_fe_kg	non_heme_fe_liver_mg_fe_kg	ret_k_ul	plt_k_ul	ret_he_pg	plasma_iron_ug_dl	anemia_v_treatment	dex_muscle_v_treatment	rbc_morphology_z_score	rbc_quant_z_score

save(divnet_analysis, file = 'divnet_analysis.save')
load('divnet_analysis.save')



source('functions.R')
library(ggplot2)
library(magrittr)
data(Lee) # DivNet


divnet_analysis <- list(
  'import_dir' = list(
    'feat' = '/home/adrian/Desktop/qiime/full_data_server/nx_data_q/', 
    'tree' = '/home/adrian/Desktop/qiime/full_data_server/nx_phylogeny/',
    'tax' = '/home/adrian/Desktop/qiime/full_data_server/nx_taxonomy/',
    'metadata' = '/home/adrian/Desktop/qiime/metadata_enhanced_phylo.tsv'),
  'formulas' = list(
    'group' = as.formula("~ group"),
    'anemia_v_treatment' = as.formula("~ anemia_v_treatment"),
    'dex_muscle_v_treatment' = as.formula("~ dex_muscle_v_treatment")
    ),
  'sample_to_remove_for_formulas' = 
    list(
      'group' = NA,
      'anemia_v_treatment' = NA,
      'dex_muscle_v_treatment' = c('A1', 'A2', 'A3', 'A4', 'A5')),
  'cores' = 4,
  'tuning' = 'fast',
  'divnet_dir' = 'divnet',
  'vis_dir' = 'divnet/vis',
  'vis_break_dir' = 'divnet/break',
  'vis_alpha_dir' = 'divnet/alpha',
  'vis_beta_dir' = 'divnet/simplifybeta',
  'alpha_metrics' = c('shannon', 'simpson')
)

dir.create(divnet_analysis$divnet_dir)
dir.create(divnet_analysis$vis_dir)
dir.create(divnet_analysis$vis_break_dir)
dir.create(divnet_analysis$vis_alpha_dir)
dir.create(divnet_analysis$vis_beta_dir)





divnet_analysis$datasets <- list(
  'dada2' = qiime2R::qza_to_phyloseq(
    features = paste0(divnet_analysis$import_dir$feat, 'dada2.FeatureTableFrequency.qza'), 
    tree = paste0(divnet_analysis$import_dir$tree, 'sepp_dada2.PhylogenyRooted.qza'),
    taxonomy = paste0(divnet_analysis$import_dir$tax, 'for_divnet_vsearch_dada2.FeatureDataTaxonomy.qza'),
    metadata = divnet_analysis$import_dir$metadata),
  'deblur' = qiime2R::qza_to_phyloseq(
    features = paste0(divnet_analysis$import_dir$feat, 'deblur.FeatureTableFrequency.qza'), 
    tree = paste0(divnet_analysis$import_dir$tree, 'sepp_deblur.PhylogenyRooted.qza'),
    taxonomy = paste0(divnet_analysis$import_dir$tax, 'for_divnet_vsearch_deblur.FeatureDataTaxonomy.qza'),
    metadata = divnet_analysis$import_dir$metadata)
)







temp <- purrr::map2(
  .x = divnet_analysis$datasets,
  .y = names(divnet_analysis$datasets),
  .f = function(dataset, dataset_name){
  
    temp_batch <- MMUPHin::adjust_batch(
      feature_abd = dataset@otu_table, 
      batch = 'mother', 
      data = data.frame(dataset@sam_data),
      control = list('zero_inflation' = T, 'pseudo_count' = NULL, 'diagnostic_plot' = paste0(divnet_analysis$vis_dir, '/', dataset_name, '_batch_adj.pdf')))[['feature_abd_adj']]
  
    return(
    phyloseq::phyloseq(
      phyloseq::otu_table(temp_batch, taxa_are_rows = T), 
      phyloseq::phy_tree(dataset@phy_tree), 
      phyloseq::tax_table(dataset@tax_table), 
      phyloseq::sample_data(dataset@sam_data))
    )
})
names(temp) <- paste0('batch_corr_', names(divnet_analysis$datasets))

divnet_analysis$datasets <- rlist::list.flatten(list(temp, divnet_analysis$datasets))


divnet_analysis$datasets_taxa <- purrr::map(
  .x = divnet_analysis$datasets,
  .f = function(dataset) {
    
    purrr::map(
      .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'),
      .f = function(tax_level) {
        
        phyloseq::tax_glom(physeq = dataset, taxrank = tax_level) 
      })
  })
divnet_analysis$datasets_asv_taxa <- rlist::list.flatten(list(divnet_analysis$datasets, divnet_analysis$datasets_taxa))









divnet_analysis$datasets_no_single <- list(
  'dada2' = qiime2R::qza_to_phyloseq(
    features = paste0(divnet_analysis$import_dir$feat, 'filtered_sepp_dada2.FeatureTableFrequency.qza'), 
    tree = paste0(divnet_analysis$import_dir$tree, 'sepp_dada2.PhylogenyRooted.qza'),
    taxonomy = paste0(divnet_analysis$import_dir$tax, 'vsearch_filtered_sepp_dada2.FeatureDataTaxonomy.qza'),
    metadata = divnet_analysis$import_dir$metadata),
  'deblur' = qiime2R::qza_to_phyloseq(
    features = paste0(divnet_analysis$import_dir$feat, 'filtered_sepp_deblur.FeatureTableFrequency.qza'), 
    tree = paste0(divnet_analysis$import_dir$tree, 'sepp_deblur.PhylogenyRooted.qza'),
    taxonomy = paste0(divnet_analysis$import_dir$tax, 'vsearch_filtered_sepp_deblur.FeatureDataTaxonomy.qza'),
    metadata = divnet_analysis$import_dir$metadata)
)

temp_2 <- purrr::map2(
  .x = divnet_analysis$datasets_no_single,
  .y = names(divnet_analysis$datasets_no_single),
  .f = function(dataset, dataset_name){
    
    temp_batch <- MMUPHin::adjust_batch(
      feature_abd = dataset@otu_table, 
      batch = 'mother', 
      data = data.frame(dataset@sam_data),
      control = list('zero_inflation' = T, 'pseudo_count' = NULL, 'diagnostic_plot' = paste0(divnet_analysis$vis_dir, '/', dataset_name, '_batch_adj.pdf')))[['feature_abd_adj']]
    
    return(
      phyloseq::phyloseq(
        phyloseq::otu_table(temp_batch, taxa_are_rows = T), 
        phyloseq::phy_tree(dataset@phy_tree), 
        phyloseq::tax_table(dataset@tax_table), 
        phyloseq::sample_data(dataset@sam_data))
    )
  })
names(temp_2) <- paste0('batch_corr_', names(divnet_analysis$datasets_no_single))
divnet_analysis$datasets_no_single <- rlist::list.flatten(list(temp_2, divnet_analysis$datasets_no_single))


divnet_analysis$datasets_no_single_taxa <- purrr::map(
  .x = divnet_analysis$datasets_no_single,
  .f = function(dataset) {
    
    purrr::map(
      .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'),
      .f = function(tax_level) {
        
        phyloseq::tax_glom(physeq = dataset, taxrank = tax_level) 
      })
  })
divnet_analysis$datasets_no_single_asv_taxa <- rlist::list.flatten(list(divnet_analysis$datasets_no_single, divnet_analysis$datasets_no_single_taxa))











### It returns the significance of the covariates in explaining diversity and a hypothesis test for heterogeneity.
### check also filtered data for large disparities with data with singletons. The differences should not be huge breakaway with singletons vs breakaway nof without singletons. breakaway estmate richness, not evennness
### !!! if something interesting comes up, try also breakaway::objective_bayes_*
### Chyba można plotować wynik breakaway i czarne diamenty musza iść w tym plocie tam gdzie white circles, look here for how it should look https://adw96.github.io/breakaway/articles/intro-diversity-estimation.html
### try summary on breakaway objects
### I dont know whats happening with batch corrected datasets with breakaway. They have no error bars, and I dont get it
divnet_analysis$breakaway <- purrr::pmap(
  list(
    divnet_analysis$datasets_asv_taxa,
    names(divnet_analysis$datasets_asv_taxa),
    divnet_analysis$datasets_no_single_asv_taxa,
    names(divnet_analysis$datasets_no_single_asv_taxa)),
  .f = function(dataset_, dataset_names, dataset_no, dataset_no_names){
    
    breakaway_ <- breakaway::breakaway(dataset_)
    breakaway_plot <- plot(breakaway_, physeq = dataset_, color = 'group')
    
    ggsave(filename = paste0(divnet_analysis$vis_break_dir, '/', dataset_names, '_breakaway_plot.png'), plot = breakaway_plot, device = 'png', width = 9.6, height = 5.4)
    
    meta <- as.data.frame( phyloseq::sample_data(dataset_) )
    meta$sample_names <- rownames(meta)
    meta <- tibble::as_tibble(meta)
    
    combined_richness <- merge(meta, summary(breakaway_), by = "sample_names")
    
    breakaway_bettas <- purrr::map2(
      .x = divnet_analysis$formulas,
      .y = names(divnet_analysis$formulas),
      .f = function(formula, column){
        
        combined_richness <- subset(x = combined_richness, subset = !is.na(combined_richness[[column]]))
        
        
        tryCatch(
          {
            betta <- breakaway::betta(chats = combined_richness$estimate,
                                      ses = combined_richness$error,
                                      X = model.matrix(formula, data = combined_richness))$table
            
            return(betta)
          },
          error=function(cond) {
            message(paste0(cond, '\t for \t', column))

            return(NA)
          })
      })
    
    breakaway_nof1_ <- breakaway::breakaway_nof1(dataset_no)
    breakaway_nof1_plot <- plot(breakaway_nof1_, physeq = dataset_no, color = 'group')
    
    ggsave(filename = paste0(divnet_analysis$vis_break_dir, '/', dataset_no_names, '_nof1_plot.png'), plot = breakaway_nof1_plot, device = 'png', width = 9.6, height = 5.4)
    
    return(list('breakaway' = breakaway_, 'breakaway_bettas' = breakaway_bettas, 'nof1' = breakaway_nof1_))
  })

### divnet_analysis[["breakaway"]][["batch_corr_dada2.Species"]][["breakaway_bettas"]][["group"]] - groupnanofe_dextran_50nm_oral    0.041










divnet_analysis$datasets_taxa_single_and_no_single <- rlist::list.flatten(list('singl' = divnet_analysis$datasets_taxa, 'no_singl' = divnet_analysis$datasets_no_single_taxa))

divnet_analysis$divnet_taxa_single_and_no_single <- purrr::map(
  .x = divnet_analysis$datasets_taxa_single_and_no_single,
  .f = function(dataset){
    
    purrr::map2(
      .x = divnet_analysis$formulas,
      .y = names(divnet_analysis$formulas),
      .f = function(formula_, cov){
        
        if (cov == 'dex_muscle_v_treatment') {
          dataset <- phyloseq::subset_samples(physeq = dataset, !is.na(dex_muscle_v_treatment))
        }
        
        DivNet::divnet(W = dataset, X = cov, ncores = divnet_analysis$cores, tuning = divnet_analysis$tuning)
      })
  })









divnet_analysis$alpha_from_server <- purrr::map(
  .x = divnet_analysis$alpha_metrics,
  .f = function(metric){
    
    purrr::pmap(.l = list(
      divnet_analysis$divnet_from_server,
      divnet_analysis$datasets_asv_taxa,
      names(divnet_analysis$datasets_asv_taxa)
    ),
    .f = function(divnet_dataset, original_dataset, dataset_name){
      
      purrr::map2(
        .x = divnet_dataset,
        .y = names(divnet_dataset),
        .f = function(col, col_name){
          
          if (col_name == 'dex_muscle_v_treatment') {
            original_dataset <- phyloseq::subset_samples(physeq = original_dataset, !is.na(dex_muscle_v_treatment))
          }
          
          plot_name <- paste0(divnet_analysis_2$vis_alpha_dir, '/', col_name, '_', dataset_name, '_', metric, '.png')
          temp <- summary(col[[metric]])
          
          temp$mother <- recode_values_based_on_key(
            to_recode_chrvec = temp$sample_names, 
            replace_this_chrvec = rownames(original_dataset@sam_data), 
            with_this_chrvec = original_dataset@sam_data$mother)
          
          temp <- temp[order(temp$mother),] 
          temp$sample_names_factor <- factor(x = temp$sample_names, levels = temp$sample_names)
          
          main_plot <- ggplot(data = temp, mapping = aes(x = sample_names_factor, y = estimate, color = mother)) +
            geom_point() +
            geom_errorbar(aes(ymin = lower, ymax = upper), width=.1)
          
          ggsave(
            filename = plot_name, 
            plot = main_plot,
            width = 8.5, 
            height = 5.4)
          
          
          
          ### Main stat
          estimates <- col[[metric]] %>% summary %$% estimate
          ses <- sqrt(col[[ paste0(metric, '-variance') ]])
          X <- breakaway::make_design_matrix(original_dataset, col_name)
          betta <- breakaway::betta(estimates, ses, X)$table
          
          betta_df <- as.data.frame(betta[-1,])
          
          if (length(betta_df$`p-values`) > 1) {
            betta_df$`adj_p-values` <- p.adjust(p = betta_df$`p-values`, method = 'BH')
          }
          
          return(list('betta' = betta, 'betta_df' = betta_df))
        })
    })
  })
names(divnet_analysis$alpha_from_server) <- divnet_analysis$alpha_metrics


### Beta is pointless, cause no statistical method of analysis is given, and im not doing adonis again, when alpha metrics are of suspect anyway
divnet_analysis$beta <- purrr::pmap(
  .l = list(
    divnet_analysis$divnet_from_server,
    divnet_analysis$datasets_asv_taxa,
    names(divnet_analysis$divnet_from_server)),
  .f = function(divnet_dataset, original_dataset, dataset_name){
    
    purrr::map2(
      .x = divnet_dataset,
      .y = names(divnet_dataset),
      .f = function(cols, cols_names){
        
        png_name <- paste0(divnet_analysis$vis_beta_dir, '/', cols_names, '_', dataset_name, '_bray_curtis.png')
        
        
        if (cols_names == 'dex_muscle_v_treatment') {
          original_dataset = phyloseq::subset_samples(original_dataset, !is.na(dex_muscle_v_treatment))
        }
        if (cols_names == 'group') {
          ggplot2::ggsave(filename = png_name,
                          plot = factoextra::fviz_pca_ind(
                            prcomp(cols$`bray-curtis`),
                            repel = T,
                            habillage = data.frame(original_dataset@sam_data)$group),
                          width = 9.6,
                          height = 5.4)
        }
        
        simplifyBeta <- DivNet::simplifyBeta(dv = cols, physeq = original_dataset, measure = "bray-curtis", x = cols_names)
        
        
        return(simplifyBeta)
      })
  })




























### Before taking that result too seriously, let’s check for model misspecification. The DivNet model performs well when there is not significant curvature in the log-ratios. To investigate this, let’s pull out the observed log-ratios (observed_yy) as well as the fitted log-ratios (fitted_yy). So do this only for potentially interesting resluts

### Here we draw alpha divnet pics
### Assessing model fit - You do not need to do this if you are fitting DivNet with categorical covariates, only if you are fitting DivNet with continuous covariates.
divnet_analysis$alpha <- purrr::map(
  .x = divnet_analysis$alpha_metrics,
  .f = function(metric){
    
    purrr::pmap(.l = list(
      divnet_analysis$divnet,
      divnet_analysis$datasets_taxa_no_sp,
      names(divnet_analysis$datasets_taxa_no_sp)
    ),
    .f = function(divnet_dataset, original_dataset, dataset_name){
      
      purrr::map2(
        .x = divnet_dataset,
        .y = names(divnet_dataset),
        .f = function(col, col_name){
          
          if (col_name == 'dex_muscle_v_treatment') {
            original_dataset <- phyloseq::subset_samples(physeq = original_dataset, !is.na(dex_muscle_v_treatment))
          }
          
          plot_name <- paste0(divnet_analysis$vis_alpha_dir, '/', col_name, '_', dataset_name, '_', metric, '.png')
          
          main_plot <- col[[metric]] %>%
            plot(original_dataset, color = col_name) +
            xlab(col_name) +
            ylab(metric)
          
          ggsave(
            filename = plot_name, 
            plot = main_plot)
          
          
          
          ### Main stat
          estimates <- col[[metric]] %>% summary %$% estimate
          ses <- sqrt(col[[ paste0(metric, '-variance') ]])
          X <- breakaway::make_design_matrix(original_dataset, col_name)
          betta <- breakaway::betta(estimates, ses, X)$table
          
          betta_df <- as.data.frame(betta[-1,])
          
          if (length(betta_df$`p-values`) > 1) {
            betta_df$`adj_p-values` <- p.adjust(p = betta_df$`p-values`, method = 'BH')
          }
          
          return(list('betta' = betta, 'betta_df' = betta_df))
        })
    })
  })
names(divnet_analysis$alpha) <- divnet_analysis$alpha_metrics
### Plot values depend on model chosen for given gene. Therefore it is highly unreasonable to compare them using a plot
### Values for shannon and simpson are in correct range (https://www.biorxiv.org/content/10.1101/305045v1.full.pdf)
### Values for shannon are similar to what I see in qiime-based shannon estimate, but p-values are insanely low
### Shannon and Simpson agree as to which groups are where compared to other, as they are negatively correlated
### Co zrobić z nierównymi grupami anemia/dex?
### https://arxiv.org/pdf/1506.05710.pdf

### follow this https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-020-01723-9

### https://peerj.com/articles/9350/

save(divnet_analysis, file = 'divnet_analysis_with_server.save')
























divnet_analysis$divnet_no_single$batch_corr_dada2 <- purrr::map2(
  .x = divnet_analysis$formulas,
  .y = names(divnet_analysis$formulas),
  .f = function(formula_, cov){
    
    dataset <- divnet_analysis$datasets_no_single$batch_corr_dada2
    
    if (cov == 'dex_muscle_v_treatment') {
      dataset <- phyloseq::subset_samples(physeq = dataset, !is.na(dex_muscle_v_treatment))
    }
    
    DivNet::divnet(W = dataset, X = cov, ncores = divnet_analysis$cores, tuning = divnet_analysis$tuning)
  })



divnet_analysis$alpha_no_single$batch_corr_dada2 <- purrr::map(
  .x = divnet_analysis$alpha_metrics,
  .f = function(metric){
    
    original_dataset <- divnet_analysis$datasets_no_single$batch_corr_dada2
    divnet_dataset <- divnet_analysis$divnet_no_single$batch_corr_dada2
    dataset_name <- names(divnet_analysis$divnet_no_single$batch_corr_dada2)
    
    
    
    purrr::map2(
      .x = divnet_dataset,
      .y = names(divnet_dataset),
      .f = function(col, col_name){
        
        if (col_name == 'dex_muscle_v_treatment') {
          original_dataset <- phyloseq::subset_samples(physeq = original_dataset, !is.na(dex_muscle_v_treatment))
        }
        
        plot_name <- paste0(divnet_analysis$vis_alpha_dir, '/no_single_', col_name, '_', metric, '.png')
        
        main_plot <- col[[metric]] %>%
          plot(original_dataset, color = col_name) +
          xlab(col_name) +
          ylab(metric)
        
        ggsave(
          filename = plot_name, 
          plot = main_plot, width = 9.6, height = 5.4)
        
        
        
        ### Main stat
        estimates <- col[[metric]] %>% summary %$% estimate
        ses <- sqrt(col[[ paste0(metric, '-variance') ]])
        X <- breakaway::make_design_matrix(original_dataset, col_name)
        betta <- breakaway::betta(estimates, ses, X)$table
        
        betta_df <- as.data.frame(betta[-1,])
        
        if (length(betta_df$`p-values`) > 1) {
          betta_df$`adj_p-values` <- p.adjust(p = betta_df$`p-values`, method = 'BH')
        }
        
        return(list('betta' = betta, 'betta_df' = betta_df))
      })
  })
names(divnet_analysis$alpha_no_single$batch_corr_dada2) <- divnet_analysis$alpha_metrics
# divnet_analysis[["alpha_no_single"]][["batch_corr_dada2"]][["shannon"]][["group"]][["betta"]]


### Beta is pointless, cause no statistical method of analysis is given, and im not doing adonis again, when alpha metrics are of suspect anyway
divnet_analysis$beta_no_single$batch_corr_dada2 <- purrr::map2(
  .x = divnet_analysis$divnet_no_single$batch_corr_dada2,
  .y = names(divnet_analysis$divnet_no_single$batch_corr_dada2),
  .f = function(cols, cols_names){
    original_dataset <- divnet_analysis$datasets_no_single$batch_corr_dada2
    png_name <- paste0(divnet_analysis$vis_beta_dir, '/beta_no_single_', cols_names, '_bray_curtis.png')
    
    
    if (cols_names == 'dex_muscle_v_treatment') {
      original_dataset = phyloseq::subset_samples(original_dataset, !is.na(dex_muscle_v_treatment))
    }
    if (cols_names == 'group') {
      ggplot2::ggsave(filename = png_name,
                      plot = factoextra::fviz_pca_ind(
                        prcomp(cols$`bray-curtis`),
                        repel = T,
                        habillage = data.frame(original_dataset@sam_data)$group),
                      width = 9.6,
                      height = 5.4)
    }
    
    simplifyBeta <- DivNet::simplifyBeta(dv = cols, physeq = original_dataset, measure = "bray-curtis", x = cols_names)
    
    
    return(simplifyBeta)
  })

divnet_analysis$std_alpha$rare <- phyloseq::rarefy_even_depth(physeq = divnet_analysis$datasets_asv_taxa$batch_corr_dada2dada2)

divnet_analysis$std_alpha$est <- phyloseq::estimate_richness(
  physeq = divnet_analysis$std_alpha$rare, 
  measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))

divnet_analysis$std_alpha$est$group <- recode_values_based_on_key(
  to_recode_chrvec = rownames(divnet_analysis$std_alpha$est), 
  replace_this_chrvec = rownames(divnet_analysis$std_alpha$rare@sam_data), 
  with_this_chrvec = divnet_analysis$std_alpha$rare@sam_data$group)

summary(aov(Shannon ~ group, data = divnet_analysis$std_alpha$est))
DescTools::PostHocTest(aov(Shannon ~ group, data = divnet_analysis$std_alpha$est), method = c("hsd"), conf.level = 0.95)

kruskal.test(Observed ~ group, data = divnet_analysis$std_alpha$est)

phyloseq::plot_richness(physeq = divnet_analysis_2$std_alpha$rare, x = 'group', color = 'group', measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) + ggplot2::geom_boxplot()
###












































######################################
### MOTHER DADA2 CENTERED ANALYSIS ###

save(divnet_analysis_2, file = 'divnet_analysis_2.save')
load('divnet_analysis_2.save')



source('functions.R')
library(ggplot2)
library(magrittr)



divnet_analysis_2 <- list(
  'import_dir' = list(
    'feat' = '/home/adrian/Desktop/qiime/full_data_server/nx_data_q/', 
    'tree' = '/home/adrian/Desktop/qiime/full_data_server/nx_phylogeny/',
    'tax' = '/home/adrian/Desktop/qiime/full_data_server/nx_taxonomy/',
    'metadata' = '/home/adrian/Desktop/qiime/metadata_enhanced_phylo.tsv'),
  'formulas' = list(
    'mother' = as.formula("~ mother")
  ),
  'sample_to_remove_for_formulas' = 
    list(
      'mother' = NA),
  'cores' = 4,
  'tuning' = 'fast',
  'divnet_dir' = 'divnet_2',
  'vis_dir' = 'divnet_2/vis',
  'vis_break_dir' = 'divnet_2/break',
  'vis_alpha_dir' = 'divnet_2/alpha',
  'vis_beta_dir' = 'divnet_2/simplifybeta',
  'alpha_metrics' = c('shannon', 'simpson')
)



dir.create(divnet_analysis_2$divnet_dir)
dir.create(divnet_analysis_2$vis_dir)
dir.create(divnet_analysis_2$vis_break_dir)
dir.create(divnet_analysis_2$vis_alpha_dir)
dir.create(divnet_analysis_2$vis_beta_dir)



divnet_analysis_2$datasets <- list(
  'dada2' = qiime2R::qza_to_phyloseq(
    features = paste0(divnet_analysis_2$import_dir$feat, 'dada2.FeatureTableFrequency.qza'), 
    tree = paste0(divnet_analysis_2$import_dir$tree, 'sepp_dada2.PhylogenyRooted.qza'),
    taxonomy = paste0(divnet_analysis_2$import_dir$tax, 'for_divnet_vsearch_dada2.FeatureDataTaxonomy.qza'),
    metadata = divnet_analysis_2$import_dir$metadata)
)

divnet_analysis_2$datasets$dada2 <- phyloseq:::subset_samples(divnet_analysis_2$datasets$dada2, mother != "7635_16") ### Remove mother with only one piglet

temp_batch <- MMUPHin::adjust_batch(
  feature_abd = divnet_analysis_2$datasets$dada2@otu_table, 
  batch = 'group', 
  data = data.frame(divnet_analysis_2$datasets$dada2@sam_data),
  control = list('zero_inflation' = T, 'pseudo_count' = NULL, 'diagnostic_plot' = paste0(divnet_analysis_2$vis_dir, '/dada2_mother_batch_adj.pdf')))[['feature_abd_adj']]

divnet_analysis_2$datasets$batch_corr_mother_dada2 <- phyloseq::phyloseq(
  phyloseq::otu_table(temp_batch, taxa_are_rows = T), 
  phyloseq::phy_tree(divnet_analysis_2$datasets$dada2@phy_tree), 
  phyloseq::tax_table(divnet_analysis_2$datasets$dada2@tax_table), 
  phyloseq::sample_data(divnet_analysis_2$datasets$dada2@sam_data))

divnet_analysis_2$datasets_taxa <- purrr::map(
  .x = divnet_analysis_2$datasets,
  .f = function(dataset) {
    
    purrr::map(
      .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'),
      .f = function(tax_level) {
        
        phyloseq::tax_glom(physeq = dataset, taxrank = tax_level) 
      })
  })
divnet_analysis_2$datasets_asv_taxa <- rlist::list.flatten(list(divnet_analysis_2$datasets, divnet_analysis_2$datasets_taxa))




# divnet_analysis_2$which_datasets_to_analyze <- stringr::str_detect(string = names(divnet_analysis_2$datasets_asv_taxa), pattern = 'Class|Family')

# divnet_analysis_2$divnet <- purrr::map(
#   .x = divnet_analysis_2$datasets_asv_taxa[divnet_analysis_2$which_datasets_to_analyze],
#   .f = function(dataset){
#     
#     purrr::map2(
#       .x = divnet_analysis_2$formulas,
#       .y = names(divnet_analysis_2$formulas),
#       .f = function(formula_, cov){
# 
#         DivNet::divnet(W = dataset, X = cov, ncores = divnet_analysis_2$cores, tuning = divnet_analysis_2$tuning)
#       })
#   })




### Done on calculation server ###
divnet_analysis_2$divnet <- DivNet::divnet(
  W = divnet_analysis_2$datasets_asv_taxa$batch_corr_mother_dada2, 
  X = "mother", 
  ncores = divnet_analysis_2$cores, 
  tuning = divnet_analysis_2$tuning)

load('divnet_analysis_2_calculated.save')

divnet_analysis_2$divnet$batch_corr_mother_dada2 <- list()

divnet_analysis_2$divnet$batch_corr_mother_dada2$mother <- divnet_analysis_2_calculated$divnet
divnet_analysis_2$divnet[[1]] <- NULL

divnet_analysis_2$which_datasets_to_analyze <- names(divnet_analysis_2$datasets_asv_taxa) == 'batch_corr_mother_dada2'

rm(divnet_analysis_2_calculated)
### Done on calculation server ###





divnet_analysis_2$alpha <- purrr::map(
  .x = divnet_analysis_2$alpha_metrics,
  .f = function(metric){
    
    
    
    purrr::pmap(.l = list(
      divnet_analysis_2$divnet,
      divnet_analysis_2$datasets_asv_taxa[divnet_analysis_2$which_datasets_to_analyze],
      names(divnet_analysis_2$datasets_asv_taxa[divnet_analysis_2$which_datasets_to_analyze])
    ),
    .f = function(divnet_dataset, original_dataset, dataset_name){
      
      
      
      purrr::map2(
        .x = divnet_dataset,
        .y = names(divnet_dataset),
        .f = function(col, col_name){
          
          plot_name <- paste0(divnet_analysis_2$vis_alpha_dir, '/', col_name, '_', dataset_name, '_', metric, '.png')
          temp <- summary(col[[metric]])
          
          temp$mother <- recode_values_based_on_key(
            to_recode_chrvec = temp$sample_names, 
            replace_this_chrvec = rownames(original_dataset@sam_data), 
            with_this_chrvec = original_dataset@sam_data$mother)
          
          temp <- temp[order(temp$mother),] 
          temp$sample_names_factor <- factor(x = temp$sample_names, levels = temp$sample_names)
          
          main_plot <- ggplot(data = temp, mapping = aes(x = sample_names_factor, y = estimate, color = mother)) +
            geom_point() +
            geom_errorbar(aes(ymin = lower, ymax = upper), width=.1)
          
          ggsave(
            filename = plot_name, 
            plot = main_plot,
            width = 8.5, 
            height = 5.4)
          
          main_stat <- purrr::map(
            .x = levels(original_dataset@sam_data$mother), 
            .f = function(level)
            {
              original_dataset@sam_data$mother <- relevel(x = original_dataset@sam_data$mother, level)
              
              estimates <- col[[metric]] %>% summary %$% estimate
              ses <- sqrt(col[[ paste0(metric, '-variance') ]])
              X <- breakaway::make_design_matrix(original_dataset, col_name)
              betta <- breakaway::betta(estimates, ses, X)$table
              
              betta_df <- as.data.frame(betta[-1,])
              
              if (length(betta_df$`p-values`) > 1) {
                betta_df$`adj_p-values` <- p.adjust(p = betta_df$`p-values`, method = 'BH')
              }
              
              return(list('betta' = betta, 'betta_df' = betta_df))
              
            })
          names(main_stat) <- levels(original_dataset@sam_data$mother)
          
          return(main_stat)
        })
    })
  })
names(divnet_analysis_2$alpha) <- divnet_analysis_2$alpha_metrics
### !!! does divnet needs rarefaction?


### Beta is pointless, cause no statistical method of analysis is given, and im not doing adonis again, when alpha metrics are of suspect anyway
divnet_analysis_2$beta <- purrr::pmap(
  .l = list(
    divnet_analysis_2$divnet,
    divnet_analysis_2$datasets_asv_taxa[divnet_analysis_2$which_datasets_to_analyze],
    names(divnet_analysis_2$datasets_asv_taxa[divnet_analysis_2$which_datasets_to_analyze])),
  .f = function(divnet_dataset, original_dataset, dataset_name){
    
    
    
    purrr::map2(
      .x = divnet_dataset,
      .y = names(divnet_dataset),
      .f = function(cols, cols_names){
        
        png_name <- paste0(divnet_analysis_2$vis_beta_dir, '/', cols_names, '_', dataset_name, '_bray_curtis.png')
        
          ggplot2::ggsave(filename = png_name,
                          plot = factoextra::fviz_pca_var(
                            prcomp(cols$`bray-curtis`),
                            repel = T,
                            habillage = data.frame(original_dataset@sam_data)$mother),
                          width = 9.6,
                          height = 5.4)

        simplifyBeta <- DivNet::simplifyBeta(dv = cols, physeq = original_dataset, measure = "bray-curtis", x = cols_names)
        
        
        return(simplifyBeta)
      })
  })
### !!! tutaj mamy dystans b-c, więc na nich się mds robi, nie pca. tylko pytanie, czy b-c to euklideański dystans? metric MDS is fine for b-c: https://www.researchgate.net/figure/Metric-multidimensional-scaling-MDS-analysis-using-Bray-Curtis-dissimilarity-A-and_fig4_325303905


divnet_analysis_2$std_alpha$rare <- phyloseq::rarefy_even_depth(physeq = divnet_analysis_2$datasets_asv_taxa$batch_corr_mother_dada2)

divnet_analysis_2$std_alpha$est <- phyloseq::estimate_richness(
  physeq = divnet_analysis_2$std_alpha$rare, 
  measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))

divnet_analysis_2$std_alpha$est$mother <- recode_values_based_on_key(
  to_recode_chrvec = rownames(divnet_analysis_2$std_alpha$est), 
  replace_this_chrvec = rownames(divnet_analysis_2$std_alpha$rare@sam_data), 
  with_this_chrvec = divnet_analysis_2$std_alpha$rare@sam_data$mother)

summary(aov(Shannon ~ mother, data = divnet_analysis_2$std_alpha$est))
DescTools::PostHocTest(aov(Shannon ~ mother, data = divnet_analysis_2$std_alpha$est), method = c("hsd"), conf.level = 0.95)

kruskal.test(Fisher ~ mother, data = divnet_analysis_2$std_alpha$est)

phyloseq::plot_richness(physeq = divnet_analysis_2$std_alpha$rare, x = 'mother', color = 'mother', measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) + ggplot2::geom_boxplot()



divnet_analysis_2$datasets_asv_taxa$dada2

divnet_analysis_2$deicode_raw <- qiime2R::read_qza('/home/adrian/Desktop/qiime/full_data_server/nx_diversity/beta/deicode_dada2.DistanceMatrix.qza')

divnet_analysis_2$b_c_raw <- qiime2R::read_qza('/home/adrian/Desktop/qiime/full_data_server/nx_diversity/beta/b-c_rarefied_tmp_dada2.DistanceMatrix.qza')





phyloseq::plot_ordination(physeq = divnet_analysis_2$datasets_asv_taxa$dada2, ordination = divnet_analysis_2$deicode_raw$data, color = 'mother') + 
  stat_ellipse(aes(group = mother))

phyloseq::plot_net(physeq = divnet_analysis_2$datasets_asv_taxa$dada2, color = 'mother')

plot_clust(physeq = divnet_analysis_2$datasets_asv_taxa$dada2)

phyloseq::plot_heatmap(physeq = divnet_analysis_2$datasets_asv_taxa$dada2, sample.label = 'mother', )

phyloseq::plot_heatmap(physeq = divnet_analysis_2$datasets_asv_taxa$dada2, distance = divnet_analysis_2$deicode_raw$data, sample.label = 'mother')

test <- hclust(divnet_analysis_2$deicode_raw$data, method = "complete", members = NULL)

test$labels <- recode_values_based_on_key(
  to_recode_chrvec = test$labels, 
  replace_this_chrvec = rownames(original_dataset@sam_data), 
  with_this_chrvec = original_dataset@sam_data$mother)

plot(test)

ggplot() +
geom_segment(divnet_analysis_2$deicode_raw$data)


limma::plotMDS(x = divnet_analysis_2$deicode_raw$data)

ggplot(data = )
plot(tmp[,1], tmp[,2])
plot(x = tmp[,1], y = tmp[,2], xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="n")

factoextra::fviz_pca_var(
  cmdscale(cols$`bray-curtis`),
  repel = T,
  habillage = data.frame(original_dataset@sam_data)$mother)

### https://forum.qiime2.org/t/how-to-make-pcoa-biplot-in-r-using-q2-deicode-ordination/8377/8
### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320001/
### On the basis of our observations described here, the horseshoe effect appears in dimensionality reduction techniques due to the saturation property of distance metrics. While we have tested only a few distance metrics, it is suspected that a vast majority of these distance metrics exhibit the same property, which would also explain why horseshoes are encountered so frequently across many different fields. The saturation property has also been observed in multiple other fields, and other studies from different disciplines have led to similar conclusions (5). In spite of the saturation property of distance metrics, identifying horseshoes is still highly useful for identifying patterns concerning niche differentiation. These insights can ultimately guide additional statistical analyses, such as network analyses and indicator taxon analyses, to facilitate the targeted characterization of microbial niches.
### https://www.gdc-docs.ethz.ch/MDA/handouts/MDA20_PhyloseqFormation_Mahendra_Mariadassou.pdf
### Spatial Characteristics and Temporal Evolution of Chemical and Biological Freshwater Status as Baseline Assessment on the Tropical Island San Cristóbal - To extenuate the horseshoe effect, a chord or Hellinger transformation is often performed before running PCA
### https://mb3is.megx.net/gustame/reference/transformations
### https://forum.qiime2.org/t/horseshoe-effect-question/3571/2
### MOTHER DADA2 CENTERED ANALYSIS ###
######################################







######################################
### ADDITIONAL GROUPXMOTHER DIVNET ###
divnet_analysis_3 <- divnet_analysis_2

save(divnet_analysis_3, file = 'divnet_analysis_3.save')

### Done on calculation server ### WILL THIS WORK?
divnet_analysis_3$divnet <- DivNet::divnet(
  W = divnet_analysis$datasets_asv_taxa$dada2, 
  X = "group + mother", 
  ncores = divnet_analysis$cores, 
  tuning = divnet_analysis$tuning)



load('divnet_analysis_2_calculated.save')

divnet_analysis_2$divnet$batch_corr_mother_dada2 <- list()

divnet_analysis_2$divnet$batch_corr_mother_dada2$mother <- divnet_analysis_2_calculated$divnet
divnet_analysis_2$divnet[[1]] <- NULL

divnet_analysis_2$which_datasets_to_analyze <- names(divnet_analysis_2$datasets_asv_taxa) == 'batch_corr_mother_dada2'

rm(divnet_analysis_2_calculated)
### Done on calculation server ###

### ADDITIONAL GROUPXMOTHER DIVNET ###
######################################

# Beta diversity was assessed using the resulting Bray-Curtis distance
# matrix in a principal component analysis using the prcomp() function in R and visualized
# in a PCA using ggplot2. The distance matrix was also used to test for significant difference
# among samples by depth grouping and month using the adonis() function in the R vegan
# package with the following command: adonis(formula = bc ~ Depth + Month, data =
#                                                       795 metadata, permutations = 999).
# https://www.biorxiv.org/content/10.1101/2020.10.18.342550v2.full.pdf


### asv level is absolutely absurd when it comes to computation power needed
# divnet_analysis_2$datasets_taxa_no_sp <- rlist::list.flatten(purrr::map(
#   .x = divnet_analysis_2$datasets_taxa,
#   .f = function(tax_levels){
#     
#     tax_levels[names(tax_levels) != 'Species']
#   })
# )



# divnet_analysis_2$datasets_no_single_taxa_no_sp <- rlist::list.flatten(purrr::map(
#   .x = divnet_analysis_2$datasets_no_single_taxa,
#   .f = function(tax_levels){
#     
#     tax_levels[names(tax_levels) != 'Species']
#   })
# )
# 
# 
# divnet_analysis_2$divnet_no_single <- purrr::map(
#   .x = divnet_analysis_2$datasets_no_single_taxa_no_sp,
#   .f = function(dataset){
#     
#     purrr::map2(
#       .x = divnet_analysis_2$formulas,
#       .y = names(divnet_analysis_2$formulas),
#       .f = function(formula_, cov){
#         
#         if (cov == 'dex_muscle_v_treatment') {
#           dataset <- phyloseq::subset_samples(physeq = dataset, !is.na(dex_muscle_v_treatment))
#         }
#         
#         DivNet::divnet(W = dataset, X = cov, ncores = divnet_analysis_2$cores, tuning = divnet_analysis_2$tuning)
#       })
#   })
