devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
BiocManager::install("biomformat")
remotes::install_github("vmikk/metagMisc")
install.packages("ggpubr")
install.packages("cowplot")

save(phyloseq, file = 'phyloseq.save')
load('phyloseq.save')

phyloseq <- list(
  'dir' = 'phyloseq',
  'qualit_dir' = 'phyloseq/qualitative/',
  'deicode' = list(
    'input' = list(
      'dada2' = qiime2R::read_qza('/home/adrian/Desktop/qiime/full_data_server/nx_diversity/beta/deicode_dada2.DistanceMatrix.qza'),
      'dada2_no_singular_mother' = qiime2R::read_qza('/home/adrian/Desktop/qiime/full_data_server/no_singular_mother_dada2.DistanceMatrix.qza'),
      'deblur' = qiime2R::read_qza('/home/adrian/Desktop/qiime/full_data_server/nx_diversity/beta/deicode_deblur.DistanceMatrix.qza'))
    ),
  'datasets' = list())

dir.create(phyloseq$dir)
dir.create(phyloseq$qualit_dir)





#################
### LOAD DATA ###

phyloseq$datasets$dada2 <- list(
  'phylo' = divnet_analysis$datasets_asv_taxa$dada2,
  'divnet' = divnet_analysis$divnet_from_server$dada2,
  'alpha' = list(
    'shannon' = divnet_analysis$alpha_from_server$shannon$dada2,
    'simpson' = divnet_analysis$alpha_from_server$simpson$dada2)
  
)
phyloseq::sample_data(phyloseq$datasets$dada2$phylo)$sample_id <- rownames(phyloseq::sample_data(phyloseq$datasets$dada2$phylo))

temp_cluster_1 <- readr::read_tsv(file = 'clusters_7635_16.tsv') ### Based on PCoA of DEICODE and hclust

phyloseq::sample_data(phyloseq$datasets$dada2$phylo)$cluster <- recode_values_based_on_key(
  to_recode_chrvec = rownames(phyloseq::sample_data(phyloseq$datasets$dada2$phylo)), 
  replace_this_chrvec = temp_cluster_1$sample, 
  with_this_chrvec = temp_cluster_1$cluster)



phyloseq$datasets$mother_effect_rm_dada2 <- list(
  'phylo' = divnet_analysis$datasets_asv_taxa$batch_corr_dada2,
  'divnet' = divnet_analysis$divnet_from_server$batch_corr_dada2,
  'alpha' = list(
    'shannon' = divnet_analysis$alpha_from_server$shannon$batch_corr_dada2,
    'simpson' = divnet_analysis$alpha_from_server$simpson$batch_corr_dada2)
  
)



phyloseq$datasets$group_effect_rm_dada2 <- list(
  'phylo' = divnet_analysis_2$datasets_asv_taxa$batch_corr_mother_dada2,
  'divnet' = divnet_analysis_2$divnet$batch_corr_mother_dada2,
  'alpha' = list(
    'shannon' = divnet_analysis_2$alpha$shannon$batch_corr_mother_dada2,
    'simpson' = divnet_analysis_2$alpha$simpson$batch_corr_mother_dada2)
  
)
phyloseq::sample_data(phyloseq$datasets$group_effect_rm_dada2$phylo)$sample_id <- rownames(phyloseq::sample_data(phyloseq$datasets$group_effect_rm_dada2$phylo))

temp_cluster <- readr::read_tsv(file = 'clusters.tsv') ### Based on PCoA of DEICODE and hclust

phyloseq::sample_data(phyloseq$datasets$group_effect_rm_dada2$phylo)$cluster <- recode_values_based_on_key(
  to_recode_chrvec = rownames(phyloseq::sample_data(phyloseq$datasets$group_effect_rm_dada2$phylo)), 
  replace_this_chrvec = temp_cluster$sample, 
  with_this_chrvec = temp_cluster$cluster)



phyloseq$datasets$deblur <- list(
  'phylo' = divnet_analysis$datasets_asv_taxa$deblur,
  'divnet' = divnet_analysis$divnet_from_server$deblur,
  'alpha' = list(
    'shannon' = divnet_analysis$alpha_from_server$shannon$deblur,
    'simpson' = divnet_analysis$alpha_from_server$simpson$deblur)
  
)

### LOAD DATA ###
#################








############################################################
### BATCH CORRECTED FEATUREFREQS TO QIIME2 FOR DEICODING ###

# dir.create('additional_deicode_batch_corrected')
# 
# purrr::walk2(
#   .x = list(phyloseq$datasets$mother_effect_rm_dada2$phylo, phyloseq$datasets$group_effect_rm_dada2$phylo), 
#   .y = c('mother_effect_rm_dada2', 'group_effect_rm_dada2'), 
#   .f = function(dataset, dataset_name){
#     
#     otu_biom <- biomformat::make_biom(
#       data = as(phyloseq::otu_table(dataset),"matrix"))
#     
#     biomformat::write_biom(otu_biom, paste0('additional_deicode_batch_corrected/', dataset_name, "_otu.biom"))
#   })

# qiime tools import \
#   --input-path mother_effect_rm_dada2_otu.biom \
#   --type 'FeatureTable[Frequency]' \
#   --input-format BIOMV100Format \
#   --output-path mother_effect_rm_dada2.FeatureTable.qza
# 
# qiime tools import \
# --input-path group_effect_rm_dada2_otu.biom \
# --type 'FeatureTable[Frequency]' \
# --input-format BIOMV100Format \
# --output-path group_effect_rm_dada2.FeatureTable.qza

# qiime deicode auto-rpca \
#   --i-table mother_effect_rm_dada2.FeatureTable.qza \
#   --o-biplot mother_effect_rm_dada2.PCoAResults.qza \
#   --o-distance-matrix mother_effect_rm_dada2.DistanceMatrix.qza
# 
# qiime deicode auto-rpca \
#   --i-table group_effect_rm_dada2.FeatureTable.qza \
#   --o-biplot group_effect_rm_dada2.PCoAResults.qza \
#   --o-distance-matrix group_effect_rm_dada2.DistanceMatrix.qza

phyloseq$deicode$input$group_effect_rm_dada2 <- qiime2R::read_qza('/home/adrian/Desktop/qiime/R/additional_deicode_batch_corrected/group_effect_rm_dada2.DistanceMatrix.qza')

phyloseq$deicode$input$mother_effect_rm_dada2 <- qiime2R::read_qza('/home/adrian/Desktop/qiime/R/additional_deicode_batch_corrected/mother_effect_rm_dada2.DistanceMatrix.qza')

### BATCH CORRECTED FEATUREFREQS TO QIIME2 FOR DEICODING ###
############################################################






##################
### DRAW ALPHA ###

purrr::walk(.x = c('shannon', 'simpson'), .f = function(metric){

  purrr::walk2(
    .x = phyloseq$datasets, 
    .y = names(phyloseq$datasets), 
    .f = function(phyloseq_dataset, phyloseq_dataset_name){

      draw_alpha(
        dir_ = phyloseq$dir,
        analysis_names_ = names(phyloseq_dataset$divnet),
        dataset_name_ = phyloseq_dataset_name, 
        metric_ = metric, 
        divnet_ = phyloseq_dataset$divnet, 
        original_dataset_ = phyloseq_dataset$phylo)
  })
})
### DRAW ALPHA ###
##################






#################
### PERMANOVA ###
### NIe da się obliczyć interakcji pomiędzy mother a group bo każda kombinacja mother-group to jeden osobnik

phyloseq$deicode$adonis <- purrr::map2(
  .x = phyloseq$deicode$input, 
  .y = list(
    phyloseq$datasets$dada2$phylo, 
    phyloseq:::subset_samples(phyloseq$datasets$dada2$phylo, mother != "7635_16"),
    phyloseq$datasets$deblur$phylo,
    phyloseq$datasets$group_effect_rm_dada2$phylo,
    phyloseq$datasets$mother_effect_rm_dada2$phylo), 
  .f = function(deicode_data, phyloseq_object){
    
    adonis_ <- vegan::adonis(deicode_data$data ~ phyloseq::sample_data(phyloseq_object)$group + phyloseq::sample_data(phyloseq_object)$mother)
    
    pairwise_adonis_group <- pairwiseAdonis::pairwise.adonis(
      x = deicode_data$data, 
      factors = phyloseq::sample_data(phyloseq_object)$group, 
      p.adjust.m = 'fdr')
 
    pairwise_adonis_mother <- pairwiseAdonis::pairwise.adonis(
      x = deicode_data$data, 
      factors = phyloseq::sample_data(phyloseq_object)$mother, 
      p.adjust.m = 'fdr')
    
    permadisp_group <- anova(
      vegan::betadisper(
        d = deicode_data$data, 
        group = phyloseq::sample_data(phyloseq_object)$group)
    )
    
    permadisp_mother <- anova(
      vegan::betadisper(
        d = deicode_data$data, 
        group = phyloseq::sample_data(phyloseq_object)$mother)
    )
    
    
    return(list('adonis' = adonis_, 'pairwise_adonis_group' = pairwise_adonis_group, 'pairwise_adonis_mother' = pairwise_adonis_mother, 'permadisp_group' = permadisp_group, 'permadisp_mother' = permadisp_mother))
})

phyloseq$deicode$adonis$group_effect_rm_dada2$adonis_only_mother <- vegan::adonis(phyloseq$deicode$input$group_effect_rm_dada2$data ~ phyloseq::sample_data(phyloseq$datasets$group_effect_rm_dada2$phylo)$mother)

### PERMANOVA ###
#################





##################################
### PRINT DIVERSITY STATISTICS ###

temp_list <- list(
  'original_data_suppl_shannon' = phyloseq[["datasets"]][["dada2"]][["alpha"]][["shannon"]][["group"]][["betta_df"]],
  'original_data_suppl_simpson' = phyloseq[["datasets"]][["dada2"]][["alpha"]][["simpson"]][["group"]][["betta_df"]],
  'updated_data_mother_shannon' = phyloseq[["datasets"]][["group_effect_rm_dada2"]][["alpha"]][["shannon"]][["mother"]][["2983_15"]][["betta_df"]],
  'updated_data_mother_simpson' = phyloseq[["datasets"]][["group_effect_rm_dada2"]][["alpha"]][["simpson"]][["mother"]][["2983_15"]][["betta_df"]],
  'original_data_deicode' = phyloseq[["deicode"]][["adonis"]][["dada2"]][["adonis"]][["aov.tab"]],
  'original_data_suppl_permdisp' = phyloseq[["deicode"]][["adonis"]][["dada2"]][["permadisp_group"]],
  'original_data_mother_permdisp' = phyloseq[["deicode"]][["adonis"]][["dada2"]][["permadisp_mother"]],
  'updated_data_deicode' = phyloseq$deicode$adonis$group_effect_rm_dada2$adonis_only_mother[["aov.tab"]],
  'updated_data_permdisp' = phyloseq[["deicode"]][["adonis"]][["group_effect_rm_dada2"]][["permadisp_mother"]]
)
openxlsx::write.xlsx(x = temp_list, file = 'diversity_statistics.xlsx', rowNames = T)

### PRINT DIVERSITY STATISTICS ###
##################################







##################################
### GET  DIFFERENTIATING GENES ###

library(phyloseq)
phyloseq$datasets$group_effect_rm_dada2$phylo_filtered <- metagMisc::phyloseq_filter_prevalence(
  physeq = phyloseq$datasets$group_effect_rm_dada2$phylo, 
  prev.trh = 0.33,
  abund.trh = 10,
  threshold_condition = 'AND')

phyloseq$datasets$group_effect_rm_dada2$deseq <- purrr::map(
  .x = list('mother' = 'mother', 'cluster' = 'cluster'), 
  .f = function(factor){
    
    matrix <- phyloseq::phyloseq_to_deseq2(
      physeq = phyloseq$datasets$group_effect_rm_dada2$phylo_filtered,
      design = as.formula(paste0('~ ', factor)))

    deseq <- DESeq2::DESeq(matrix, test = "LRT", reduced = ~ 1)
    
    lrt <- DESeq2::results(deseq)
    
    results <- get_nice_results_from_lrt_deseq2(results_lrt = lrt, phyloseq_object = phyloseq[["datasets"]][["group_effect_rm_dada2"]][["phylo_filtered"]])
    
    return(list('matrix' = matrix, 'deseq' = deseq, 'lrt' = lrt, 'results' = results))
  })



phyloseq$datasets$group_effect_rm_dada2$deseq_tax_cluster <- purrr::map(
  .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'), 
  .f = function(taxon_level){
    glom_ <- phyloseq::tax_glom(physeq = phyloseq$datasets$group_effect_rm_dada2$phylo_filtered, taxrank = taxon_level)
  
    matrix_ <- phyloseq::phyloseq_to_deseq2(
      physeq = glom_,
      design = ~ cluster)
    
    deseq_ <- DESeq2::DESeq(matrix_, test = "LRT", reduced = ~ 1)
    
    lrt_ <- DESeq2::results(deseq_)
    
    result_ = lrt_[which(lrt_$padj < 0.05), ]
    result_ = cbind(as(result_, "data.frame"), as(tax_table(glom_)[rownames(result_), ], "matrix"))
    
    return('result' = result_)
})

### https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
###https://support.bioconductor.org/p/73172/#78824
### GET  DIFFERENTIATING GENES ###
##################################


###########################################
### DRAW CLUSTERS DIFFERENTIATING GENES ###




### DRAW CLUSTERS DIFFERENTIATING GENES ###
###########################################




#################
### DRAW DATA ###

temp_clustering <- hclust(
  phyloseq$deicode$input$group_effect_rm_dada2$data, 
  method = "average") ### Average would be better i think

temp_clustering_data <- ggdendro::dendro_data(temp_clustering, type = "rectangle")

temp_clustering_data$labels$sample <- temp_clustering_data$labels$label

temp_clustering_data$labels$label <- recode_values_based_on_key(
  to_recode_chrvec = temp_clustering_data$labels$label, 
  replace_this_chrvec = rownames(phyloseq$datasets$dada2$phylo@sam_data), 
  with_this_chrvec = phyloseq$datasets$dada2$phylo@sam_data$mother)

temp_clustering_data$labels$color <- recode_values_based_on_key(
  to_recode_chrvec = temp_clustering_data$labels$label, 
  replace_this_chrvec = unique(temp_clustering_data$labels$label), 
  with_this_chrvec = RColorBrewer::brewer.pal(n = 6, name = 'Set2'))

ggdendrogram_ <- ggdendro::ggdendrogram(temp_clustering_data) +
  theme(axis.text.x = element_text(hjust = 1, size = 10, face = 'bold', colour = temp_clustering_data$labels$color))




phyloseq$datasets$group_effect_rm_dada2$deseq$for_heatmap <- phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy_metadata

phyloseq$datasets$group_effect_rm_dada2$deseq$for_heatmap$frequency <-phyloseq$datasets$group_effect_rm_dada2$deseq$for_heatmap$frequency + 1

phyloseq$datasets$group_effect_rm_dada2$deseq$for_heatmap$order_in_clust <- as.factor(
  as.numeric(
    recode_values_based_on_key(
      to_recode_chrvec = phyloseq$datasets$group_effect_rm_dada2$deseq$for_heatmap$sample,
      replace_this_chrvec = temp_clustering_data$labels$sample,
      with_this_chrvec = temp_clustering_data$labels$x)))

phyloseq$datasets$group_effect_rm_dada2$deseq$for_heatmap$sample <- factor(
  x = phyloseq$datasets$group_effect_rm_dada2$deseq$for_heatmap$sample, 
  levels = temp_clustering_data$labels$sample)

heatmap_ <- ggplot(
  data = phyloseq$datasets$group_effect_rm_dada2$deseq$for_heatmap,
  mapping = aes(x = paste0(mother, "__", sample), y = full_name, fill = log(frequency))) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = 'bold'))









pcoa_ <- phyloseq::plot_ordination(
  physeq = phyloseq$datasets$group_effect_rm_dada2$phylo,
  ordination = phyloseq::ordinate(
    physeq = phyloseq$datasets$group_effect_rm_dada2$phylo, 
    method = "PCoA",
    distance = dist(phyloseq$deicode$input$group_effect_rm_dada2$data)),
  color = 'mother') +
  geom_text(mapping = aes(label = sample_id), size = 4, position = position_jitter(width = 0.2, height = 0.2)) +
  ggpubr::stat_chull(aes(fill = mother), geom = "polygon", alpha = 0.05) ### convex hull polygons; stat_ellipse(aes(group = mother))





phyloseq::plot_ordination(
  physeq = phyloseq$datasets$dada2$phylo,
  ordination = phyloseq::ordinate(
    physeq = phyloseq$datasets$dada2$phylo, 
    method = "PCoA",
    distance = dist(phyloseq$deicode$input$dada2$data)),
  color = 'mother') +
  #geom_text(mapping = aes(label = sample_id), size = 4, position = position_jitter(width = 0.2, height = 0.2)) +
  ggpubr::stat_chull(aes(fill = mother), geom = "polygon", alpha = 0.05)

  


# metagMisc::shepard_plot(
#   dis = phyloseq$deicode$input$group_effect_rm_dada2$data, 
#   ord = phyloseq::ordinate(
#     physeq = phyloseq$datasets$group_effect_rm_dada2$phylo, 
#     method = "PCoA", 
#     distance = phyloseq$deicode$input$group_effect_rm_dada2$data)
#   )

# NMDS may give stress zero error due to: When you have a few samples that are just very different from every other sample, their distances will just dominate the first axis. The result is not incorrect. Even when just inspecting the abundance table manually, it's clear that TTF0D and TTF0R are just very different. But I don't know about the Past software, so I cannot say why it creates different results. There are different implementations of the nMDS algorithm, maybe it prefilters or applies some transformation, I don't know. But euclidean distance is very different from Bray-Curtis or others, and it's generally not suited for community data without transformation. https://github.com/MadsAlbertsen/ampvis2/issues/78

### DRAW DATA ###
#################








































#######################################################################
### CHECK IF FACTORS ARE RELATED TO ANY METADATA FACTORS/COVARIATES ###
### Neither in dada2 not in both batch corrected dada2sets no correlation was found except for mother

test <- phyloseq::ordinate(
  physeq = phyloseq$datasets$mother_effect_rm_dada2$phylo, 
  method = "NMDS", 
  distance = phyloseq$deicode$input$mother_effect_rm_dada2$data)

test2 <- phyloseq::sample_data(phyloseq$datasets$mother_effect_rm_dada2$phylo)

vegan::envfit(
  ord = test$points, 
  env = as.data.frame(test2$rbc_quant_z_score),
  na.rm = TRUE)

### CHECK IF FACTORS ARE RELATED TO ANY METADATA FACTORS/COVARIATES ###
#######################################################################







######################################
### QUALITATIVE ANALYSIS - FIGURES ###

dir.create(phyloseq$qualit_dir)
phyloseq$qualititation <- NULL

phyloseq$qualititation$group_effect_rm_dada2 <- purrr::map(
  .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'),
  .f = function(tax_level) {
    
    phyloseq::tax_glom(physeq = phyloseq$datasets$group_effect_rm_dada2$phylo, taxrank = tax_level)
  })

phyloseq$qualititation$dada2 <- purrr::map(
  .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'),
  .f = function(tax_level) {
    
    phyloseq::tax_glom(physeq = phyloseq$datasets$dada2$phylo, taxrank = tax_level) 
  })



purrr::walk2(
  .x = phyloseq$qualititation$dada2,
  .y = names(phyloseq$qualititation$dada2),
  .f = function(phyloseq_object, tax_level){
    
    heat_ <- phyloseq::plot_heatmap(
      physeq = phyloseq_object,
      taxa.label = tax_level)
    
    ggplot2::ggsave(
      filename = paste0(phyloseq$qualit_dir, tax_level, '_dada2.png'), 
      plot = heat_, 
      width = 8.5, 
      height = 12.75)
    
    
    
    heat_2 <- phyloseq::plot_heatmap(
      physeq = phyloseq_object,
      taxa.label = tax_level,
      sample.order = phyloseq::sample_names(phyloseq_object))
    
    ggplot2::ggsave(
      filename = paste0(phyloseq$qualit_dir, 'original_order_', tax_level, '_dada2.png'), 
      plot = heat_2, 
      width = 8.5, 
      height = 12.75)

    
    
    purrr::walk(
      .x = c('mother', 'group'), 
      .f = function(factor_){
        
        dir.create(paste0(phyloseq$qualit_dir, factor_))
  
        heat <- phyloseq::plot_heatmap(
          physeq = phyloseq_object,
          taxa.label = tax_level,
          sample.label = factor_)

        ggplot2::ggsave(
          filename = paste0(phyloseq$qualit_dir, factor_, '/', factor_, '_', tax_level, "_", factor_, '_dada2.png'), 
          plot = heat, 
          width = 8.5, 
          height = 12.75)
        
        
        
        heat2 <- phyloseq::plot_heatmap(
          physeq = phyloseq_object,
          taxa.label = tax_level,
          sample.label = factor_,
          sample.order = factor_)
        
        ggplot2::ggsave(
          filename = paste0(phyloseq$qualit_dir, factor_, '/original_order_', factor_, '_', tax_level, "_", factor_, '_dada2.png'), 
          plot = heat2, 
          width = 8.5, 
          height = 12.75)
      })
})


### Clusterings/oridinance of both group_effect_rm_dada2 and dada2 on taxa levels do not show anything interesting for neither unifrac nor jaccard
# 'Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'
# temp <- 'Genus'
# temp_phyl <- phyloseq$qualititation$group_effect_rm_dada2[[temp]]
# temp_ord <- phyloseq::distance(physeq = temp_phyl, method = 'jaccard')
# temp_clustering <- hclust(
#   temp_ord,
#   method = "average")
# 
# temp_clustering_data <- ggdendro::dendro_data(temp_clustering, type = "rectangle")
# 
# temp_clustering_data$labels$sample <- temp_clustering_data$labels$label
# 
# temp_clustering_data$labels$label <- recode_values_based_on_key(
#   to_recode_chrvec = temp_clustering_data$labels$label,
#   replace_this_chrvec = rownames(phyloseq$datasets$group_effect_rm_dada2$phylo@sam_data),
#   with_this_chrvec = phyloseq$datasets$group_effect_rm_dada2$phylo@sam_data$mother)
# 
# temp_clustering_data$labels$color <- recode_values_based_on_key(
#   to_recode_chrvec = temp_clustering_data$labels$label,
#   replace_this_chrvec = unique(temp_clustering_data$labels$label),
#   with_this_chrvec = RColorBrewer::brewer.pal(n = 6, name = 'Set2'))
# 
# ggdendro::ggdendrogram(temp_clustering_data) +
#   theme(axis.text.x = element_text(hjust = 1, size = 10, face = 'bold', colour = temp_clustering_data$labels$color))
# 
# 
# 
# phyloseq::plot_ordination(
#   physeq = temp_phyl,
#   ordination = phyloseq::ordinate(physeq = temp_phyl, method = "DCA", distance = "bray"),
#   color = 'mother') +
# stat_ellipse(aes(group = mother))

### QUALITATIVE ANALYSIS - FIGURES ###
######################################





####################################
### QUALITATIVE ANALYSIS - TABLE ###

phyloseq$qualititation$dada2_taxon_qual <- purrr::map2(
  .x = phyloseq$qualititation$dada2,
  .y = names(phyloseq$qualititation$dada2),
  .f = function(phyloseq_object, tax_level){
    
    mean_sd_percent_by_rows_phyloseq(phyloseq_object, tax_level)
  })

openxlsx::write.xlsx(x = phyloseq$qualititation$dada2_taxon_qual, file = 'taxon_tables.xlsx')


### QUALITATIVE ANALYSIS - TABLE ###
####################################






##############################
### PLOT RAREFACTION CURVE ###

vegan::rarecurve(t(otu_table(phyloseq$datasets$dada2$phylo)), step=1000, cex=0.5)

### PLOT RAREFACTION CURVE ###
##############################




########################
### CHECK TOP TAXONS ###
# No interesting results for top taxons nor asvs based on plots below
# 'Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'
# 
temp_tax_lev <- 'Class'
temp_top <- 10
temp_phyl <- phyloseq::prune_taxa(names(sort(phyloseq::taxa_sums(phyloseq$qualititation$group_effect_rm_dada2[[temp_tax_lev]]), TRUE))[1:temp_top], phyloseq$qualititation$group_effect_rm_dada2[[temp_tax_lev]])
temp_group <- 'cluster'

temp_phyl_asv <- phyloseq::prune_taxa(names(sort(phyloseq::taxa_sums(phyloseq$datasets$group_effect_rm_dada2$phylo), TRUE))[1:500], phyloseq$datasets$group_effect_rm_dada2$phylo)

phyloseq::plot_heatmap(
  physeq = temp_phyl,
  taxa.label = temp_tax_lev)

phyloseq::plot_heatmap(
  physeq = temp_phyl,
  taxa.label = temp_tax_lev,
  sample.label = temp_group)

phyloseq::plot_heatmap(
  physeq = temp_phyl,
  taxa.label = temp_tax_lev,
  sample.order = phyloseq::sample_names(temp_phyl))

phyloseq::plot_heatmap(
  physeq = temp_phyl,
  taxa.label = temp_tax_lev,
  sample.label = temp_group,
  sample.order = 'cluster')


phyloseq::plot_ordination(
  physeq = temp_phyl,
  ordination = phyloseq::ordinate(physeq = temp_phyl, method = "DCA", distance = "bray"),
  color = temp_group) +
stat_ellipse(aes(group = !!as.symbol(temp_group)))



temp_method <- 'jaccard'
temp_method <- 'unifrac'
temp_ord <- phyloseq::distance(physeq = temp_phyl, method = temp_method)
temp_clustering <- hclust(
  temp_ord,
  method = "average")

temp_clustering_data <- ggdendro::dendro_data(temp_clustering, type = "rectangle")

temp_clustering_data$labels$sample <- temp_clustering_data$labels$label

temp_clustering_data$labels$label <- recode_values_based_on_key(
  to_recode_chrvec = temp_clustering_data$labels$label,
  replace_this_chrvec = rownames(phyloseq$datasets$group_effect_rm_dada2$phylo@sam_data),
  with_this_chrvec = phyloseq$datasets$group_effect_rm_dada2$phylo@sam_data$mother)

temp_clustering_data$labels$color <- recode_values_based_on_key(
  to_recode_chrvec = temp_clustering_data$labels$label,
  replace_this_chrvec = unique(temp_clustering_data$labels$label),
  with_this_chrvec = RColorBrewer::brewer.pal(n = 6, name = 'Set2'))

ggdendro::ggdendrogram(temp_clustering_data) +
  theme(axis.text.x = element_text(hjust = 1, size = 10, face = 'bold', colour = temp_clustering_data$labels$color))

### CHECK TOP TAXONS ###
########################




#################################################
### CHECK DEICODE DATA WITHOUT 5142_17 MOTHER ###

phyloseq$no_5142_17 <- list('phylo' = phyloseq::subset_samples(physeq = phyloseq$datasets$group_effect_rm_dada2$phylo, mother != "5142_17"))

dir.create('deicode_no_5142_17')

otu_biom <- biomformat::make_biom(data = as(phyloseq::otu_table(phyloseq$no_5142_17$phylo), "matrix"))

biomformat::write_biom(otu_biom, paste0('deicode_no_5142_17/group_effect_rm_dada2_for_no_5142_17_otu.biom'))

# qiime tools import \
# --input-path group_effect_rm_dada2_for_no_5142_17_otu.biom \
# --type 'FeatureTable[Frequency]' \
# --input-format BIOMV100Format \
# --output-path group_effect_rm_dada2_for_no_5142_17_otu.FeatureTable.qza
# 
# qiime deicode auto-rpca \
# --i-table group_effect_rm_dada2_for_no_5142_17_otu.FeatureTable.qza \
# --o-biplot group_effect_rm_dada2_for_no_5142_17_otu.PCoAResults.qza \
# --o-distance-matrix group_effect_rm_dada2_for_no_5142_17_otu.DistanceMatrix.qza

phyloseq$no_5142_17$deicode <- qiime2R::read_qza('/home/adrian/Desktop/qiime/R/deicode_no_5142_17/group_effect_rm_dada2_for_no_5142_17_otu.DistanceMatrix.qza')


phyloseq$no_5142_17$adonis <- vegan::adonis(phyloseq$no_5142_17$deicode$data ~ phyloseq::sample_data(phyloseq$no_5142_17$phylo)$group + phyloseq::sample_data(phyloseq$no_5142_17$phylo)$mother)

### CHECK DEICODE DATA WITHOUT 5142_17 MOTHER ###
#################################################














########################
### CLUSTER ANALYSIS ###
phyloseq$cluster_analysis <- list()
phyloseq$cluster_analysis$adonis <- vegan::adonis(phyloseq$deicode$input$dada2$data ~ phyloseq::sample_data(phyloseq$datasets$dada2$phylo)$cluster)

phyloseq$cluster_analysis$permadisp <- anova(
      vegan::betadisper(
        d = phyloseq$deicode$input$dada2$data, 
        group = phyloseq::sample_data(phyloseq$datasets$dada2$phylo)$cluster)
    )



library(phyloseq)
phyloseq$datasets$dada2$phylo_filtered <- metagMisc::phyloseq_filter_prevalence(
  physeq = phyloseq$datasets$dada2$phylo, 
  prev.trh = 0.33,
  abund.trh = 10,
  threshold_condition = 'AND')

phyloseq$cluster_analysis$matrix <- phyloseq::phyloseq_to_deseq2(
  physeq = phyloseq$datasets$dada2$phylo_filtered,
  design = ~ cluster)
    
phyloseq$cluster_analysis$deseq <- DESeq2::DESeq(phyloseq$cluster_analysis$matrix, test = "LRT", reduced = ~ 1)

phyloseq$cluster_analysis$lrt <- DESeq2::results(phyloseq$cluster_analysis$deseq)

phyloseq$cluster_analysis$results <- get_nice_results_from_lrt_deseq2(results_lrt = phyloseq$cluster_analysis$lrt, phyloseq_object = phyloseq[["datasets"]][["dada2"]][["phylo_filtered"]])

phyloseq[["cluster_analysis"]][["results"]][["lrt_df_signif_full"]] <- subset(
  x = phyloseq[["cluster_analysis"]][["results"]][["lrt_df_signif_full"]], 
  subset = phyloseq[["cluster_analysis"]][["results"]][["lrt_df_signif_full"]]$padj < 0.05)

phyloseq[["cluster_analysis"]][["results"]]$lrt_df_signif_full_tidy_metadata <- subset(
  x = phyloseq[["cluster_analysis"]][["results"]][["lrt_df_signif_full_tidy_metadata"]], 
  subset = phyloseq[["cluster_analysis"]][["results"]][["lrt_df_signif_full_tidy_metadata"]]$padj < 0.05)

openxlsx::write.xlsx(x = list('adonis' = phyloseq$cluster_analysis$adonis$aov.tab, 'permadisp' = phyloseq$cluster_analysis$permadisp, 'DA' = phyloseq[["cluster_analysis"]][["results"]][["lrt_df_signif_full"]]), file = 'cluster_analysis.xlsx')




phyloseq$cluster_analysis$clustering <- hclust(
  phyloseq$deicode$input$dada2$data, 
  method = "average") ### Average would be better i think

temp_clustering_data <- ggdendro::dendro_data(phyloseq$cluster_analysis$clustering, type = "rectangle")

temp_for_drawing <- phyloseq[["cluster_analysis"]][["results"]]$lrt_df_signif_full_tidy_metadata

temp_for_drawing$frequency <-temp_for_drawing$frequency + 1

temp_for_drawing$order_in_clust <- as.factor(
  as.numeric(
    recode_values_based_on_key(
      to_recode_chrvec = temp_for_drawing$sample,
      replace_this_chrvec = phyloseq$cluster_analysis$clustering$labels$sample,
      with_this_chrvec = phyloseq$cluster_analysis$clustering$labels$x)))

temp_for_drawing$sample <- factor(
  x = temp_for_drawing$sample, 
  levels = phyloseq$cluster_analysis$clustering$labels$sample)

ggplot(
  data = temp_for_drawing,
  mapping = aes(x = paste0(cluster, "__", sample), y = full_name, fill = log(frequency))) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = 'bold'))




phyloseq::plot_ordination(
  physeq = phyloseq$datasets$dada2$phylo,
  ordination = phyloseq::ordinate(physeq = phyloseq$datasets$dada2$phylo, method = "PCoA", distance = "jaccard"),
  color = 'cluster') 

# temp <- phyloseq::prune_taxa(names(sort(phyloseq::taxa_sums(phyloseq$datasets$dada2$phylo_filtered), TRUE))[1:500], phyloseq$datasets$dada2$phylo_filtered)
# 
# phyloseq::plot_heatmap(
#   physeq = temp,
#   sample.label = 'cluster',
#   sample.order = 'cluster')


### CLUSTER ANALYSIS ###
########################












































### MOże nie ma sensu zbijać heatmapy i dendrogramu w ten sam rysunek...


figure <- ggpubr::ggarrange(heatmap_, ggdendrogram_, 
                            labels = c("A", "B"),
                            align = 'v',
                            nrow = 2,
                            ncol = 2)

cowplot::ggdraw() +
  cowplot::draw_plot(heatmap_, x = 0, y = 0, width = 1, height = 0.8) +
  cowplot::draw_plot(ggdendrogram_, x = 0.62, y = 0.8, width = 0.31, height = 0.2)

cowplot::plot_grid(
  plotlist = list(heatmap_, ggdendrogram_),
  nrow = 2,
  scale = c(1, 0.4),
  hjust = c(-0.5, -0.5)
  
)
cowplot::axis_canvas(plot = ggdendrogram_, axis = 'x') +
  ggdendrogram_

cowplot::ggdraw(cowplot::insert_xaxis_grob(heatmap_, cowplot::axis_canvas(plot = ggdendrogram_, axis = 'x'), grid::unit(.07, "null"),
                                           position = "bottom"))

heatmap_
cowplot::insert_xaxis_grob(plot = heatmap_, grob = ggdendrogram_)


ggpubr::ggarrange(heatmap_, 
                  ggpubr::ggarrange(
                    ggdendrogram_,
                    width = 1),
                  nrow = 2)




phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy_metadata %>%
  ggplot(mapping = aes(x = mother, y = frequency)) +
  geom_boxplot()



root(1)

phyloseq::plot_heatmap(
  physeq = phyloseq$datasets$dada2$phylo, 
  distance = phyloseq$deicode$input$data, 
  sample.label = 'mother')



phyloseq::plot_ordination(
  physeq = phyloseq$datasets$dada2$phylo, 
  ordination = phyloseq$deicode$input$data, 
  color = 'group') 



temp_clustering <- hclust(
  divnet_analysis_2$deicode_raw$data, 
  method = "average")


temp_clustering_data$labels$label <- recode_values_based_on_key(
  to_recode_chrvec = temp_clustering_data$labels$label, 
  replace_this_chrvec = rownames(phyloseq$datasets$dada2$phylo@sam_data), 
  with_this_chrvec = phyloseq$datasets$dada2$phylo@sam_data$group)

temp_clustering_data$labels$color <- recode_values_based_on_key(
  to_recode_chrvec = temp_clustering_data$labels$label, 
  replace_this_chrvec = unique(temp_clustering_data$labels$label), 
  with_this_chrvec = RColorBrewer::brewer.pal(n = 7, name = 'Set2'))

ggdendro::ggdendrogram(temp_clustering_data) +
  theme(axis.text.x = element_text(hjust = 1, size = 10, face = 'bold', colour = temp_clustering_data$labels$color))




temp_ord <- phyloseq::ordinate(physeq = temp_phyl, method = 'MDS')

phyloseq::plot_ordination(
  physeq = temp_phyl, 
  ordination = temp_ord, 
  color = 'mother',
  type = 'samples') + 
stat_ellipse(aes(group = mother))

rank_names(temp_phyl)

temp_dist <- phyloseq::distance(physeq = temp_phyl, method = 'unifrac')





phyloseq::plot_net(physeq = divnet_analysis_2$datasets_asv_taxa$dada2, color = 'mother')

plot_clust(physeq = divnet_analysis_2$datasets_asv_taxa$dada2)

phyloseq::plot_heatmap(physeq = divnet_analysis_2$datasets_asv_taxa$dada2, sample.label = 'mother', )





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
















# phyloseq$datasets$deblur_no_single <- list(
#   'phylo' = divnet_analysis_no_single$datasets$deblur_no_single,
#   'divnet' = divnet_analysis_no_single$divnet_from_server$deblur_no_single,
#   'alpha' = list(
#     'shannon' = divnet_analysis_no_single$alpha_from_server$shannon$deblur_no_single,
#     'simpson' = divnet_analysis_no_single$alpha_from_server$simpson$deblur_no_single)
#   
# )
# 
# phyloseq$datasets$mother_effect_rm_deblur_no_single <- list(
#   'phylo' = divnet_analysis_no_single$datasets$mother_effect_rm_deblur_no_single,
#   'divnet' = divnet_analysis_no_single$divnet_from_server$mother_effect_rm_deblur_no_single,
#   'alpha' = list(
#     'shannon' = divnet_analysis_no_single$alpha_from_server$shannon$mother_effect_rm_deblur_no_single,
#     'simpson' = divnet_analysis_no_single$alpha_from_server$simpson$mother_effect_rm_deblur_no_single)
#   
# )
# 
# phyloseq$datasets$group_effect_rm_deblur_no_single <- list(
#   'phylo' = divnet_analysis_no_single$datasets$group_effect_rm_deblur_no_single,
#   'divnet' = divnet_analysis_no_single$divnet$group_effect_rm_deblur_no_single,
#   'alpha' = list(
#     'shannon' = divnet_analysis_no_single$alpha$shannon$group_effect_rm_deblur_no_single,
#     'simpson' = divnet_analysis_no_single$alpha$simpson$group_effect_rm_deblur_no_single)
#   
# )














# ############################
# ### GET FORMATED DA DATA ###
# phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df <- as.data.frame(temp_lrt_results)
# phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif <- subset(x = temp_lrt_results_df, subset = temp_lrt_results_df$pvalue < 0.05)
# 
# phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full <- merge(
#   x = phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif, 
#   y = data.frame(phyloseq[["datasets"]][["group_effect_rm_dada2"]][["phylo_filtered"]]@tax_table@.Data), 
#   by = 'row.names', 
#   all.x = T)
# 
# rownames(phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full) <- phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full$Row.names
# phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full$Row.names <- NULL
# 
# phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full <- merge(
#   x = phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full, 
#   y = data.frame(phyloseq[["datasets"]][["group_effect_rm_dada2"]][["phylo_filtered"]]@otu_table), 
#   by = 'row.names', 
#   all.x = T)
# 
# phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy <- tidyr::pivot_longer(
#   data = phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full,
#   cols = c(16:47), 
#   names_to = 'sample', 
#   values_to = 'frequency')
# 
# phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy <- dplyr::select(
#   .data = phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy,
#   Row.names, log2FoldChange, pvalue, padj, Kingdom:Species, sample, frequency)
# 
# phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy_metadata <- merge(
#   x = phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy, 
#   y = data.frame(phyloseq[["datasets"]][["group_effect_rm_dada2"]][["phylo_filtered"]]@sam_data), 
#   by.x = 'sample', 
#   by.y = 0, 
#   all.x = T)
# 
# for (row_nb in seq_along(rownames(phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy_metadata))) {
#   
#   phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy_metadata$full_name[[row_nb]] <- paste(
#     phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy_metadata[row_nb,6:12],
#     collapse = '__')
# }
# 
# phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy_metadata$full_name <- stringr::str_remove(string = phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full_tidy_metadata$full_name, pattern = '^d__')
# 
# openxlsx::write.xlsx(
#   x = subset(
#     x = phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full, 
#     subset = phyloseq$datasets$group_effect_rm_dada2$deseq$lrt_df_signif_full$pvalue <0.05),
#   file = 'diff_abund_significant_asvs.xlsx')
# 
# ### GET FORMATED DA DATA ###
# ############################

