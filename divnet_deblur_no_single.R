# install.packages(c("devtools", "RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))
# devtools::install_github('schuyler-smith/phylosmith')
# install.packages('sf')

divnet_analysis_no_single <- list(
  'import_dir' = list(
    'feat' = '/home/adrian/Desktop/qiime/full_data_server/deblur_singletons/', 
    'tree' = '/home/adrian/Desktop/qiime/full_data_server/deblur_singletons/',
    'tax' = '/home/adrian/Desktop/qiime/full_data_server/deblur_singletons/',
    'metadata' = '/home/adrian/Desktop/qiime/metadata_enhanced_phylo.tsv'),
  'dir' = 'deblur_singletons',
  'alpha_metrics' = c('shannon', 'simpson')
)

dir.create(divnet_analysis_no_single$dir)


temp <- qiime2R::read_qza(file = paste0(divnet_analysis_no_single$import_dir$feat, 'deblur_singletons_rarefied.FeatureTableFrequency.qza'))

temp_2 <- data.frame(temp$data)

temp_2 <- temp_2[c("A1", "A2", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B6", "C2", "C3", "C4", "C5", "C6", "D1", "D2", "D3", "D4", "D5", "D6", "E1", "E2", "E3", "E5", "F1", "F3", "F4", "F6", "G1", "G3", "G5", "G6")]

otu_biom <- biomformat::make_biom(data = temp_2)

biomformat::write_biom(otu_biom, paste0("deblur_singletons_rarefied_otu_for_correcting_sample_order.biom"))

# qiime tools import \
#   --input-path deblur_singletons_rarefied_otu_for_correcting_sample_order.biom \
#   --type 'FeatureTable[Frequency]' \
#   --input-format BIOMV100Format \
#   --output-path deblur_singletons_rarefied_otu_for_corrected_sample_order.FeatureTable.qza



divnet_analysis_no_single$datasets <- list(
  'deblur_no_single' = qiime2R::qza_to_phyloseq(
    features = paste0(divnet_analysis_no_single$import_dir$feat, 'deblur_singletons.FeatureTableFrequency.qza'),
    tree = paste0(divnet_analysis_no_single$import_dir$tree, 'sepp_deblur_singletons.PhylogenyRooted.qza'),
    taxonomy = paste0(divnet_analysis_no_single$import_dir$tax, 'for_divnet_vsearch_deblur_singletons.FeatureDataTaxonomy.qza'),
    metadata = divnet_analysis_no_single$import_dir$metadata)
)

divnet_analysis_no_single$datasets$mother_eff_rm_deblur_no_single <- phyloseq::phyloseq(
  phyloseq::otu_table(
    MMUPHin::adjust_batch(
      feature_abd = divnet_analysis_no_single$datasets$deblur_no_single@otu_table, 
      batch = 'mother', 
      data = data.frame(divnet_analysis_no_single$datasets$deblur_no_single@sam_data),
      control = list('zero_inflation' = T, 'pseudo_count' = NULL, 'diagnostic_plot' = paste0(divnet_analysis_no_single$dir, '/mother_eff_rm.pdf')))[['feature_abd_adj']], 
    taxa_are_rows = T), 
  phyloseq::phy_tree(divnet_analysis_no_single$datasets$deblur_no_single@phy_tree), 
  phyloseq::tax_table(divnet_analysis_no_single$datasets$deblur_no_single@tax_table), 
  phyloseq::sample_data(divnet_analysis_no_single$datasets$deblur_no_single@sam_data))



divnet_analysis_no_single$datasets$group_eff_rm_deblur_no_single <-   phyloseq::phyloseq(
  phyloseq::otu_table(
    MMUPHin::adjust_batch(
      feature_abd = divnet_analysis_no_single$datasets$deblur_no_single@otu_table, 
      batch = 'group', 
      data = data.frame(divnet_analysis_no_single$datasets$deblur_no_single@sam_data),
      control = list('zero_inflation' = T, 'pseudo_count' = NULL, 'diagnostic_plot' = paste0(divnet_analysis_no_single$dir, '/group_eff_rm.pdf')))[['feature_abd_adj']], 
    taxa_are_rows = T), 
  phyloseq::phy_tree(divnet_analysis_no_single$datasets$deblur_no_single@phy_tree), 
  phyloseq::tax_table(divnet_analysis_no_single$datasets$deblur_no_single@tax_table), 
  phyloseq::sample_data(divnet_analysis_no_single$datasets$deblur_no_single@sam_data))

divnet_analysis_no_single$datasets$group_eff_rm_deblur_no_single <- phyloseq:::subset_samples(divnet_analysis_no_single$datasets$group_eff_rm_deblur_no_single, mother != "7635_16") 




divnet_analysis_no_single$datasets$deblur_no_single_rare <- phyloseq::rarefy_even_depth(physeq = divnet_analysis_no_single$datasets$deblur_no_single, sample.size = 5000, rngseed = 1)

divnet_analysis_no_single$datasets$mother_eff_rm_deblur_no_single_rare <- phyloseq::rarefy_even_depth(physeq = divnet_analysis_no_single$datasets$mother_eff_rm_deblur_no_single, sample.size = 5000, rngseed = 1)

divnet_analysis_no_single$datasets$group_eff_rm_deblur_no_single_rare <- phyloseq::rarefy_even_depth(physeq = divnet_analysis_no_single$datasets$group_eff_rm_deblur_no_single, sample.size = 5000, rngseed = 1)





save(divnet_analysis_no_single, file = 'divnet_analysis_no_single.save')

### To do on server
### load('divnet_analysis_no_single.save')

divnet_analysis_no_single$divnet_from_server <- purrr::map2(
  .x = divnet_analysis_no_single$datasets,
  .y = c("group", "group", "mother"),
  .f = function(dataset, cov){
    
        DivNet::divnet(W = dataset, X = cov, ncores = 28, tuning = 'fast')
})
### To do on server

divnet_analysis_no_single$divnet_from_server$divnet_group_eff_rm_deblur_no_single_rare <- DivNet::divnet(W = divnet_analysis_no_single$datasets$group_eff_rm_deblur_no_single_rare, X = 'mother', ncores = 28, tuning = 'fast')

divnet_analysis_no_single$datasets$mother_eff_rm_deblur_no_single





divnet_analysis_no_single$alpha_from_server <- purrr::map(
  .x = divnet_analysis_no_single$alpha_metrics,
  .f = function(metric){
    
    purrr::pmap(.l = list(
      divnet_analysis_no_single$divnet_from_server,
      divnet_analysis_no_single$datasets,
      names(divnet_analysis_no_single$datasets)
    ),
    .f = function(divnet_dataset, original_dataset, dataset_name){
      
      purrr::map2(
        .x = divnet_dataset,
        .y = names(divnet_dataset),
        .f = function(col, col_name){

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
names(divnet_analysis_no_single$alpha_from_server) <- divnet_analysis_no_single$alpha_metrics





divnet_analysis_no_single$datasets_taxa <- purrr::map(
  .x = divnet_analysis_no_single$datasets,
  .f = function(dataset) {
    
    purrr::map(
      .x = list('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family', 'Genus' = 'Genus', 'Species' = 'Species'),
      .f = function(tax_level) {
        
        phyloseq::tax_glom(physeq = dataset, taxrank = tax_level) 
      })
  })
divnet_analysis_no_single$datasets_asv_taxa <- rlist::list.flatten(list(divnet_analysis_no_single$datasets, divnet_analysis_no_single$datasets_taxa))