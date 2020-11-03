biplot_analysis <- list('freq_dir' = '/home/adrian/Desktop/nx_temp/nx_data_q/',
                       'freq_patt' = '^filtered.*FeatureTable.*',
                       'tax_dir' = '/home/adrian/Desktop/nx_temp/nx_taxonomy/',
                       'tax_patt' = '^vsearch.*Taxonomy\\.qza$',
                       'metadata' = qiime2R::read_q2metadata(file = '../qiime_metadata_nf.tsv'),
                       'formulas' = list('group' = '~ group'),
                       'seq_col_name' = 'seqs')





test <- qiime2R::read_qza('/home/adrian/Desktop/nx_temp/nx_differential/gneiss_phyl/phyl_filtered_sepp_dada2.FeatureTableBalance.qza')

BiocManager::install("phyloseq")
devtools::install_github("adw96/breakaway")
devtools::install_github("adw96/DivNet")


# qza_to_phyloseq() - Imports multiple artifacts to produce a phyloseq object.