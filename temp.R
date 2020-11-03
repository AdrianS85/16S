purrr::walk(
  .x = unique(metadata_analysis$metadata_tidy$variable), 
  .f = function(col_name){
    plot_data <- metadata_analysis$metadata_tidy %>%
      dplyr::filter(variable == col_name)
    
    
    plot_ <- ggplot(plot_data, aes(x = value)) +
      geom_density()
    
    ggsave(
      filename = paste0(metadata_analysis$hist, '/', col_name, '.png'), 
      plot = plot_)

    # png(filename = paste0(metadata_analysis$dir, '/', col_name))
    # plot(hist(metadata_analysis$metadata[[col_name]]))
    # dev.off()
  })
