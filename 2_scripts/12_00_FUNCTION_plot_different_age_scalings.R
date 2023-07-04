plot_different_age_scalings_FUN <- 
function(taxa_file_raw){

  # NO OFFSET
  plot_age_uncertainty_youngest_date_zero <- 
    ggplot2::ggplot(data = dplyr::arrange(subset_taxa, max) %>% dplyr::top_n(n=-10), 
                    ggplot2::aes(y = reorder(taxon, -max), 
                                 x=max,
                                 xmin = oneSigma_rangeMax, 
                                 xmax = oneSigma_rangeMin,
                                 color = max)) + 
    ggplot2::geom_errorbar(width = 0.5,
                           lty = "dashed") +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_reverse() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Age calBP") +
    ggplot2::ylab("Artefactnames") +
    ggplot2::geom_vline(xintercept = 0,
                        lty = "dashed",
                        color = "red") +
    ggplot2::ggtitle("Starting ages with age uncertainty AND all age ranges\n(age uncertainties cannot be <0)")
  ##################################
  
  # with OFFSET
  plot_age_uncertainty_WOFFSET <- 
    ggplot2::ggplot(data = dplyr::arrange(subset_taxa_WOFFSET, max) %>% dplyr::top_n(n=-10), 
                    ggplot2::aes(y = reorder(taxon, -max), 
                                 x= max,
                                 xmin = oneSigma_rangeMax, 
                                 xmax = oneSigma_rangeMin,
                                 color = max)) + 
    ggplot2::geom_errorbar(width = 0.5,
                           lty = "dashed") +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_reverse() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Age calBP") +
    ggplot2::ylab("Artefactnames") +
    ggplot2::geom_vline(xintercept = 0,
                        lty = "dashed",
                        color = "red") +
    ggplot2::ggtitle("Starting ages with age uncertainty\n(for offset-script)")
  ##################################
  
  ##################################
  # Youngest age fixed to 0 for SKYLINE script
  plot_age_uncertainty_youngest_date_fixed <- 
    ggplot2::ggplot(data = dplyr::arrange(subset_taxa_skyline, max) %>% dplyr::top_n(n=-10), 
                    ggplot2::aes(y = reorder(taxon, -max), 
                                 x=max,
                                 xmin = oneSigma_rangeMax, 
                                 xmax = oneSigma_rangeMin,
                                 color = max)) + 
    ggplot2::geom_errorbar(width = 0.5,
                           lty = "dashed") +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_reverse() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Age calBP") +
    ggplot2::ylab("Artefactnames") +
    ggplot2::geom_vline(xintercept = 0,
                        lty = "dashed",
                        color = "red") +
    ggplot2::ggtitle("Starting ages with youngest age fixed\n(for skyline model with age uncertainty)")
  ##################################
  
  ##################################
  cowplot::plot_grid(plot_age_uncertainty_youngest_date_zero,
                     plot_age_uncertainty_youngest_date_fixed,
                     plot_age_uncertainty_WOFFSET,
                     ncol = 1)
}
