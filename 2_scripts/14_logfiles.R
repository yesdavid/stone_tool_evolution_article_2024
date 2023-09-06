
library(coda)
library(dplyr)
library(magrittr)

done_output_folder_paths <- 
  list.dirs(file.path(
    "2_scripts",
    "new_xmls",
    "done"),
    recursive = F,
    full.names = F)

done_taxa_pc_combination_list <- list()
for(i in done_output_folder_paths){
  cat(i, "\n")
  current <- file.path("2_scripts",
                       "new_xmls",
                       "done",
                       i)
  
  output_combined_chains <-
    list.files(current,
               pattern = "*COMBINED.log",
               full.names = T)

  current_model_list <- list()
  for (ii in output_combined_chains) {
    
    current_model <- 
      gsub(strsplit(ii, split = "/")[[1]][length(strsplit(ii, split = "/")[[1]])],
           pattern = ".log", replacement = "")
    print(current_model)
    
    logfile <- read.table(ii, header = T)
    
    variables_log <- colnames(logfile)
    
    curent_model_variables_log_list <- list()
    for(iii in 2:length(variables_log)){
      curent_model_variables_log_list[[iii]] <- 
        data.frame(TAXA_PC = i,
                   Taxa = strsplit(i, split = "_")[[1]][2],
                   traits = strsplit(i, split = "_")[[1]][4],
                   current_model = current_model,
                   variable = variables_log[iii],
                   mean = mean(coda::as.mcmc(logfile[,variables_log[iii]])),
                   median = median(coda::as.mcmc(logfile[,variables_log[iii]])),
                   HPDinterval = coda::HPDinterval(coda::as.mcmc(logfile[,variables_log[iii]])),
                   ESS = if(length(unique(coda::as.mcmc(logfile[,variables_log[iii]])))>1){
                     coda::effectiveSize(coda::as.mcmc(logfile[,variables_log[iii]]))
                   } else {
                     NA
                   })
      # coda::traceplot(coda::as.mcmc(logfile[,variables_log[iii]]))
    }
    current_model_list[[current_model]] <- 
      do.call(rbind.data.frame, curent_model_variables_log_list)
    
    rm("logfile")
  } 
  
  done_taxa_pc_combination_list[[i]] <- 
    do.call(rbind.data.frame, current_model_list)
}

done_all_combinations_results_table <- 
  do.call(rbind.data.frame, done_taxa_pc_combination_list)

done_all_combinations_results_table$current_model <- gsub(pattern = "_COMBINED",
                                                          replacement = "",
                                                          x = done_all_combinations_results_table$current_model)

readr::write_csv(done_all_combinations_results_table,
                 file = file.path("3_output",
                                  "done_all_combinations_results_table.csv"))

###############
taxa_fun <- function(x){as.integer(strsplit(x, split = "_")[[1]][2])}
trait_fun <- function(x){as.integer(strsplit(x, split = "_")[[1]][4])}

##############################
# nCat2
nCat2_rate_results <- 
done_all_combinations_results_table %>% 
  # dplyr::filter(current_model == "FBD_constant_rates") %>%
  dplyr::filter(variable %in% c("rateValues.1", "rateValues.2")) %>% 
  dplyr::select("TAXA_PC", "variable", "mean", "median", "current_model") %>% 
  dplyr::mutate(current_model = factor(current_model, levels = c("FBD_constant_rates",
                                                                 "FBD_constant_rate_age_uncertainty",
                                                                 "FBD_constant_rate_age_uncertainty_woffset",
                                                                 "FBD_skyline",
                                                                 "FBD_skyline_age_uncertainty")),
                taxa = sapply(TAXA_PC, FUN = taxa_fun),
                traits = sapply(TAXA_PC, FUN = trait_fun)) %>% 
  dplyr::arrange(., taxa, traits) %>% 
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa)),
                traits = factor(traits, levels = unique(traits))) %>% 
  reshape(., idvar = c("TAXA_PC", "current_model"), timevar = "variable", direction = "wide") %>% 
  dplyr::mutate(difference_between_means = mean.rateValues.1 - mean.rateValues.2,
                difference_between_medians = median.rateValues.1 - median.rateValues.2)

nCat2_rate_median_plot <- 
  nCat2_rate_results %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = taxa.rateValues.2, 
                                  y = traits.rateValues.2, 
                                  fill = difference_between_medians),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "medianRate1-medianRate2",
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 2,
                      dir = "v")

ggplot2::ggsave(nCat2_rate_median_plot,
       filename = file.path("3_output", "nCat2_rate_median_plot.png"),
       width = 22, height = 20, units = "cm", device = "png")
ggplot2::ggsave(nCat2_rate_median_plot,
       filename = file.path("3_output", "nCat2_rate_median_plot.tif"),
       width = 22, height = 20, units = "cm", device = "tiff")


nCat2_rate_mean_plot <- 
  nCat2_rate_results %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = taxa.rateValues.2, 
                                  y = traits.rateValues.2, 
                                  fill = difference_between_means),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "meanRate1-meanRate2",
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 2,
                      dir = "v")

ggplot2::ggsave(nCat2_rate_mean_plot,
       filename = file.path("3_output", "nCat2_rate_mean_plot.png"),
       width = 22, height = 20, units = "cm", device = "png")
ggplot2::ggsave(nCat2_rate_mean_plot,
       filename = file.path("3_output", "nCat2_rate_mean_plot.tif"),
       width = 22, height = 20, units = "cm", device = "tiff")

##############################
# ULNC
ULNC_rate_results <- 
  done_all_combinations_results_table %>% 
  # dplyr::filter(current_model == "FBD_constant_rates") %>%
  dplyr::filter(variable %in% c("rate.mean", "rate.variance")) %>% 
  dplyr::select("TAXA_PC", "variable", "mean", "median", "current_model") %>% 
  dplyr::mutate(current_model = factor(current_model, levels = c("ULNC_FBD_constant_rates",
                                                                 "ULNC_FBD_constant_rate_age_uncertainty",
                                                                 "ULNC_FBD_constant_rate_age_uncertainty_woffset",
                                                                 "ULNC_FBD_skyline",
                                                                 "ULNC_FBD_skyline_age_uncertainty")),
                taxa = sapply(TAXA_PC, FUN = taxa_fun),
                traits = sapply(TAXA_PC, FUN = trait_fun)) %>% 
  dplyr::arrange(., taxa, traits) %>% 
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa)),
                traits = factor(traits, levels = unique(traits)))

ULNC_rate_plot <- 
  ULNC_rate_results %>% 
  tibble::as.tibble() %>% 
  subset(., variable == "rate.mean") %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = taxa, 
                                  y = traits, 
                                  fill = median),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "median(rate.mean)",
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 2,
                      dir = "v")


ggplot2::ggsave(ULNC_rate_plot,
       filename = file.path("3_output", "ULNC_rate_plot.png"),
       width = 22, height = 20, units = "cm", device = "png")
ggplot2::ggsave(ULNC_rate_plot,
       filename = file.path("3_output", "ULNC_rate_plot.tif"),
       width = 22, height = 20, units = "cm", device = "tiff")



cowplot_ULNC_NCATmeanMedian <- 
cowplot::plot_grid(nCat2_rate_mean_plot +
                     ggplot2::facet_wrap(~current_model,
                                         ncol = 1) + 
                     ggplot2::theme(legend.position="bottom"),
                   nCat2_rate_median_plot +
                     ggplot2::facet_wrap(~current_model,
                                         ncol = 1) + 
                     ggplot2::theme(legend.position="bottom"),
                   ULNC_rate_plot +
                     ggplot2::facet_wrap(~current_model,
                                         ncol = 1) + 
                     ggplot2::theme(legend.position="bottom"),
                   labels = "AUTO",
                   ncol = 3,
                   byrow = T)

ggplot2::ggsave(cowplot_ULNC_NCATmeanMedian,
                filename = file.path("3_output", "cowplot_ULNC_NCATmeanMedian.png"),
                width = 30, height = 30, units = "cm", device = "png", bg = "white")
ggplot2::ggsave(cowplot_ULNC_NCATmeanMedian,
                filename = file.path("3_output", "cowplot_ULNC_NCATmeanMedian.tif"),
                width = 30, height = 30, units = "cm", device = "tiff")

