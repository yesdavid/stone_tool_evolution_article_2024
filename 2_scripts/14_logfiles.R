
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

      clockmodel_ULNC_nCAT <- 
        if(length(strsplit(current_model, split = "ULNC")[[1]]) == 1){
          "nCat"
        } else {
          "ULNC"
        }
      
      if(clockmodel_ULNC_nCAT == "nCat"){

        curent_model_variables_log <-
          data.frame(TAXA_PC = i,
                     Taxa = strsplit(i, split = "_")[[1]][2],
                     traits = strsplit(i, split = "_")[[1]][4],
                     current_model = current_model,
                     clockmodel = clockmodel_ULNC_nCAT,
                     medianRate = median(c(coda::as.mcmc(logfile[,"rateValues.1"], coda::as.mcmc(logfile[,"rateValues.2"])))), #use the average of the two rate categories for the nCat model
                     varianceRate = abs(median(coda::as.mcmc(logfile[,"rateValues.1"])) - median(coda::as.mcmc(logfile[,"rateValues.2"])))
                     )
        
      } else if(clockmodel_ULNC_nCAT == "ULNC"){

        curent_model_variables_log <- 
            data.frame(TAXA_PC = i,
                       Taxa = strsplit(i, split = "_")[[1]][2],
                       traits = strsplit(i, split = "_")[[1]][4],
                       current_model = current_model,
                       clockmodel = clockmodel_ULNC_nCAT,
                       medianRate = median(coda::as.mcmc(logfile[,"rate.mean"])),
                       varianceRate = median(coda::as.mcmc(logfile[,"rate.variance"]))
                       )
      }
      
    current_model_list[[current_model]] <- 
      curent_model_variables_log
    
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


##############################
# RW: Actually, thinking about it more, aren’t we showing two different things here?
#   Rate 1 – Rate 2 reflects the difference between the rate categories, whereas rate.mean is the overall rate. 
# 
# I think I would separate these two aspects of the clock model (which are both relevant): 
#   1. Rates
# - use the average of the two rate categories for the nCat model
# - use rate.mean for the UCLN clock model (as you do in this fig)
# 2. Variance
# - use medianRate1 – medianRate2 (as you do in this fig)
# - use the estimated variance / standard deviation of the UCLN clock model

done_all_combinations_results_table <- 
  readr::read_csv(file = file.path("3_output",
                                  "done_all_combinations_results_table.csv"))

#################### Rates
# nCat2
nCat2_rate_median_plot <- 
  done_all_combinations_results_table %>%
  dplyr::filter(clockmodel == "nCat") %>%
  dplyr::mutate(Taxa = factor(Taxa),
                traits = factor(traits),
                current_model = factor(current_model, levels = c("FBD_constant_rates",
                                                                 "FBD_constant_rate_age_uncertainty",
                                                                 "FBD_constant_rate_age_uncertainty_woffset",
                                                                 "FBD_skyline",
                                                                 "FBD_skyline_age_uncertainty"))) %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = Taxa,
                                  y = traits,
                                  fill = medianRate),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "Median rate", # (median(rate1, rate2))
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 1,
                      dir = "v") + 
  ggplot2::theme(legend.position="bottom")

nCat2_rate_median_plot


# ULNC
ULNC_rate_plot <- 
  done_all_combinations_results_table %>%
  dplyr::filter(clockmodel == "ULNC") %>%
  dplyr::mutate(Taxa = factor(Taxa),
                traits = factor(traits),
                current_model = factor(current_model, levels = c("ULNC_FBD_constant_rates",
                                                                 "ULNC_FBD_constant_rate_age_uncertainty",
                                                                 "ULNC_FBD_constant_rate_age_uncertainty_woffset",
                                                                 "ULNC_FBD_skyline",
                                                                 "ULNC_FBD_skyline_age_uncertainty"))) %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = Taxa,
                                  y = traits,
                                  fill = medianRate),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "Median rate", #(median(rate.mean))
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 1,
                      dir = "v") + 
  ggplot2::theme(legend.position="bottom")

ULNC_rate_plot


# ## cowplot rates
# cowplot_ULNC_NCATmeanMedian <- 
#   cowplot::plot_grid(nCat2_rate_median_plot +
#                        ggplot2::facet_wrap(~current_model,
#                                            ncol = 1) + 
#                        ggplot2::theme(legend.position="bottom"),
#                      ULNC_rate_plot +
#                        ggplot2::facet_wrap(~current_model,
#                                            ncol = 1) + 
#                        ggplot2::theme(legend.position="bottom"),
#                      labels = "AUTO",
#                      ncol = 2,
#                      byrow = T)
# 
# cowplot_ULNC_NCATmeanMedian
# 
# ggplot2::ggsave(cowplot_ULNC_NCATmeanMedian,
#                 filename = file.path("3_output", "cowplot_ULNC_NCATmeanMedian.png"),
#                 width = 20, height = 30, units = "cm", device = "png", bg = "white")
# ggplot2::ggsave(cowplot_ULNC_NCATmeanMedian,
#                 filename = file.path("3_output", "cowplot_ULNC_NCATmeanMedian.tif"),
#                 width = 20, height = 30, units = "cm", device = "tiff")



#################### Variances

# nCat2
nCat2_variance_median_plot <- 
  done_all_combinations_results_table %>%
  dplyr::filter(clockmodel == "nCat") %>%
  dplyr::mutate(Taxa = factor(Taxa),
                traits = factor(traits),
                current_model = factor(current_model, levels = c("FBD_constant_rates",
                                                                 "FBD_constant_rate_age_uncertainty",
                                                                 "FBD_constant_rate_age_uncertainty_woffset",
                                                                 "FBD_skyline",
                                                                 "FBD_skyline_age_uncertainty"))) %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = Taxa,
                                  y = traits,
                                  fill = varianceRate),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "Median variance", # (abs(median(rate1)-median(rate2)))
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 1,
                      dir = "v") + 
  ggplot2::theme(legend.position="bottom")

nCat2_variance_median_plot


# ULNC
ULNC_variance_plot <- 
  done_all_combinations_results_table %>%
  dplyr::filter(clockmodel == "ULNC") %>%
  dplyr::mutate(Taxa = factor(Taxa),
                traits = factor(traits),
                current_model = factor(current_model, levels = c("ULNC_FBD_constant_rates",
                                                                 "ULNC_FBD_constant_rate_age_uncertainty",
                                                                 "ULNC_FBD_constant_rate_age_uncertainty_woffset",
                                                                 "ULNC_FBD_skyline",
                                                                 "ULNC_FBD_skyline_age_uncertainty"))) %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = Taxa,
                                  y = traits,
                                  fill = varianceRate),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "Median variance", # (median(rate.variance))
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 1,
                      dir = "v") + 
  ggplot2::theme(legend.position="bottom")

ULNC_variance_plot


# ## cowplot variances
# cowplot_ULNC_NCAT_variance <- 
#   cowplot::plot_grid(nCat2_variance_median_plot +
#                        ggplot2::facet_wrap(~current_model,
#                                            ncol = 1) + 
#                        ggplot2::theme(legend.position="bottom"),
#                      ULNC_variance_plot +
#                        ggplot2::facet_wrap(~current_model,
#                                            ncol = 1) + 
#                        ggplot2::theme(legend.position="bottom"),
#                      labels = "AUTO",
#                      ncol = 2,
#                      byrow = T)
# 
# cowplot_ULNC_NCAT_variance
# 
# ggplot2::ggsave(cowplot_ULNC_NCAT_variance,
#                 filename = file.path("3_output", "cowplot_ULNC_NCAT_variance.png"),
#                 width = 20, height = 30, units = "cm", device = "png", bg = "white")
# ggplot2::ggsave(cowplot_ULNC_NCAT_variance,
#                 filename = file.path("3_output", "cowplot_ULNC_NCAT_variance.tif"),
#                 width = 20, height = 30, units = "cm", device = "tiff")


############################################################################################
####################### cowplot rates+variances nCat2

cowplot_NCAT_meanMedian_variance <- 
  cowplot::plot_grid(nCat2_rate_median_plot + 
                       ggplot2::theme(legend.key.size = ggplot2::unit(1, "cm")),
                   nCat2_variance_median_plot + 
                     ggplot2::theme(legend.key.size = ggplot2::unit(1, "cm")),
                   labels = "AUTO",
                   ncol = 2,
                   byrow = T)

ggplot2::ggsave(cowplot_NCAT_meanMedian_variance,
                filename = file.path("3_output", "Figures", "cowplot_NCAT_meanMedian_variance.png"),
                width = 20, height = 30, units = "cm", device = "png", bg = "white")
ggplot2::ggsave(cowplot_NCAT_meanMedian_variance,
                filename = file.path("3_output", "Figures", "cowplot_NCAT_meanMedian_variance.tif"),
                width = 20, height = 30, units = "cm", device = "tiff")

############################################################################################
####################### cowplot rates+variances ULNC

cowplot_ULNC_meanMedian_variance <- 
  cowplot::plot_grid(ULNC_rate_plot + 
                       ggplot2::theme(legend.key.size = ggplot2::unit(1, "cm")),
                   ULNC_variance_plot + 
                     ggplot2::theme(legend.key.size = ggplot2::unit(1, "cm")),
                   labels = "AUTO",
                   ncol = 2,
                   byrow = T)

ggplot2::ggsave(cowplot_ULNC_meanMedian_variance,
                filename = file.path("3_output", "Figures", "cowplot_ULNC_meanMedian_variance.png"),
                width = 20, height = 30, units = "cm", device = "png", bg = "white")
ggplot2::ggsave(cowplot_ULNC_meanMedian_variance,
                filename = file.path("3_output", "Figures", "cowplot_ULNC_meanMedian_variance.tif"),
                width = 20, height = 30, units = "cm", device = "tiff")







