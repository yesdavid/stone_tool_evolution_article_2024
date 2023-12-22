
library(ggplot2)
library(magrittr)

######################################
######################################
# RATES
######################################
######################################

# logfile <- beastio::readLog("2_scripts/new_xmls/done/TAXA_16_PCs_9/independent_run_1/output/ULNC_FBD_skyline.log",
#                                 burnin=0.1)

selected_taxatrait_combinations <- 
  c("TAXA_16_PCs_2",
    "TAXA_16_PCs_9",
    "TAXA_16_PCs_44",
    "TAXA_32_PCs_2",
    "TAXA_32_PCs_9",
    "TAXA_32_PCs_44",
    "TAXA_87_PCs_2",
    "TAXA_87_PCs_9",
    "TAXA_87_PCs_44")

# list_of_plots <- list()
birth_death_div_turn_long_list <- list()
rates_list <- list()
for(i in selected_taxatrait_combinations){
  
  current_logfile_path <- paste0("2_scripts/new_xmls/done/", i, "/ULNC_FBD_skyline_age_uncertainty_COMBINED.log")
  
  if(file.exists(current_logfile_path)){
    
    logfile <- beastio::readLog(current_logfile_path)
    # Episodic Diversification Analysis
    # specify the output files
    deathRates <-
      coda::as.mcmc(logfile) %>%
      beastio::getLogFileSubset(., 
                                "deathRateFBD.",
                                start = T) 
    deathRates_HPD <- 
      deathRates %>% 
      beastio::getHPDMedian() %>% 
      as.data.frame() %>%
      dplyr::mutate(., item = as.factor("Death"))
    deathRates_HPD$Event <- c("GS-2", "GI-1", "GS-1", "Holocene")
    
    
    birthRates <-
      coda::as.mcmc(logfile) %>%
      beastio::getLogFileSubset(., 
                                "birthRateFBD.",
                                start = T) 
    birthRates_HPD <- 
      birthRates %>% 
      beastio::getHPDMedian() %>% 
      as.data.frame() %>%
      dplyr::mutate(., item = as.factor("Birth"))
    birthRates_HPD$Event <- c("GS-2", "GI-1", "GS-1", "Holocene")
    
    # diversification = birth - death
    # if it goes up, the overall population of tools increases
    # if flat, there has been a constant number of tools over time; the environmental events would not have had any influence
    # The diversification rate represent the rate at which the species diversity increases
    
    logfile_diversification <-
      beastio::getLogFileSubset(coda::as.mcmc(logfile),
                                "birthRateFBD.",
                                start = T) -
      beastio::getLogFileSubset(coda::as.mcmc(logfile),
                                "deathRateFBD.",
                                start = T)
    divRates <-
      logfile_diversification %>%
      beastio::getHPDMedian() %>%
      as.data.frame(row.names = paste0("divRateFBD.", c(1:4))) %>%
      dplyr::mutate(., item = as.factor("Diversification"))
    divRates$Event <- c("GS-2", "GI-1", "GS-1", "Holocene")
    
    # turnover = death / birth
    # The turnover rate is the rate at which one species is replaced by another species due to a birth plus death event. 
    # Hence, the turnover rate represent the longevity of a species.
    
    logfile_turnover <-
      beastio::getLogFileSubset(coda::as.mcmc(logfile),
                                "deathRateFBD.",
                                start = T) /
      beastio::getLogFileSubset(coda::as.mcmc(logfile),
                                "birthRateFBD.",
                                start = T)
    turnoverRates <-
      logfile_turnover %>%
      beastio::getHPDMedian() %>%
      as.data.frame(row.names = paste0("turnoverRateFBD.", c(1:4))) %>%
      dplyr::mutate(., item = as.factor("Turnover"))
    turnoverRates$Event <- c("GS-2", "GI-1", "GS-1", "Holocene")
    
    rm(list = c("logfile"))
    gc()
    ### combine 
    rates <- rbind(deathRates_HPD, 
                   birthRates_HPD, 
                   divRates,
                   turnoverRates)
    rates$item <- as.factor(rates$item)
    rates$Event <- factor(rates$Event, levels = c("GS-2", "GI-1", "GS-1", "Holocene"))
    rates$taxa <- factor(strsplit(i, "_")[[1]][2],
                         levels = c("16","32", "87"))
    rates$traits <- factor(strsplit(i, "_")[[1]][4],
                           levels = c("2", "9", "44"))
    
    rates_list[[i]] <- 
      rates
    
    rm(list = c("deathRates_HPD", "birthRates_HPD","divRates", "turnoverRates"))
    gc()
    # data
    deathRates_long <-
      as.data.frame(deathRates) %>%
      tidyr::pivot_longer(cols = colnames(deathRates),
                          names_to = "rate") %>%
      tibble::add_column(type = "death")
    
    birthRates_long <-
      as.data.frame(birthRates) %>%
      tidyr::pivot_longer(cols = colnames(birthRates),
                          names_to = "rate") %>%
      tibble::add_column(type = "birth")
    
    logfile_diversification_long <-
      as.data.frame(logfile_diversification) %>%
      tidyr::pivot_longer(cols = colnames(logfile_diversification),
                          names_to = "rate") %>%
      tibble::add_column(type = "diversification")
    
    logfile_turnover_long <-
      as.data.frame(logfile_turnover) %>%
      tidyr::pivot_longer(cols = colnames(logfile_turnover),
                          names_to = "rate") %>%
      tibble::add_column(type = "turnover")
    
    birth_death_div_turn_long_list[[i]] <-
      rbind(deathRates_long,
            birthRates_long,
            logfile_diversification_long,
            logfile_turnover_long) %>%
      tibble::add_column(taxaTraits = i)
  }
  
}


################################################
################################################
rates_df <- 
  do.call(rbind.data.frame, rates_list)
head(rates_df)
# readr::write_csv(rates_df,
#                  file = file.path("3_output", "rates_df_Taxa_16_32_87_Traits_2_9_44_ULNC_FBD_skyline_age_uncertainty.csv"))
# rates_df <- 
#   readr::read_csv(file.path("3_output", "rates_df_Taxa_16_87_Traits_2_9_44_ULNC_FBD_skyline_age_uncertainty.csv"))
rates_df <-
  readr::read_csv(file.path("3_output", "rates_df_Taxa_16_32_87_Traits_2_9_44_ULNC_FBD_skyline_age_uncertainty.csv"))

rates_df$Event <- factor(rates_df$Event,
                         levels = c("GS-2","GI-1","GS-1","Holocene"))

death_panel <- 
  rates_df %>%
  subset(., item == "Death") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(ggplot2::aes(x = Event,
                                   y = med),
                      size = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(x = Event,
                                      ymin = lower, ymax = upper),
                         width = 0.1) +
  ggplot2::geom_line(ggplot2::aes(x = Event,
                                  y = med, 
                                  group = item)) +
  ggplot2::facet_wrap(taxa~traits,
                      # scales = "free_y",
                      nrow = 3,
                      dir = "h") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("Rate") +
  ggplot2::xlab("Greenland Stadials/Interstadials") +
  ggplot2::ggtitle(label = "Death rate") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

birth_panel <- 
  rates_df %>%
  subset(., item == "Birth") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(ggplot2::aes(x = Event,
                                   y = med),
                      size = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(x = Event,
                                      ymin = lower, ymax = upper),
                         width = 0.1) +
  ggplot2::geom_line(ggplot2::aes(x = Event,
                                  y = med, 
                                  group = item)) +
  ggplot2::facet_wrap(taxa~traits,
                      # scales = "free_y",
                      nrow = 3,
                      dir = "h") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("Rate") +
  ggplot2::xlab("Greenland Stadials/Interstadials") +
  ggplot2::ggtitle(label = "Birth rate") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

div_panel <- 
  rates_df %>%
  subset(., item == "Diversification") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(ggplot2::aes(x = Event,
                                   y = med),
                      size = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(x = Event,
                                      ymin = lower, ymax = upper),
                         width = 0.1) +
  ggplot2::geom_line(ggplot2::aes(x = Event,
                                  y = med, 
                                  group = item)) +
  ggplot2::facet_wrap(taxa~traits,
                      # scales = "free_y",
                      nrow = 3,
                      dir = "h") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("Rate") +
  ggplot2::xlab("Greenland Stadials/Interstadials") +
  ggplot2::ggtitle(label = "Diversification rate") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

turn_panel <- 
  rates_df %>%
  subset(., item == "Turnover") %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(ggplot2::aes(x = Event,
                                   y = med),
                      size = 2) +
  ggplot2::geom_errorbar(ggplot2::aes(x = Event,
                                      ymin = lower, ymax = upper),
                         width = 0.1) +
  ggplot2::geom_line(ggplot2::aes(x = Event,
                                  y = med, 
                                  group = item)) +
  ggplot2::facet_wrap(taxa~traits,
                      # scales = "free_y",
                      nrow = 3,
                      dir = "h") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::ylab("Rate") +
  ggplot2::xlab("Greenland Stadials/Interstadials") +
  ggplot2::ggtitle(label = "Turnover rate") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


skyline_results_cowplot <- 
  cowplot::plot_grid(death_panel,
                     birth_panel,
                     div_panel,
                     turn_panel,
                     nrow = 2,
                     labels = "AUTO")
skyline_results_cowplot


ggplot2::ggsave(skyline_results_cowplot,
                filename = file.path("3_output", "ULNC_FBD_skyline_age_uncertainty_results_plot_grob.png"),
                width = 35, height = 25, units = "cm", device = "png")
ggplot2::ggsave(skyline_results_cowplot,
                filename = file.path("3_output", "ULNC_FBD_skyline_age_uncertainty_results_plot_grob.eps"),
                width = 35, height = 25, units = "cm", device = "eps")


