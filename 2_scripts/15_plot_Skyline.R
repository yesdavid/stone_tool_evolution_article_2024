
library(ggplot2)

######################################
######################################
# RATES
######################################
######################################

# logfile <- beastio::readLog("2_scripts/new_xmls/done/TAXA_16_PCs_9/independent_run_1/output/FBD_skyline.log",
#                                 burnin=0.1)

selected_taxatrait_combinations <- 
  c("TAXA_87_PCs_2",
    "TAXA_87_PCs_9",
    "TAXA_87_PCs_44")

list_of_plots <- list()
birth_death_div_turn_long_list <- list()
for(i in selected_taxatrait_combinations){
  
  current_logfile_path <- paste0("2_scripts/new_xmls/done/", i, "/FBD_skyline_age_uncertainty_COMBINED.log")
  
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
    
    
    ### combine and plot
    
    rates <- rbind(deathRates_HPD, 
                   birthRates_HPD, 
                   divRates,
                   turnoverRates)
    rates$item <- as.factor(rates$item)
    rates$Event <- factor(rates$Event, levels = c("GS-2", "GI-1", "GS-1", "Holocene"))
    
    list_of_plots[[i]] <- 
      rates %>% 
      ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = Event,
                              y = med)) +
      ggplot2::geom_errorbar(ggplot2::aes(x = Event,
                                 ymin = lower, ymax = upper),
                             width = 0.1) +
      ggplot2::geom_line(ggplot2::aes(x = Event,
                             y = med,
                             group = item)) +
      ggplot2::facet_wrap(~item,
                 scales = "free_y",
                 nrow = 4) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ylab("Rate") #+
      # ggplot2::ggtitle(label = i)
    
    
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

skyline_results_plot_grob <- 
  cowplot::plot_grid(plotlist = list_of_plots,
                        nrow = 1,
                        labels = "AUTO")
skyline_results_plot_grob

ggplot2::ggsave(skyline_results_plot_grob,
                filename = file.path("3_output", "FBD_skyline_age_uncertainty_results_plot_grob.png"),
                width = 30, height = 25, units = "cm", device = "png")
ggplot2::ggsave(skyline_results_plot_grob,
                filename = file.path("3_output", "FBD_skyline_age_uncertainty_results_plot_grob.eps"),
                width = 30, height = 25, units = "cm", device = "eps")

list_of_plots$TAXA_87_PCs_2 +
  ggplot2::ggtitle(label = "")

ggplot2::ggsave(list_of_plots$TAXA_87_PCs_2 +
                  ggplot2::ggtitle(label = ""),
                filename = file.path("3_output", "FBD_skyline_age_uncertainty_results_plot_8702.png"),
                width = 10, height = 15, units = "cm", device = "png")
ggplot2::ggsave(list_of_plots$TAXA_87_PCs_2 +
                  ggplot2::ggtitle(label = ""),
                filename = file.path("3_output", "FBD_skyline_age_uncertainty_results_plot_8702.eps"),
                width = 10, height = 15, units = "cm", device = "eps")


################################################
################################################
birth_death_div_turn_long_df <- 
  do.call(rbind.data.frame, birth_death_div_turn_long_list)


# kruskal-wallis test
## within a taxa-trait set
comparing_rates_within_taxaTraitSet_kruskal <- 
birth_death_div_turn_long_df %>% 
  group_by(taxaTraits, type) %>% 
  dplyr::do(broom::tidy(kruskal.test(x = .$value,
                                     g = .$rate,
                                     p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() 
subset(comparing_rates_within_taxaTraitSet_kruskal, p.value >= 0.05)
readr::write_csv(comparing_rates_within_taxaTraitSet_kruskal,
                 file = file.path("3_output", "comparing_rates_within_taxaTraitSet_kruskal.csv"))

## between taxa-trait sets
comparing_rates_between_taxaTraitSet_kruskal <- 
birth_death_div_turn_long_df %>% 
  group_by(type, rate) %>% 
  dplyr::do(broom::tidy(kruskal.test(x = .$value,
                                     g = .$taxaTraits,
                                     p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() 
subset(comparing_rates_between_taxaTraitSet_kruskal, p.value >= 0.05)
readr::write_csv(comparing_rates_between_taxaTraitSet_kruskal,
                 file = file.path("3_output", "comparing_rates_between_taxaTraitSet_kruskal.csv"))

# pairwise wilcox-test
## comparing birth, death, turnover, and diversification between time-bin 1, 2, 3, 4
comparing_rates_within_taxaTraitSet_pairwise <- 
birth_death_div_turn_long_df %>% 
  group_by(taxaTraits, type) %>% 
  dplyr::do(broom::tidy(pairwise.wilcox.test(x = .$value,
                                             g = .$rate,
                                             p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() 
subset(comparing_rates_within_taxaTraitSet_pairwise, p.value >= 0.05)
readr::write_csv(comparing_rates_within_taxaTraitSet_pairwise,
                 file = file.path("3_output", "comparing_rates_within_taxaTraitSet_pairwise.csv"))

## comparing birth, death, turnover, and diversification between taxa-trait combinations
comparing_rates_between_taxaTraitSet_pairwise <- 
birth_death_div_turn_long_df %>% 
  group_by(type, rate) %>% 
  dplyr::do(broom::tidy(pairwise.wilcox.test(x = .$value,
                                             g = .$taxaTraits,
                                             p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() 
subset(comparing_rates_between_taxaTraitSet_pairwise, p.value >= 0.05)
readr::write_csv(comparing_rates_between_taxaTraitSet_pairwise,
                 file = file.path("3_output", "comparing_rates_between_taxaTraitSet_pairwise.csv"))




# ################################################
# ################################################
# library(RevGadgets)
# 
# logfile_rates <- 
#   logfile %>% 
#   coda::as.mcmc() %>% 
#   as.data.frame() %>% 
#   dplyr::select(c("birthRateFBD.1", "birthRateFBD.2", "birthRateFBD.3", "birthRateFBD.4",
#                   "deathRateFBD.1", "deathRateFBD.2", "deathRateFBD.3", "deathRateFBD.4", 
#                   "samplingRateFBD.1", "samplingRateFBD.2", "samplingRateFBD.3", "samplingRateFBD.4")) %>% 
#   dplyr::mutate(diversificationRateFBD.1 = birthRateFBD.1 - deathRateFBD.1,
#                 diversificationRateFBD.2 = birthRateFBD.2 - deathRateFBD.2,
#                 diversificationRateFBD.3 = birthRateFBD.3 - deathRateFBD.3,
#                 diversificationRateFBD.4 = birthRateFBD.4 - deathRateFBD.4,
#                 turnoverRateFBD.1 = deathRateFBD.1 / birthRateFBD.1,
#                 turnoverRateFBD.2 = deathRateFBD.2 / birthRateFBD.2,
#                 turnoverRateFBD.3 = deathRateFBD.3 / birthRateFBD.3,
#                 turnoverRateFBD.4 = deathRateFBD.4 / birthRateFBD.4) %>% 
#   tidyr::pivot_longer(cols = 1:ncol(.),
#                       values_to = "rate")
# 
# logfile_rates_name_type <- 
#   logfile_rates %>% 
#   dplyr::left_join(.,
#                    data.frame(name = unique(logfile_rates$name)) %>% 
#                      dplyr::mutate(., type = sapply(name, function(x){
#                        strsplit(x, split = "Rate")[[1]][1]})) %>% 
#                      dplyr::mutate(., timeBin = sapply(name, function(x){
#                        strsplit(x, split = "RateFBD.")[[1]][2]})))
# 
# 
# logfile_rates_name_type %>% 
#   ggplot2::ggplot(data = .) +
#   ggplot2::geom_boxplot(ggplot2::aes(y = (rate), x = timeBin)) + 
#   ggplot2::facet_wrap(~factor(type,
#                      levels = c("birth", "death", "sampling",
#                                 "diversification", "turnover")),
#              scales = "free",
#              ncol = 2) +
#   ggplot2::theme_bw()
