############################################################
############################################################
library(magrittr)
############################################################
############################################################


############################################################
# convert trees to MCC trees
############################################################
converged_runs <- file.path("2_scripts", "new_xmls", "done")

tree_paths_raw <- 
  list.files(converged_runs,
             pattern = "\\COMBINED.trees",
             recursive = T,
             full.names = T)
length(tree_paths_raw)

all_runs <- 
  tree_paths_raw %>% 
  lapply(., FUN = function(x){
    strsplit(x, split = ".trees")[[1]][1]
  }) %>% 
  unlist() 

already_MCC <- 
list.files(converged_runs,
                pattern = "\\.mcc.tre",
                recursive = T,
           full.names = T) %>% 
  lapply(., FUN = function(x){
    strsplit(x, split = ".mcc")[[1]][1]
  }) %>% 
  unlist()

tree_paths <- all_runs[which(!(all_runs %in% already_MCC))]
tree_paths # trees that have not been converted to MCC trees yet

cat("#!/bin/sh\n\n",
    file = file.path("run_treeannotator_all.sh"),
    append = F)
for(i in tree_paths){
  sys_cmd <- 
    paste("~/Downloads/BEAST.v2.6.2.Linux/beast/bin/treeannotator -heights keep", 
          file.path(getwd(),
                    paste0(i,
                           ".trees")),
          file.path(getwd(),
                    paste0(strsplit(i,
                                    split = "\\.")[[1]][1],
                           ".mcc.tre")))
  # system(sys_cmd)
  cat(sys_cmd, "\n\n",
      file = file.path("run_treeannotator_all.sh"),
      append = T)
}
# system("bash run_treeannotator_all.sh")
# or exec run_treeannotator_all.sh via terminal


all_MCC <- list.files(converged_runs,
                      pattern = "\\.mcc.tre",
                      recursive = T,
                      full.names = T)
length(all_MCC)

############################################################
############################################################
############################################################
library(ggtree)
library(ggplot2)
library(Momocs)

outlines_centered_scaled_subset_PCA_raw <- 
  readRDS(file = file.path("1_data",
                           "Outlines",
                           "final_subset_outlines_centered_scaled_seed1_PCA.RDS"))
############################################################
############################################################
############################################################


############################################################
# cowplot_tree_posterior_length_rate NCAT"
############################################################
tree_data_list <- list()
for(current_set in list.files(file.path("2_scripts", "new_xmls", "done"))){
  if("FBD_skyline_age_uncertainty_COMBINED.mcc.tre" %in% list.files(file.path("2_scripts", "new_xmls", "done", current_set))){
    file <- file.path("2_scripts", "new_xmls", "done", current_set, "FBD_skyline_age_uncertainty_COMBINED.mcc.tre")
    tree <- treeio::read.beast(file)
    
    tree@data$data_set <- current_set
    tree@data$taxa <- strsplit(current_set, split = "_")[[1]][2]
    tree@data$traits <- strsplit(current_set, split = "_")[[1]][4]
    tree@data$height_0.95_HPD_absolute <- unlist(lapply(tree@data$height_0.95_HPD, function(x){abs(x[2]-x[1])}))
    tree_data_list[[current_set]] <- tree@data
  }
}
tree_data_df <- 
  do.call(rbind.data.frame, tree_data_list)
tree_data_df$taxa <- factor(tree_data_df$taxa, levels = c("16", "32", "60", "87"))
tree_data_df$traits <- factor(tree_data_df$traits, levels = c("2", "3", "6", "9", "10", "20", "44", "87"))

plot_posteriors <- 
  ggplot(data = tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = posterior,
                   group = traits#,
                   #fill = taxa
                   )) +
  facet_wrap(~taxa,
             nrow = 4#,
             # scale = "free_y"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Posterior clade probability")+
  xlab("Traits")

plot_rate_median <- 
  ggplot(data = tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = rate_median,
                   group = traits#,
                   #fill = taxa
  )) +
  facet_wrap(~taxa,
             nrow = 4#,
             # scale = "free_y"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Branch rate (median)")+
  xlab("Traits")

plot_length_median <- 
  ggplot(data = tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = length_median,
                   group = traits#,
                   #fill = taxa
  )) +
  facet_wrap(~taxa,
             nrow = 4#,
             # scale = "free_y"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Branch length (median)") +
  xlab("Traits")

plot_height_median <- 
  ggplot(data = tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = height_0.95_HPD_absolute,
                   group = traits#,
                   #fill = taxa
  )) +
  facet_wrap(~taxa,
             nrow = 4#,
             # scale = "free_y"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Length of 95% HPD age uncertainty range") +
  xlab("Traits")

cowplot_tree_posterior_length_rate <- 
cowplot::plot_grid(plot_posteriors,
                   plot_length_median,
                   plot_rate_median,
                   plot_height_median,
                   ncol = 4,
                   labels = "AUTO")
cowplot_tree_posterior_length_rate

ggsave(cowplot_tree_posterior_length_rate,
       filename = file.path("3_output", "cowplot_tree_posterior_length_rate.png"),
       width = 30, height = 20, units = "cm", device = "png")
ggsave(cowplot_tree_posterior_length_rate,
       filename = file.path("3_output", "cowplot_tree_posterior_length_rate.eps"),
       width = 30, height = 20, units = "cm", device = "eps")


############################################################
# cowplot_ULNC_tree_posterior_length_rate   ULNC
############################################################
ULNC_tree_data_list <- list()
for(current_set in list.files(file.path("2_scripts", "new_xmls", "done"))){
  if("ULNC_FBD_skyline_age_uncertainty_COMBINED.mcc.tre" %in% list.files(file.path("2_scripts", "new_xmls", "done", current_set))){
    file <- file.path("2_scripts", "new_xmls", "done", current_set, "ULNC_FBD_skyline_age_uncertainty_COMBINED.mcc.tre")
    tree <- treeio::read.beast(file)
    
    tree@data$data_set <- current_set
    tree@data$taxa <- strsplit(current_set, split = "_")[[1]][2]
    tree@data$traits <- strsplit(current_set, split = "_")[[1]][4]
    tree@data$height_0.95_HPD_absolute <- unlist(lapply(tree@data$height_0.95_HPD, function(x){abs(x[2]-x[1])}))
    ULNC_tree_data_list[[current_set]] <- tree@data
  }
}
ULNC_tree_data_df <- 
  do.call(rbind.data.frame, ULNC_tree_data_list)
ULNC_tree_data_df$taxa <- factor(ULNC_tree_data_df$taxa, levels = c("16", "32", "60", "87"))
ULNC_tree_data_df$traits <- factor(ULNC_tree_data_df$traits, levels = c("2", "3", "6", "9", "10", "20", "44", "87"))

plot_ULNC_posteriors <- 
  ggplot(data = ULNC_tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = posterior,
                   group = traits#,
                   #fill = taxa
  )) +
  facet_wrap(~taxa,
             nrow = 4#,
             # scale = "free_y"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Posterior clade probability")+
  xlab("Traits")

plot_ULNC_rate_median <- 
  ggplot(data = ULNC_tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = rate_median,
                   group = traits#,
                   #fill = taxa
  )) +
  facet_wrap(~taxa,
             nrow = 4#,
             # scale = "free_y"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Branch rate (median)")+
  xlab("Traits")

plot_ULNC_length_median <- 
  ggplot(data = ULNC_tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = length_median,
                   group = traits#,
                   #fill = taxa
  )) +
  facet_wrap(~taxa,
             nrow = 4#,
             # scale = "free_y"
             ) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Branch length (median)") +
  xlab("Traits")

plot_ULNC_height_median <- 
  ggplot(data = ULNC_tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = height_0.95_HPD_absolute,
                   group = traits#,
                   #fill = taxa
  )) +
  facet_wrap(~taxa,
             nrow = 4#,
             # scale = "free_y"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Length of 95% HPD age uncertainty range") +
  xlab("Traits")

cowplot_ULNC_tree_posterior_length_rate <- 
  cowplot::plot_grid(plot_ULNC_posteriors,
                     plot_ULNC_length_median,
                     plot_ULNC_rate_median,
                     plot_ULNC_height_median,
                     ncol = 4,
                     labels = "AUTO")
cowplot_ULNC_tree_posterior_length_rate

ggsave(cowplot_ULNC_tree_posterior_length_rate,
       filename = file.path("3_output", "cowplot_ULNC_tree_posterior_length_rate.png"),
       width = 30, height = 20, units = "cm", device = "png")
ggsave(cowplot_ULNC_tree_posterior_length_rate,
       filename = file.path("3_output", "cowplot_ULNC_tree_posterior_length_rate.eps"),
       width = 30, height = 20, units = "cm", device = "eps")


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################



variables <- c("posterior", "rate_median", "length_median", "height_0.95_HPD_absolute")

per_taxa_list <- list()
for(current_taxa in unique(ULNC_tree_data_df$taxa)){
  current_subset_taxa <- 
    ULNC_tree_data_df %>% 
    subset(., taxa == current_taxa)
  current_variables_list <- list()
  for(current_variable in variables){
    current_variables_list[[current_variable]] <- 
      data.frame(taxa = factor(current_taxa, levels = c("16", "32", "60", "87")),
                 variable = current_variable,
                 median = median(dplyr::pull(current_subset_taxa, var = current_variable),
                                 na.rm = T),
                 mean = mean(dplyr::pull(current_subset_taxa, var = current_variable),
                             na.rm = T),
                 sd = sd(dplyr::pull(current_subset_taxa, var = current_variable),
                         na.rm = T),
                 HPD_lower = coda::HPDinterval(coda::as.mcmc(current_subset_taxa[,current_variable]))[1,"lower"],
                 HPD_upper = coda::HPDinterval(coda::as.mcmc(current_subset_taxa[,current_variable]))[1,"upper"])
  }
  per_taxa_list[[current_taxa]] <-
    do.call(rbind.data.frame,current_variables_list)
}
per_taxa_df <- 
  do.call(rbind.data.frame, per_taxa_list)

per_trait_list <- list()
for(current_traits in unique(ULNC_tree_data_df$traits)){
  current_subset_trait <- 
    ULNC_tree_data_df %>% 
    subset(., traits == current_traits)
  current_variables_list <- list()
  for(current_variable in variables){
    current_variables_list[[current_variable]] <- 
      data.frame(traits = factor(current_traits, levels = c("2", "3", "6", "9", "10", "20", "44", "87")) ,
                 variable = current_variable,
                 median = median(dplyr::pull(current_subset_trait, var = current_variable),
                                 na.rm = T),
                 mean = mean(dplyr::pull(current_subset_trait, var = current_variable),
                             na.rm = T),
                 sd = sd(dplyr::pull(current_subset_trait, var = current_variable),
                         na.rm = T),
                 HPD_lower = coda::HPDinterval(coda::as.mcmc(current_subset_trait[,current_variable]))[1,"lower"],
                 HPD_upper = coda::HPDinterval(coda::as.mcmc(current_subset_trait[,current_variable]))[1,"upper"])
  }
  per_trait_list[[current_traits]] <-
    do.call(rbind.data.frame,current_variables_list)
}
per_trait_df <- 
  do.call(rbind.data.frame, per_trait_list)


per_taxa_traits_list <- list()
for(current_taxa in unique(ULNC_tree_data_df$taxa)){
  current_subset_taxa <- 
    ULNC_tree_data_df %>% 
    subset(., taxa == current_taxa)
  
  per_traits_list <- list()
  for(current_traits in unique(current_subset_taxa$traits)){
    current_subset_taxa_traits <- 
      current_subset_taxa %>% 
      subset(., traits == current_traits)
    
    current_variables_list <- list()
    for(current_variable in variables){
      current_variables_list[[current_variable]] <- 
        data.frame(taxa = factor(current_taxa, levels = c("16", "32", "60", "87")),
                   traits = factor(current_traits, levels = c("2", "3", "6", "9", "10", "20", "44", "87")),
                   variable = current_variable,
                   median = median(dplyr::pull(current_subset_taxa_traits, var = current_variable),
                                   na.rm = T),
                   mean = mean(dplyr::pull(current_subset_taxa_traits, var = current_variable),
                               na.rm = T),
                   sd = sd(dplyr::pull(current_subset_taxa_traits, var = current_variable),
                           na.rm = T),
                   HPD_lower = coda::HPDinterval(coda::as.mcmc(current_subset_taxa_traits[,current_variable]))[1,"lower"],
                   HPD_upper = coda::HPDinterval(coda::as.mcmc(current_subset_taxa_traits[,current_variable]))[1,"upper"])
      
    }
    per_traits_list[[current_traits]] <-
      do.call(rbind.data.frame, current_variables_list)
    
  }
  per_taxa_traits_list[[current_taxa]] <-
    do.call(rbind.data.frame, per_traits_list)
  
}
per_taxa_traits_df <- 
  do.call(rbind.data.frame, per_taxa_traits_list)







plot_ULNC_posteriors <- 
  per_taxa_traits_df %>% 
  subset(., variable == "posterior") %>% 
  ggplot(., ) +
  geom_point(aes(x = traits,
                 y = median)) + 
  geom_point(aes(x = traits,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = traits,
                y = median,
                group = taxa)) +
  geom_errorbar(aes(x = traits,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  facet_wrap(~taxa,
             nrow = 1) +
  theme_bw() +
  ylab("Posterior clade probability")+
  xlab("Traits")

plot_ULNC_rate_median <- 
  per_taxa_traits_df %>% 
  subset(., variable == "rate_median") %>% 
  ggplot(., ) +
  geom_point(aes(x = traits,
                 y = median)) + 
  geom_point(aes(x = traits,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = traits,
                y = median,
                group = taxa)) +
  geom_errorbar(aes(x = traits,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  facet_wrap(~taxa,
             nrow = 1) +
  theme_bw()  +
  ylab("Median branch rate")+
  xlab("Traits")

plot_ULNC_length_median <- 
  per_taxa_traits_df %>% 
  subset(., variable == "length_median") %>% 
  ggplot(., ) +
  geom_point(aes(x = traits,
                 y = median)) + 
  geom_point(aes(x = traits,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = traits,
                y = median,
                group = taxa)) +
  geom_errorbar(aes(x = traits,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  facet_wrap(~taxa,
             nrow = 1) +
  theme_bw()  +
  ylab("Median branch length") +
  xlab("Traits")

plot_ULNC_height_median <- 
  per_taxa_traits_df %>% 
  subset(., variable == "height_0.95_HPD_absolute") %>% 
  ggplot(., ) +
  geom_point(aes(x = traits,
                 y = median)) + 
  geom_point(aes(x = traits,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = traits,
                y = median,
                group = taxa)) +
  geom_errorbar(aes(x = traits,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  facet_wrap(~taxa,
             nrow = 1) +
  theme_bw() +
  ylab("Age uncertainty range length") +
  xlab("Traits")

cowplot_ULNC_tree_posterior_length_rate <- 
  cowplot::plot_grid(plot_ULNC_posteriors,
                     plot_ULNC_length_median,
                     plot_ULNC_rate_median,
                     plot_ULNC_height_median,
                     ncol = 1,
                     labels = "AUTO")
cowplot_ULNC_tree_posterior_length_rate




############ posterior

per_taxa_plot_posterior <- 
  per_taxa_df %>% 
  subset(., variable == "posterior") %>% 
  ggplot(., ) +
  geom_point(aes(x = taxa,
                 y = median)) + 
  geom_point(aes(x = taxa,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = taxa,
                y = median,
                group = variable)) +
  geom_errorbar(aes(x = taxa,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  theme_bw()  +
  theme_bw() +
  ylab("Posterior clade probability")+
  xlab("Taxa")


per_trait_plot_posterior <- 
  per_trait_df %>% 
  subset(., variable == "posterior") %>% 
  ggplot(., ) +
  geom_point(aes(x = traits,
                 y = median)) + 
  geom_point(aes(x = traits,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = traits,
                y = median,
                group = variable)) +
  geom_errorbar(aes(x = traits,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  theme_bw()  +
  theme_bw() +
  ylab("Posterior clade probability")+
  xlab("Traits")

# cowplot_posteriors <- 
#   cowplot::plot_grid(plot_ULNC_posteriors,
#                      cowplot::plot_grid(per_taxa_plot_posterior,
#                                         per_trait_plot_posterior,
#                                         ncol = 2,
#                                         labels = c("B", "C")),
#                      nrow = 2,
#                      labels = c("A"))


############ height_0.95_HPD_absolute

per_taxa_plot_height_0.95_HPD_absolute <- 
  per_taxa_df %>% 
  subset(., variable == "height_0.95_HPD_absolute") %>% 
  ggplot(., ) +
  geom_point(aes(x = taxa,
                 y = median)) + 
  geom_point(aes(x = taxa,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = taxa,
                y = median,
                group = variable)) +
  geom_errorbar(aes(x = taxa,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  theme_bw()  +
  theme_bw() +
  ylab("Age uncertainty range length")+
  xlab("Taxa")


per_trait_plot_height_0.95_HPD_absolute <- 
  per_trait_df %>% 
  subset(., variable == "height_0.95_HPD_absolute") %>% 
  ggplot(., ) +
  geom_point(aes(x = traits,
                 y = median)) + 
  geom_point(aes(x = traits,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = traits,
                y = median,
                group = variable)) +
  geom_errorbar(aes(x = traits,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  theme_bw()  +
  theme_bw() +
  ylab("Age uncertainty range length")+
  xlab("Traits")

# cowplot_height_median <- 
#   cowplot::plot_grid(plot_ULNC_height_median,
#                      cowplot::plot_grid(per_taxa_plot_height_0.95_HPD_absolute,
#                                         per_trait_plot_height_0.95_HPD_absolute,
#                                         ncol = 2,
#                                         labels = c("B", "C")),
#                      nrow = 2,
#                      labels = c("A"))

############ length_median

per_taxa_plot_length_median <- 
  per_taxa_df %>% 
  subset(., variable == "length_median") %>% 
  ggplot(., ) +
  geom_point(aes(x = taxa,
                 y = median)) + 
  geom_point(aes(x = taxa,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = taxa,
                y = median,
                group = variable)) +
  geom_errorbar(aes(x = taxa,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  theme_bw()  +
  theme_bw() +
  ylab("Median branch length")+
  xlab("Taxa")


per_trait_plot_length_median <- 
  per_trait_df %>% 
  subset(., variable == "length_median") %>% 
  ggplot(., ) +
  geom_point(aes(x = traits,
                 y = median)) + 
  geom_point(aes(x = traits,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = traits,
                y = median,
                group = variable)) +
  geom_errorbar(aes(x = traits,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  theme_bw()  +
  theme_bw() +
  ylab("Median branch length")+
  xlab("Traits")

# cowplot_length_median <- 
#   cowplot::plot_grid(plot_ULNC_length_median,
#                      cowplot::plot_grid(per_taxa_plot_length_median,
#                                         per_trait_plot_length_median,
#                                         ncol = 2,
#                                         labels = c("B", "C")),
#                      nrow = 2,
#                      labels = c("A"))

############ rate_median

per_taxa_plot_rate_median <- 
  per_taxa_df %>% 
  subset(., variable == "rate_median") %>% 
  ggplot(., ) +
  geom_point(aes(x = taxa,
                 y = median)) + 
  geom_point(aes(x = taxa,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = taxa,
                y = median,
                group = variable)) +
  geom_errorbar(aes(x = taxa,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  theme_bw()  +
  theme_bw() +
  ylab("Median branch rate")+
  xlab("Taxa")


per_trait_plot_rate_median <- 
  per_trait_df %>% 
  subset(., variable == "rate_median") %>% 
  ggplot(., ) +
  geom_point(aes(x = traits,
                 y = median)) + 
  geom_point(aes(x = traits,
                 y = mean),
             color = "red") + 
  geom_line(aes(x = traits,
                y = median,
                group = variable)) +
  geom_errorbar(aes(x = traits,
                    y = median,
                    ymin=HPD_lower, 
                    ymax=HPD_upper), 
                width=.2) +
  theme_bw()  +
  theme_bw() +
  ylab("Median branch rate")+
  xlab("Traits")

# cowplot_rate_median <- 
#   cowplot::plot_grid(plot_ULNC_rate_median,
#                      cowplot::plot_grid(per_taxa_plot_rate_median,
#                                         per_trait_plot_rate_median,
#                                         ncol = 2,
#                                         labels = c("B", "C")),
#                      nrow = 2,
#                      labels = c("A"))




#################################
# mega_cowplot_ULNC_tree_posterior_length_rate <- 
#   cowplot::plot_grid(cowplot::plot_grid(plot_ULNC_posteriors,
#                                         cowplot::plot_grid(per_taxa_plot_posterior,
#                                                            per_trait_plot_posterior,
#                                                            ncol = 2,
#                                                            labels = c("Ab", "Ac")),
#                                         nrow = 2,
#                                         labels = c("Aa")),
#                      cowplot::plot_grid(plot_ULNC_height_median,
#                                         cowplot::plot_grid(per_taxa_plot_height_0.95_HPD_absolute,
#                                                            per_trait_plot_height_0.95_HPD_absolute,
#                                                            ncol = 2,
#                                                            labels = c("Bb", "Bc")),
#                                         nrow = 2,
#                                         labels = c("Ba")),
#                      cowplot::plot_grid(plot_ULNC_length_median,
#                                         cowplot::plot_grid(per_taxa_plot_length_median,
#                                                            per_trait_plot_length_median,
#                                                            ncol = 2,
#                                                            labels = c("Cb", "Cc")),
#                                         nrow = 2,
#                                         labels = c("Ca")),
#                      cowplot::plot_grid(plot_ULNC_rate_median,
#                                         cowplot::plot_grid(per_taxa_plot_rate_median,
#                                                            per_trait_plot_rate_median,
#                                                            ncol = 2,
#                                                            labels = c("Db", "Dd")),
#                                         nrow = 2,
#                                         labels = c("Da")),
#                      ncol = 2)
# 
# mega_cowplot_ULNC_tree_posterior_length_rate

# ggsave(mega_cowplot_ULNC_tree_posterior_length_rate,
#        filename = file.path("3_output", "mega_cowplot_ULNC_tree_posterior_length_rate.png"),
#        width = 30, height = 30, units = "cm", device = "png")
# ggsave(mega_cowplot_ULNC_tree_posterior_length_rate,
#        filename = file.path("3_output", "mega_cowplot_ULNC_tree_posterior_length_rate.eps"),
#        width = 30, height = 20, units = "cm", device = "eps")

#################################
mega_cowplot_ULNC_tree_posterior_length_rate_v2 <- 
  cowplot::plot_grid(cowplot::plot_grid(plot_ULNC_posteriors,
                                        per_taxa_plot_posterior,
                                        per_trait_plot_posterior,
                                        nrow = 1,
                                        labels = c("Aa","Ab", "Ac"),
                                        rel_widths = c(2,1,1)),
                     cowplot::plot_grid(plot_ULNC_height_median,per_taxa_plot_height_0.95_HPD_absolute,
                                        per_trait_plot_height_0.95_HPD_absolute,
                                        nrow = 1,
                                        labels = c("Ba","Bb", "Bc"),
                                        rel_widths = c(2,1,1)),
                     cowplot::plot_grid(plot_ULNC_length_median,per_taxa_plot_length_median,
                                        per_trait_plot_length_median,
                                        nrow = 1,
                                        labels = c("Ca","Cb", "Cc"),
                                        rel_widths = c(2,1,1)),
                     cowplot::plot_grid(plot_ULNC_rate_median,per_taxa_plot_rate_median,
                                        per_trait_plot_rate_median,
                                        nrow = 1,
                                        labels = c("Da","Db", "Dd"),
                                        rel_widths = c(2,1,1)),
                     nrow = 4)
mega_cowplot_ULNC_tree_posterior_length_rate_v2

ggsave(mega_cowplot_ULNC_tree_posterior_length_rate_v2,
       filename = file.path("3_output", "mega_cowplot_ULNC_tree_posterior_length_rate_v2.png"),
       width = 40, height = 30, units = "cm", device = "png")
ggsave(mega_cowplot_ULNC_tree_posterior_length_rate_v2,
       filename = file.path("3_output", "mega_cowplot_ULNC_tree_posterior_length_rate_v2.eps"),
       width = 40, height = 20, units = "cm", device = "eps")

#################################
# cowplot::plot_grid(per_taxa_plot_posterior,
#                    per_taxa_plot_height_0.95_HPD_absolute,
#                    per_taxa_plot_length_median,
#                    per_taxa_plot_rate_median,
#                    per_trait_plot_posterior,
#                    per_trait_plot_height_0.95_HPD_absolute,
#                    per_trait_plot_length_median,
#                    per_trait_plot_rate_median,
#                    nrow = 2)


