############################################################
############################################################

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

# library(RevGadgets) # https://revbayes.github.io/tutorials/intro/revgadgets
# library(ggtree)
# library(ggplot2)

tree_data_list <- list()
for(current_set in list.files(file.path("2_scripts", "new_xmls", "done"))){
  if("FBD_skyline_age_uncertainty_COMBINED.mcc.tre" %in% list.files(file.path("2_scripts", "new_xmls", "done", current_set))){
    file <- file.path("2_scripts", "new_xmls", "done", current_set, "FBD_skyline_age_uncertainty_COMBINED.mcc.tre")
    tree <- treeio::read.beast(file)
    
    tree@data$data_set <- current_set
    tree@data$taxa <- strsplit(current_set, split = "_")[[1]][2]
    tree@data$traits <- strsplit(current_set, split = "_")[[1]][4]
    tree@data$height_0.95_HPD_absolute <- unlist(lapply(tree@data$height_0.95_HPD, function(x){x[2]-x[1]}))
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
                   group = traits,
                   fill = taxa)) +
  facet_wrap(~taxa,
             nrow = 4,
             scale = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Posterior clade probability")+
  xlab("Traits")

plot_rate_median <- 
  ggplot(data = tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = rate_median,
                   group = traits,
                   fill = taxa)) +
  facet_wrap(~taxa,
             nrow = 4,
             scale = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Branch rate (median)")+
  xlab("Traits")

plot_length_median <- 
  ggplot(data = tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = length_median,
                   group = traits,
                   fill = taxa)) +
  facet_wrap(~taxa,
             nrow = 4,
             scale = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Branch length (median)") +
  xlab("Traits")

plot_height_median <- 
  ggplot(data = tree_data_df)+
  geom_boxplot(aes(x = traits,
                   y = height_0.95_HPD_absolute,
                   group = traits,
                   fill = taxa)) +
  facet_wrap(~taxa,
             nrow = 4,
             scale = "free_y") +
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

###################
# pairwise wilcox-test

# my_comparisons <- list(c("2", "3"), c("3", "6"), c("6", "9"), c("9", "10"), c("10", "20"))
# ggpubr::ggboxplot(#subset(
#   tree_data_df, #taxa == 16), 
#                   x = "traits", y = "posterior",
#           color = "traits", palette = "jco")+ 
#   ggpubr::stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   ggpubr::stat_compare_means(label.y = 1.50) +    # Add global p-value
#   facet_wrap(~taxa,
#              nrow = 4,
#              scale = "free_y")

## posterior
wilcox_posterior <- 
tree_data_df %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::do(broom::tidy(pairwise.wilcox.test(x = .$posterior,
                                             g = .$traits,
                                             p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() %>% 
  tibble::add_column(response = "posterior")

## length_median
wilcox_length_median <- 
  tree_data_df %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::do(broom::tidy(pairwise.wilcox.test(x = .$length_median,
                                             g = .$traits,
                                             p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() %>% 
  tibble::add_column(response = "length_median")

## rate_median
wilcox_rate_median <- 
  tree_data_df %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::do(broom::tidy(pairwise.wilcox.test(x = .$rate_median,
                                             g = .$traits,
                                             p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() %>% 
  tibble::add_column(response = "rate_median")

## height_0.95_HPD_absolute
wilcox_height_0.95_HPD_absolute <- 
  tree_data_df %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::do(broom::tidy(pairwise.wilcox.test(x = .$height_0.95_HPD_absolute,
                                             g = .$traits,
                                             p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() %>% 
  tibble::add_column(response = "height_median")
###################
rbind(wilcox_posterior,
      wilcox_length_median,
      wilcox_rate_median,
      wilcox_height_0.95_HPD_absolute) %>% 
  tibble::add_column(method = "Pairwise comparisons using Wilcoxon rank sum exact test",
                     adj = "P value adjustment method: BH") %>% 
  readr::write_csv(.,
                   file = file.path("3_output", "pairwise_wilcox_taxa_traits_response.csv"))


###################
# Nonparametric ANOVA: Kruskal-Wallis Test
# p<0.05 -> there is at least one group statistically different from the other groups

## posterior
kruskalWallis_posterior <- 
tree_data_df %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::do(broom::tidy(kruskal.test(x = .$posterior,
                                             g = .$traits,
                                             p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() %>% 
  tibble::add_column(response = "posterior")

## length_median
kruskalWallis_length_median <- 
tree_data_df %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::do(broom::tidy(kruskal.test(x = .$length_median,
                                     g = .$traits,
                                     p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() %>% 
  tibble::add_column(response = "length_median")

## rate_median
kruskalWallis_rate_median <- 
tree_data_df %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::do(broom::tidy(kruskal.test(x = .$rate_median,
                                     g = .$traits,
                                     p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() %>% 
  tibble::add_column(response = "rate_median")

## height_0.95_HPD_absolute
kruskalWallis_height_0.95_HPD_absolute <- 
tree_data_df %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::do(broom::tidy(kruskal.test(x = .$height_0.95_HPD_absolute,
                                     g = .$traits,
                                     p.adjust.method = "BH"))) %>% 
  dplyr::ungroup() %>% 
  tibble::add_column(response = "height_0.95_HPD_absolute")
###################
rbind(kruskalWallis_posterior,
      kruskalWallis_length_median,
      kruskalWallis_rate_median,
      kruskalWallis_height_0.95_HPD_absolute) %>% 
  readr::write_csv(.,
                   file = file.path("3_output", "kruskalWallis_taxa_traits_response.csv"))


#########################################################
#########################################################
#########################################################

taxa_traits <- "TAXA_16_PCs_9"#"TAXA_60_PCs_20"

# Fossilized birth-death tree
file <- file.path("2_scripts", "new_xmls", "done", taxa_traits, "FBD_skyline_age_uncertainty_COMBINED.mcc.tre")

# read in the tree 
tree <- treeio::read.beast(file)



split_fun <- 
  function(x, var, pos){
  unlist(lapply(x, FUN = function(x){
    strsplit(x,
             split = var)[[1]][pos]
  }))
}



tiplabel_df <- 
  data.frame(old = tree@phylo$tip.label,
             Site = split_fun(tree@phylo$tip.label,
                              var = "_",
                              pos = 4),
             NAT = split_fun(tree@phylo$tip.label,
                                  var = "_",
                                  pos = 2))

site_names_ID_vs_clean <- 
  readr::read_csv(file = file.path("1_data", "site_names_ID_vs_clean.csv"))

tiplabel_df_nice <- 
dplyr::left_join(tiplabel_df,
                 site_names_ID_vs_clean,
                 by = "Site")

tiplabel_df_nice$tip_label_new <- 
  paste0(tiplabel_df_nice$Site_nice, " (", tiplabel_df_nice$NAT, ")")

tree@phylo$tip.label <- 
  dplyr::pull(tiplabel_df_nice,
              tip_label_new)




tree_data <- tree@data
tree@data$rate_median_scaled <- scales::rescale(tree@data$rate_median, to = c(0,1))
tree@data$rate_scaled <- scales::rescale(tree@data$rate, to = c(0,1))

current_outlines <- 
  subset(outlines_centered_scaled_subset_PCA_raw$fac, ARTEFACTNAME %in% tree@phylo$tip.label) 
current_outlines <- 
data.frame(current_outlines, row.names = current_outlines$ARTEFACTNAME)

PCA_current_outlines <- as.data.frame(outlines_centered_scaled_subset_PCA_raw$x)
PCA_current_outlines$ARTEFACTNAME <- rownames(PCA_current_outlines)
PCA_current_outlines <- 
  subset(PCA_current_outlines, ARTEFACTNAME %in% tree@phylo$tip.label) 
PCA_current_outlines <- 
  data.frame(PCA_current_outlines, row.names = PCA_current_outlines$ARTEFACTNAME)

######################################
# plot tree with rates as branch colors
######################################

tree_plot <-
  ggtree::ggtree(tree,
         ggtree::aes(color=(rate_median_scaled)),
         size = 1) +
  ggtree::geom_tiplab(align = T,
                      as_ylab = T) +
  ggtree::geom_nodelab(ggtree::aes(x=branch, 
                                   label=round(posterior, 3)), 
                       color = "black", 
                       vjust=-.5, 
                       size=3) +
  ggtree::geom_range(range = "height_0.95_HPD") +
  ggplot2::scale_color_continuous(low="blue", high="red",
                         name = "Scaled clock rate") +
  ggplot2::theme(legend.position=c(.4, .8),
        plot.title = ggplot2::element_text(hjust = 0.5)) +
  # ggtree::geom_rootedge(rootedge = 0.2) +
  # ggplot2::ggtitle(file) +
  ggtree::geom_tippoint(pch = 16) +
  ggtree::geom_nodepoint(pch = 15)
tree_plot

ggplot2::ggsave(tree_plot,
                filename = file.path("3_output", paste0(taxa_traits, "_MCC_tree_plot.png")),
                width = 20, height = 30, units = "cm", device = "png")
ggplot2::ggsave(tree_plot,
                filename = file.path("3_output", paste0(taxa_traits, "_MCC_tree_plot.svg")),
                width = 20, height = 30, units = "cm", device = "svg")
ggplot2::ggsave(tree_plot,
                filename = file.path("3_output", paste0(taxa_traits, "_MCC_tree_plot.tif")),
                width = 20, height = 30, units = "cm", device = "tiff")

# plot rates as rates
ggplot() +
  ggplot2::annotate(xmin = c(16.000, 12.900), 
                    xmax = c(14.600, 11.700),
                    ymin = min(subset(tree@data, rate_scaled < 1)$rate_scaled)-0.1*max(subset(tree@data, rate_scaled < 1)$rate_scaled), 
                    ymax = max(subset(tree@data, rate_scaled < 1)$rate_scaled)+0.1*max(subset(tree@data, rate_scaled < 1)$rate_scaled),
                    "rect", 
                    alpha = 0.2, fill = c("grey60", "grey60")) +
  ggplot2::annotate("text",
                    y = 0.0005,
                    x = c(15.300, 13.750, 12.300, 11.100),
                    label = c("GS-2","GI-1","GS-1","Holocene"),
                    size = 6) +
  ggplot2::geom_point(data = subset(tree@data, rate_scaled < 1),
                      ggplot2::aes(x = height_median+11.206,
                 y = rate_scaled)) +
  # geom_errorbar(data = subset(tree@data, rate_scaled < 1),
  #               aes(x = height_median+11.206,
  #                   y = rate_scaled)) +
  ggplot2::scale_x_reverse(limits = c(16, 11), 
                  expand = c(0, 0)) +
  ggplot2::scale_y_continuous(expand = c(0,0))


######################################
# tree plot with metric data
library(ggtree)
tree_plot + 
  # geom_facet(panel = 'Artefact perimeter (cm)',
  #            data = dplyr::select(current_outlines,
  #                                 ARTEFACTNAME,
  #                                 perimeter_cm),
  #            geom = ggplot2::geom_bar,
  #            mapping = aes(x = perimeter_cm,
  #                          fill = perimeter_cm),
  #            color = "black",
  #            orientation = 'y',
  #            width = 0.8,
  #            stat='identity') +
  # theme_tree2(legend.position="none") +
  # geom_facet(panel = 'Latitude',
  #            data = dplyr::select(current_outlines,
  #                                 ARTEFACTNAME,
  #                                 Lat),
  #            geom = ggplot2::geom_bar,
  #            mapping = aes(x = Lat),
  #                       color = "black",
  #                       orientation = 'y',
  #                       width = 0.8,
  #                       stat='identity') +
  ggtree::geom_facet(panel = 'PC1', 
             data = dplyr::select(PCA_current_outlines, 
                                  ARTEFACTNAME, 
                                  PC1), 
             geom = ggplot2::geom_bar, 
             mapping = aes(x = PC1,
                           fill = PC1), 
                        color = "black",
                        orientation = 'y',
                        width = 0.8,
                        stat='identity') +
  ggtree::geom_facet(panel = 'PC2', 
             data = dplyr::select(PCA_current_outlines, 
                                  ARTEFACTNAME, 
                                  PC2), 
             geom = ggplot2::geom_bar, 
             mapping = aes(x = PC2,
                           fill = PC2), 
             color = "black",
             orientation = 'y',
             width = 0.8,
             stat='identity') +
  theme_tree2(legend.position="none")


######################################
# ACE, artefact Latitude
phytools::plotTree(tree@phylo,ftype="i",
                   type = "fan")

svl<-current_outlines$Lat
names(svl) <- current_outlines$ARTEFACTNAME

fit<-phytools::fastAnc(tree@phylo,
                       svl,
                       vars=TRUE,
                       CI=TRUE)
fit$CI[1,]
range(svl)
## projection of the reconstruction onto the edges of the tree
obj<-phytools::contMap(tree@phylo,
                       svl,
                       plot=F)
plot(obj,
     legend=0.7*max(phytools::nodeHeights(tree@phylo)),
     fsize=c(0.7,0.9))
######################################
# ACE, artefact perimeter_cm
phytools::plotTree(tree@phylo,ftype="i")

svl<-current_outlines$perimeter_cm
missing <- which(is.na(current_outlines[,"perimeter_cm"]))
names(svl) <- current_outlines$ARTEFACTNAME

fit<-phytools::fastAnc(ape::drop.tip(tree@phylo, current_outlines[missing, "ARTEFACTNAME"]),
                       svl[-missing],
                       vars=TRUE,
                       CI=TRUE)
fit$CI[1,]
range(svl)
## projection of the reconstruction onto the edges of the tree
obj<-phytools::contMap(ape::drop.tip(tree@phylo, current_outlines[missing, "ARTEFACTNAME"]),
                       svl[-missing],
                       plot=FALSE)
plot(obj,
     # legend=0.7*max(phytools::nodeHeights(tree@phylo)),
     fsize=c(0.7,0.9))

###################################### https://rfunctions.blogspot.com/2017/07/phylogenetic-comparative-methods-pcms.html
# Estimating Phylogenetic Signal 
# Blomberg's K
length_cm <- current_outlines$length_cm
names(length_cm) <- current_outlines$ARTEFACTNAME

phytools::phylosig(tree@phylo, 
                   length_cm, 
                   method="K", 
                   test=TRUE, 
                   nsim=999)

Lat <- current_outlines$Lat
names(Lat) <- current_outlines$ARTEFACTNAME

phytools::phylosig(tree@phylo, 
                   Lat, 
                   method="K", 
                   test=TRUE, 
                   nsim=999)

# Moran's I
phylotraits <-
  na.omit(phylobase::phylo4d(tree@phylo, current_outlines))
moran.test <- 
  adephylo::abouheif.moran(phylotraits,
                           method="Abouheif")

moran.test
plot(moran.test)

###################################### https://rfunctions.blogspot.com/2017/07/phylogenetic-comparative-methods-pcms.html
# Testing for Correlated Evolution
# i.e., Latitude vs artefact length
ggplot(data = current_outlines,
       aes(x = length_cm,
           y = Lat)) +
  geom_point()

cor.test(x = current_outlines$length_cm,
         y = current_outlines$Lat)

## Phylogenetic Independent Contrasts (PIC)
LatPIC <- ape::pic(Lat, tree@phylo)
length_cmPIC <- ape::pic(length_cm, tree@phylo)


OLSmodel <- lm(LatPIC ~ length_cmPIC)
plot(LatPIC ~ length_cmPIC)
abline(OLSmodel)


## Phylogenetic Generalized Least Squares (PGLS)
PGLSmodel <- nlme::gls(Lat ~ length_cm, 
                       correlation = ape::corBrownian(phy = tree@phylo), 
                       data = current_outlines)
PGLSmodel
summary(PGLSmodel)

plot(LatPIC ~ length_cmPIC)
abline(a = coef(PGLSmodel)[1], b = coef(PGLSmodel)[2])

plot(PGLSmodel, abline=c(0,0))
