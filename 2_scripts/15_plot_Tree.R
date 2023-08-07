############################################################
############################################################
library(Momocs)
outlines_centered_scaled_subset_PCA_raw <- 
  readRDS(file = file.path("1_data",
                           "Outlines",
                           "final_subset_outlines_centered_scaled_seed1_PCA.RDS"))
############################################################
############################################################


############################################################
# convert trees to MCC trees
############################################################
converged_runs <- file.path("2_scripts", "new_xmls", "done")

tree_paths_raw <- 
list.files(converged_runs,
           pattern = "\\.trees",
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
    paste("~/Downloads/BEAST.v2.6.2.Linux/beast/bin/treeannotator -burnin 20 -heights keep", 
          file.path(getwd(),
                    i),
          file.path(getwd(),
                    paste0(strsplit(i,
                                    split = "\\.")[[1]][1],
                           ".mcc.tre")))
  # system(sys_cmd)
  cat(sys_cmd, "\n\n",
      file = file.path("run_treeannotator_all.sh"),
      append = T)
}
# system("chmod +x run_treeannotator_all.sh")

# system("sh run_treeannotator_all.sh")
all_MCC <- list.files(converged_runs,
                      pattern = "\\.mcc.tre",
                      recursive = T,
                      full.names = T)
length(all_MCC)

############################################################
############################################################
############################################################


# library(RevGadgets) # https://revbayes.github.io/tutorials/intro/revgadgets
# library(ggtree)
# library(ggplot2)


# Fossilized birth-death tree
file <- 
  all_MCC[6]

# read in the tree 
tree <- treeio::read.beast(file)

tree_data <- tree@data
tree@data$rate_median_scaled <- scales::rescale(tree@data$rate_median, to = c(0,1))
tree@data$rate_scaled <- scales::rescale(tree@data$rate, to = c(0,1))

current_outlines <- 
  subset(outlines_centered_scaled_subset_PCA_raw$fac, ARTEFACTNAME %in% tree@phylo$tip.label) 
current_outlines <- 
data.frame(current_outlines, row.names = current_outlines$ARTEFACTNAME)

######################################
# plot tree with rates as branch colors
######################################

tree_plot <-
  ggtree(tree,
         aes(color=(rate_median_scaled)),
         size = 1) +
  geom_nodelab(aes(x=branch, label=round(posterior, 4)), color = "black", vjust=-.5, size=3) +
  geom_tiplab(aes(x=branch, label=round((rate_median_scaled), 4)), vjust=-.5, size=3) +
  scale_color_continuous(low="blue", high="red",
                         name = "Scaled rate\nof evolution") +
  theme(legend.position=c(.1, .8),
        plot.title = element_text(hjust = 0.5)) +
  geom_rootedge(rootedge = 0.2) +
  # theme_tree2()+
  # geom_tiplab(align = T,
  #             linesize=.5,
  #             offset = 0.75,
  #             color = "black") +
  ggtitle(file) +
  geom_tippoint(pch = 16) +
  geom_nodepoint(pch = 15) +
  geom_tiplab(align = T,
              as_ylab = T)
tree_plot

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
  geom_point(data = subset(tree@data, rate_scaled < 1),
             aes(x = height_median+11.206,
                 y = rate_scaled)) +
  # geom_errorbar(data = subset(tree@data, rate_scaled < 1),
  #               aes(x = height_median+11.206,
  #                   y = rate_scaled)) +
  scale_x_reverse(limits = c(16, 11), 
                  expand = c(0, 0)) +
  scale_y_continuous(expand = c(0,0))


######################################
# tree plot with metric data
tree_plot + 
  geom_facet(panel = 'Artefact length (cm)', 
             data = dplyr::select(current_outlines, 
                                  ARTEFACTNAME, 
                                  length_cm), 
             geom = geom_bar, 
             mapping = aes(x = length_cm, 
                           fill = length_cm), 
             color = "black",
             orientation = 'y', 
             width = 0.8,
             stat='identity') + 
  theme_tree2(legend.position="none") + 
  geom_facet(panel = 'Latitude', 
             data = dplyr::select(current_outlines, 
                                  ARTEFACTNAME, 
                                  Lat), 
             geom = geom_point, 
             mapping = aes(x = Lat), 
             color = "black") + 
  theme_tree2(legend.position="none")


######################################
# ACE, artefact length
phytools::plotTree(tree@phylo,ftype="i")

svl<-current_outlines$Lat
names(svl) <- current_outlines$ARTEFACTNAME

fit<-phytools::fastAnc(tree@phylo,
                       svl,
                       vars=TRUE,
                       CI=TRUE)
fit$CI[1,]
range(svl)
## projection of the reconstruction onto the edges of the tree
obj<-phytools::contMap(tree@phylo,svl,plot=FALSE)
plot(obj,legend=0.7*max(phytools::nodeHeights(tree@phylo)),
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
