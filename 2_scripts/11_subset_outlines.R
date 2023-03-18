library(outlineR)
library(magrittr)

#############################################################################################
#############################################################################################


single_outlines_list <- get_outlines(outpath = "/home/au656892/Documents/Doktor/2_projects/stone_tool_evolution_article_2022/1_data/Outlines/armatures_to_keep", 
                                     tps_file_rescale = NULL)

outlines_combined <- combine_outlines(single_outlines_list = single_outlines_list)


length(outlines_combined) #how many outlines do you have?
# stack(outlines_combined) # shows all outlines above one another(you might want to center and scale them first using Momocs)
# Momocs::panel(outlines_combined) # shows all outlines next to each other
# Momocs::inspect(outlines_combined) # shows only a single outline at a time. 

outlines_AR_fac_w_metric_measures <- readr::read_csv(file.path("1_data", "outlines_AR_fac_w_metric_measures.csv"))

fac <- 
subset(outlines_AR_fac_w_metric_measures, 
       ARTEFACTNAME %in% names(outlines_combined)) %>% 
  as.data.frame()
rownames(fac) <- fac$ARTEFACTNAME



outlines <- 
Momocs::Out(x = outlines_combined[[1]],
            fac = fac[match(names(outlines_combined[[1]]), 
                            fac$ARTEFACTNAME), ]  )


# only select outlines _with_ metric measurements
outlines <- Momocs::filter(outlines,
                           !is.na(SCALE))


## centre, scale, define first coordinate
outlines_centered <- Momocs::coo_centre(outlines) # Returns a shape centered on the origin.
outlines_centered_scaled <- Momocs::coo_scale(outlines_centered) # Scales the coordinates by the centroid size.
outlines_centered_scaled <- Momocs::coo_slidedirection(outlines_centered_scaled, 
                                                       direction = "up") # Sets the "first" coordinate to be on the top.




#############################################################################################
# combine with C14 dates, subset to outlines with C14 data available
#############################################################################################

filtered_C14_v4_calibrated_for_outlines <- readr::read_csv(file.path("3_output", 
                                                                     "C14",
                                                                     "filtered_C14_v4_calibrated_for_outlines.csv"))

C14_metadata <-
dplyr::left_join(outlines_centered_scaled$fac, 
                 filtered_C14_v4_calibrated_for_outlines, 
                 by = c("Site" = "Site_ID_split", 
                        "Long", "Lat")) 

table(C14_metadata$Macro_region)
table(C14_metadata$Site_strat)
table(C14_metadata$Site_type)


outlines_centered_scaled_C14 <- 
  Momocs::Out(x = outlines_centered_scaled[[1]],
              fac = C14_metadata[match(names(outlines_centered_scaled[[1]]), 
                                       C14_metadata$ARTEFACTNAME), ]  ) %>% 
  Momocs::filter(., !is.na(median_SPD_age_calBP)) 


#############################################################################################
#############################################################################################

# set.seed(1)
# subsample_stratified_n1_perSite <-
# splitstackshape::stratified(outlines_centered_scaled_C14$fac,
#                             group = "Site", # column "site" is site+layer; however: some of the unique site/layer combinations have the same 14C dates associated with them.
#                             size = 1)
set.seed(3)
subsample_stratified_n1_perSite <-
splitstackshape::stratified(outlines_centered_scaled_C14$fac,
                            group = c("TaxUnit_unique", "median_SPD_age_calBP"), # column "site" is site+layer
                            size = 1)

outlines_centered_scaled_subs <-
  Momocs::filter(outlines_centered_scaled_C14,
                 ARTEFACTNAME %in% subsample_stratified_n1_perSite$ARTEFACTNAME)

outlines_centered_scaled_subs$fac$Region <- factor(outlines_centered_scaled_subs$fac$Region)
outlines_centered_scaled_subs$fac$Expert_editor <- factor(outlines_centered_scaled_subs$fac$Expert_editor)

Momocs::panel(outlines_centered_scaled_subs,
              fac = "Region")
Momocs::panel(outlines_centered_scaled_subs,
              fac = "Expert_editor")

table(outlines_centered_scaled_subs$fac$Region)

# fileConn<-file(file.path("1_data",
#                          "Outlines",
#                          "artefactnames_subsample_stratified_n1_perSite_seed1.txt"))
# writeLines(paste0(outlines_centered_scaled_subs$fac$ARTEFACTNAME, ".jpg"), fileConn)
# close(fileConn)


#############################################################################################
#############################################################################################

# Momocs::hcontrib(outlines_centered_scaled_subs_efourier,
#                  harm.r = 1:10,
#                  amp.r=1:10,
#                  id = "TS1_Mag_BSN_BoisLaiterie_OtteStraus1997_AR_NA_2_pseudo_no_3")

# harmonic calibration
## Estimates the number of harmonics required for the Fourier methods implemented in Momocs.
outlines_centered_scaled_subs_harmonics <- Momocs::calibrate_harmonicpower_efourier(outlines_centered_scaled_subs, 
                                                                                         plot = T)  
outlines_centered_scaled_subs_harmonics$minh

# Elliptic Fourier Analysis (EFA)
outlines_centered_scaled_subs_efourier <- Momocs::efourier(outlines_centered_scaled_subs,
                                                                nb.h = as.matrix(outlines_centered_scaled_subs_harmonics[["minh"]])[[4,1]], # 4: harmonics for 99.9%
                                                                norm = F,
                                                      start = T) 

# Principal Components Analysis (PCA) on the EFA
outlines_centered_scaled_subs_PCA <- Momocs::PCA(outlines_centered_scaled_subs_efourier) # PCA on Coe objects, using prcomp.

length(unique(outlines_centered_scaled_subs_PCA$fac$Site))

saveRDS(outlines_centered_scaled_subs_PCA,
        file.path("1_data", 
                  "Outlines",
                  paste0("TAXA", nrow(outlines_centered_scaled_subs_PCA$x), 
                         "_PCs", ncol(outlines_centered_scaled_subs_PCA$x), 
                         "_outlines_centered_scaled_subset_seed1_PCA.RDS")))


#############################################################################################
#############################################################################################

## check which PC axis represents what part of the shapespace

Momocs::PCcontrib(outlines_centered_scaled_subs_PCA,
                  nax = 1:5,
                  sd.r = c(-3,-2,-1,0,1,2,3))
# 
Momocs::scree_plot(outlines_centered_scaled_subs_PCA)
# 
Momocs::plot_PCA(outlines_centered_scaled_subs_PCA,
                 f = "Region",
                 # axes = c(2,3),
                 morphospace = T,
                 labelpoints = F)

#############################################################################################
#############################################################################################

# combine FAD (first appearance date) with outlines.
# confidence interval? the dates can probably not be exactly the same for each outline due to how the FBDR process works (Rachel Warnock)
# give only the oldest date as terminus _post quem_. we don't know when they went "extinct" which is why we cannot give a LAD (last appearance date) (Cathrine Klein)

taxa_file <- data.frame(taxon = outlines_centered_scaled_subs_PCA$fac$ARTEFACTNAME,
                        max = outlines_centered_scaled_subs_PCA$fac$median_SPD_age_calBP,
                        min = 0,
                        oneSigma_rangeMax = outlines_centered_scaled_subs_PCA$fac$oneSigma_rangeMax,
                        oneSigma_rangeMin = outlines_centered_scaled_subs_PCA$fac$oneSigma_rangeMin)
readr::write_tsv(taxa_file,
                 path = file.path("1_data", "outlines_centered_scaled_subset_FAD_LAD_C14_oneSigmaMinMax.tsv"))

#############################################################################################
#############################################################################################


plot_of_selected_artefacts_and_ages_ageUncertainties <- 
  ggplot2::ggplot(data = taxa_file, 
                  aes(y = reorder(taxon, -max), 
                      x=max,
                      xmin = oneSigma_rangeMax, 
                      xmax = oneSigma_rangeMin)) + 
  ggplot2::geom_pointrange() 
plot_of_selected_artefacts_and_ages_ageUncertainties

ggsave(plot_of_selected_artefacts_and_ages_ageUncertainties,
       filename = file.path("3_output", "plot_of_selected_artefacts_and_ages_ageUncertainties.png"),
       width = 40, height = 25, units = "cm", device = "png")





