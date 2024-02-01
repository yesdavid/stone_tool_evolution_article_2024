library(outlineR)
library(magrittr)


# extract outlines from manually prepared subset 
final_subset_outlines <- 
  outlineR::get_outlines(outpath = file.path("1_data",
                                             "Outlines",
                                             "final_subset_outlines")) %>% 
  outlineR::combine_outlines()

# load calibrated dates available for sites
filtered_C14_data_calibrated_for_outlines <- readr::read_csv(file.path("3_output", 
                                                                       "C14",
                                                                       "filtered_C14_data_calibrated_for_outlines.csv"))

# load meta data from 1511NAC-database
outlines_AR_fac_artefactnames_w_metric_measures <- readr::read_csv(file.path("1_data",
                                                                             "1511NAC_Database_mod",
                                                                             "outlines",
                                                                             "outlines_AR_fac_artefactnames_w_metric_measures.csv"))

meta_outlines_AR_fac_artefactnames_w_metric_measures_w_dates <- 
  dplyr::left_join(outlines_AR_fac_artefactnames_w_metric_measures,
                   filtered_C14_data_calibrated_for_outlines, 
                   by = c("Site" = "Site_ID_split", 
                          "Long", "Lat")) %>% 
  dplyr::filter(ARTEFACTNAME %in% names(final_subset_outlines$coo))

# combine outlines with all meta data
final_subset_outlines <-
  Momocs::Out(final_subset_outlines$coo,
              # order of meta_outlines_AR_fac_artefactnames_w_metric_measures_w_dates$ARTEFACTNAMES has to match the order of the artefacts in final_subset_outlines
              fac = meta_outlines_AR_fac_artefactnames_w_metric_measures_w_dates[order(match(meta_outlines_AR_fac_artefactnames_w_metric_measures_w_dates$ARTEFACTNAME,names(final_subset_outlines$coo))),])

# take a picture of the outlines arranged in a panel
png(filename = file.path("3_output", "final_subset_outlines_panel.png"),
    width = 3000, height = 3000, units = "px", bg = "white")
Momocs::panel(final_subset_outlines,
              cols = "black")
dev.off()

# select a single artefact, check whether the outline plotted here matches the original image of the artefact
Momocs::filter(final_subset_outlines,
               Site == "Nadap") %>% 
  Momocs::panel(names = T)

final_subset_outlines_centered <- Momocs::coo_center(final_subset_outlines)
final_subset_outlines_centered_scaled <- Momocs::coo_scale(final_subset_outlines_centered)

#############################################################################################
#############################################################################################

# harmonic calibration
## Estimates the number of harmonics required for the Fourier methods implemented in Momocs.
final_subset_outlines_centered_scaled_harmonics <- Momocs::calibrate_harmonicpower_efourier(final_subset_outlines_centered_scaled, 
                                                                                         plot = T)  
final_subset_outlines_centered_scaled_harmonics$minh

# Elliptic Fourier Analysis (EFA)
final_subset_outlines_centered_scaled_efourier <- Momocs::efourier(final_subset_outlines_centered_scaled,
                                                                nb.h = as.matrix(final_subset_outlines_centered_scaled_harmonics[["minh"]])[[4,1]], # 4: harmonics for 99.9%
                                                                norm = F,
                                                      start = T) 

# Principal Components Analysis (PCA) on the EFA
final_subset_outlines_centered_scaled_PCA <- Momocs::PCA(final_subset_outlines_centered_scaled_efourier) # PCA on Coe objects, using prcomp.

# save the PCA data of the centered and scaled outlines
saveRDS(final_subset_outlines_centered_scaled_PCA,
        file.path("1_data", 
                  "Outlines",
                  "final_subset_outlines_centered_scaled_seed1_PCA.RDS"))

#############################################################################################
#############################################################################################

## check which PC axis represents what part of the shape space
Momocs::PCcontrib(final_subset_outlines_centered_scaled_PCA,
                  nax = 1:5,
                  sd.r = c(-2.5,-2,-1,0,1,2,2.5))

Momocs::scree_plot(final_subset_outlines_centered_scaled_PCA)

#############################################################################################
#############################################################################################

# combine FAD (first appearance date) with outlines.
# confidence interval? the dates can probably not be exactly the same for each outline due to how the FBDR process works (Rachel Warnock)
# give only the oldest date as terminus _post quem_. we don't know when they went "extinct" which is why we cannot give a LAD (last appearance date) (Cathrine Klein)

taxa_file <- data.frame(taxon = final_subset_outlines_centered_scaled_PCA$fac$ARTEFACTNAME,
                        max = final_subset_outlines_centered_scaled_PCA$fac$median_SPD_age_calBP,
                        min = 0,
                        oneSigma_rangeMax = final_subset_outlines_centered_scaled_PCA$fac$oneSigma_rangeMax,
                        oneSigma_rangeMin = final_subset_outlines_centered_scaled_PCA$fac$oneSigma_rangeMin)
readr::write_tsv(taxa_file,
                 path = file.path("1_data", "final_subset_outlines_centered_scaled_FAD_LAD_C14_oneSigmaMinMax.tsv"))

#############################################################################################

data_w_events <- 
  final_subset_outlines_centered_scaled_PCA$fac %>% 
  dplyr::mutate(Event = dplyr::case_when(median_SPD_age_calBP >= 14600 ~ "GS-2",
                                         median_SPD_age_calBP < 14600 & median_SPD_age_calBP >= 12900~ "GI-1",
                                         median_SPD_age_calBP < 12900 & median_SPD_age_calBP >= 11700~ "GS-1",
                                         median_SPD_age_calBP < 11700 ~ "Holocene")) %>% 
  dplyr::mutate(Event = factor(Event, 
                               levels = c("GS-2", "GI-1", "GS-1", "Holocene"))) %>% 
  unique() %>% 
  dplyr::arrange(desc(median_SPD_age_calBP)) %>% 
  dplyr::mutate(...1 = 1:nrow(.)) %>% 
  dplyr::rename(row_ID = ...1)

readr::write_csv(data_w_events,
                 file = file.path("1_data", "data_w_events.csv"))






# take a picture of the outlines arranged in a panel
## GS-2
png(filename = file.path("3_output", "Figures", "final_subset_outlines_panel_GS-2.png"),
    width = 3000, height = 3000, units = "px", bg = "white")
Momocs::filter(final_subset_outlines,
               ARTEFACTNAME %in% subset(data_w_events, Event == "GS-2")$ARTEFACTNAME) %>% 
  Momocs::panel(
    cols = "black")
dev.off()
## GI-1
png(filename = file.path("3_output", "Figures", "final_subset_outlines_panel_GI-1.png"),
    width = 3000, height = 3000, units = "px", bg = "white")
Momocs::filter(final_subset_outlines,
               ARTEFACTNAME %in% subset(data_w_events, Event == "GI-1")$ARTEFACTNAME) %>% 
  Momocs::panel(
    cols = "black")
dev.off()
## GS-1
png(filename = file.path("3_output", "Figures", "final_subset_outlines_panel_GS-1.png"),
    width = 3000, height = 3000, units = "px", bg = "white")
Momocs::filter(final_subset_outlines,
               ARTEFACTNAME %in% subset(data_w_events, Event == "GS-1")$ARTEFACTNAME) %>% 
  Momocs::panel(
    cols = "black")
dev.off()
## Holocene
png(filename = file.path("3_output", "Figures", "final_subset_outlines_panel_Holocene.png"),
    width = 3000, height = 3000, units = "px", bg = "white")
Momocs::filter(final_subset_outlines,
               ARTEFACTNAME %in% subset(data_w_events, Event == "Holocene")$ARTEFACTNAME) %>% 
  Momocs::panel(
    cols = "black")
dev.off()





