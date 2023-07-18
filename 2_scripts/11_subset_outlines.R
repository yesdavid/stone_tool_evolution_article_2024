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

## check which PC axis represents what part of the shapespace
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
#############################################################################################
library(ggplot2)

plot_of_selected_artefacts_and_ages_ageUncertainties <- 
  ggplot2::ggplot(data = taxa_file, 
                  aes(y = reorder(taxon, -max), 
                      x=max,
                      xmin = oneSigma_rangeMax, 
                      xmax = oneSigma_rangeMin,
                      color = -max)) + 
  ggplot2::geom_pointrange() +
  scale_x_reverse() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Age calBP") +
  ylab("Artefactnames")

plot_of_selected_artefacts_and_ages_ageUncertainties

ggsave(plot_of_selected_artefacts_and_ages_ageUncertainties,
       filename = file.path("3_output", "plot_of_selected_artefacts_and_ages_ageUncertainties.png"),
       width = 40, height = 25, units = "cm", device = "png")


plot_of_selected_artefacts_and_ages_ageUncertainties_SITES <- 
  ggplot2::ggplot(data = final_subset_outlines_centered_scaled_PCA$fac, 
                  aes(y = reorder(Site, -median_SPD_age_calBP), 
                      x=median_SPD_age_calBP,
                      xmin = oneSigma_rangeMax, 
                      xmax = oneSigma_rangeMin,
                      color = -median_SPD_age_calBP)) + 
  ggplot2::geom_pointrange() +
  scale_x_reverse() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Age calBP") +
  ylab("Sites")

plot_of_selected_artefacts_and_ages_ageUncertainties_SITES

ggsave(plot_of_selected_artefacts_and_ages_ageUncertainties_SITES,
       filename = file.path("3_output", "plot_of_selected_artefacts_and_ages_ageUncertainties_SITES.png"),
       width = 20, height = 25, units = "cm", device = "png")

####################
### distribution map
####################

world <- rgeos::gBuffer(rworldmap::getMap(resolution = "high"), byid=TRUE, width=0)

# create an extent which spans only the distribution of our samples
# potentially, the extents have to be manually in-/decreased
data_extent <- 
  as(raster::extent(min(final_subset_outlines_centered_scaled_PCA$fac$Long, #minimum longitude
                        na.rm = T)-3, 
                    max(final_subset_outlines_centered_scaled_PCA$fac$Long, #maximum longitude
                        na.rm = T)+3, 
                    min(final_subset_outlines_centered_scaled_PCA$fac$Lat, #minimum latitude
                        na.rm = T)-3, 
                    max(final_subset_outlines_centered_scaled_PCA$fac$Lat, #maximum latidude
                        na.rm = T)+3), # order: xmin, xmax, ymin, ymax
     "SpatialPolygons")

sp::proj4string(data_extent) <- sp::CRS(sp::proj4string(world)) # set the coordinate reference system of the data to be the same as the world map.

world_clip <- raster::intersect(world, data_extent) # select only those parts of the world map within our bounding box/extent 

world_clip_f <- fortify(world_clip) # transforms it into a data frame

# base map
base_map <- 
  ggplot() +
  geom_polygon(data = world_clip_f, 
               aes(x = long, 
                   y = lat, 
                   group = group),
               fill = "grey", 
               colour = "grey") +
  coord_fixed() +
  coord_quickmap() +  
  theme_classic() + 
  xlab("Longitude") +
  ylab("Latitude")  +
  # labs(color = "Country") + # capitalize 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  scale_color_gradient(low = "green", high = "brown") +
  theme(legend.position = "right", # or "none"
        text = element_text(size=20))

data_w_events <- 
  final_subset_outlines_centered_scaled_PCA$fac %>% 
  dplyr::mutate(Event = dplyr::case_when(median_SPD_age_calBP >= 14600 ~ "GS-2",
                                         median_SPD_age_calBP < 14600 & median_SPD_age_calBP >= 12900~ "GI-1",
                                         median_SPD_age_calBP < 12900 & median_SPD_age_calBP >= 11700~ "GS-1",
                                         median_SPD_age_calBP < 11700 ~ "Holocene")) %>% 
  dplyr::mutate(Event = factor(Event, 
                               levels = c("GS-2", "GI-1", "GS-1", "Holocene"))) %>% 
  unique() %>% 
  arrange(desc(median_SPD_age_calBP)) %>% 
  dplyr::mutate(...1 = 1:nrow(.)) %>% 
  rename(row_ID = ...1)

map_of_selected_artefacts_and_ages_ClimateEvents <-
  base_map +
  ggrepel::geom_text_repel(data = data_w_events,
                           aes(x = Long, y = Lat,
                               # label = paste0(Site, "\n(",median_SPD_age_calBP," calBP)")),
                               label = row_ID), # Site
                           # alpha = 0.7,
                           size = 4,
                           force = 80,
                           force_pull = 12#,
                           # min.segment.length = 1.5
                           ) +
  geom_jitter(data = data_w_events,
              aes(x = Long, y = Lat,
                  fill = -median_SPD_age_calBP),
              # alpha = 0.7,
              shape = 21,
              size = 3#,
              #width = 0.35, height = 0.5
              ) +
  facet_wrap(~Event) +
  theme(legend.position = "none")

map_of_selected_artefacts_and_ages_ClimateEvents

ggsave(map_of_selected_artefacts_and_ages_ClimateEvents,
       filename = file.path("3_output", "map_of_selected_artefacts_and_ages_ClimateEvents.png"),
       width = 30, height = 30, units = "cm", device = "png")

readr::write_csv(data_w_events,
                 file = file.path("3_output", "data_w_events.csv"))

for(i in 1:nrow(data_w_events)){
  cat(paste0(data_w_events[i,"row_ID"], ": ", data_w_events[i,"Site"], ", "),
      append = T,
      file = file.path("3_output", "data_w_events.txt"))
  
}




