library(outlineR)
library(magrittr)

#############################################################################################
#############################################################################################


# single_outlines_list <- get_outlines(outpath = file.path("1_data",
#                                                          "Outlines",
#                                                          "armatures_to_keep"), 
#                                      tps_file_rescale = NULL)
# 
# outlines_combined <- combine_outlines(single_outlines_list = single_outlines_list)
# 
# 
# length(outlines_combined) #how many outlines do you have?
# # stack(outlines_combined) # shows all outlines above one another(you might want to center and scale them first using Momocs)
# # Momocs::panel(outlines_combined) # shows all outlines next to each other
# # Momocs::inspect(outlines_combined) # shows only a single outline at a time. 
# 
# outlines_AR_fac_w_metric_measures <- 
#   readr::read_csv(file.path("..",
#                             "CLIOARCH_workshop_2_2020","final","1_data_paper","1511NAC_Database","outlines",
#                             "outlines_AR_fac_artefactnames_w_metric_measures.csv"))
# 
# fac <- 
#   subset(outlines_AR_fac_w_metric_measures, 
#          ARTEFACTNAME %in% names(outlines_combined)) %>% 
#   as.data.frame()
# rownames(fac) <- fac$ARTEFACTNAME
# 
# 
# 
# outlines <- 
# Momocs::Out(x = outlines_combined[[1]],
#             fac = fac[match(names(outlines_combined[[1]]), 
#                             fac$ARTEFACTNAME), ]  )
# 
# 
# # only select outlines _with_ metric measurements
# outlines <- Momocs::filter(outlines,
#                            !is.na(SCALE))
# 
# 
# ## centre, scale, define first coordinate
# outlines_centered <- Momocs::coo_centre(outlines) # Returns a shape centered on the origin.
# outlines_centered_scaled <- Momocs::coo_scale(outlines_centered) # Scales the coordinates by the centroid size.
# outlines_centered_scaled <- Momocs::coo_slidedirection(outlines_centered_scaled, 
#                                                        direction = "up") # Sets the "first" coordinate to be on the top.


###################
# load outlines
###################
outlines_AR_centered_w_metric_measures <- readRDS(file = file.path("1_data",
                                                                   "1511NAC_Database_mod",
                                                                   "outlines",
                                                                   "outlines_AR_centered_w_metric_measures.RDS"))

outlines_centered_scaled <- outlines_AR_centered_w_metric_measures #%>% 
  # Momocs::filter(., !is.na(SCALE))



#############################################################################################
# combine with C14 dates, subset to outlines with C14 data available
#############################################################################################

filtered_C14_data_calibrated_for_outlines <- readr::read_csv(file.path("3_output", 
                                                                     "C14",
                                                                     "filtered_C14_data_calibrated_for_outlines.csv"))

# C14_metadata <-
#   dplyr::left_join(outlines_centered_scaled$fac, 
#                    filtered_C14_data_calibrated_for_outlines, 
#                    by = c("Site" = "Site_ID_split", 
#                           "Long", "Lat")) 

C14_metadata <-
  dplyr::left_join(
    filtered_C14_data_calibrated_for_outlines, 
    outlines_centered_scaled$fac,
    by = c("Site_ID_split" = "Site", 
           "Long", "Lat"))

unique(filtered_C14_data_calibrated_for_outlines$Site_ID_split)[!(unique(filtered_C14_data_calibrated_for_outlines$Site_ID_split) %in% 
                                                                  unique(outlines_centered_scaled$fac$Site))]
unique(filtered_C14_data_calibrated_for_outlines$Site_ID_split)[(unique(filtered_C14_data_calibrated_for_outlines$Site_ID_split) %in% 
                                                                  unique(outlines_centered_scaled$fac$Site))]


table(C14_metadata$Region)


outlines_centered_scaled_C14 <- 
  Momocs::Out(x = outlines_centered_scaled[[1]],
              fac = C14_metadata[match(names(outlines_centered_scaled[[1]]), 
                                       C14_metadata$ARTEFACTNAME), ]  ) %>% 
  Momocs::filter(., !is.na(median_SPD_age_calBP)) 

saveRDS(outlines_centered_scaled_C14,
        file.path("1_data",
                  "Outlines",
                  "outlines_toKeep_centered_scaled_C14metadata.RDS"))

# #############################################################################################
# #############################################################################################
# 
# subset to control
# 
# metadata <- 
#   dplyr::left_join(outlines_centered_scaled_C14$fac,
#                    data_w_events)
# 
# split_FUN <- function(x){
#   sapply(x, strsplit(x, split = "/")[[1]][4])
# }
# 
# all_AR_images_single <-
#   data.frame(short_path = 
#                list.files("/home/au656892/Documents/Doktor/2_projects/CLIOARCH_workshop_2_2020/final/1_data_paper/1_data/artefact_images_by_region/", 
#                           recursive = T, 
#                           pattern = "_AR_"))  %>% 
#   dplyr::mutate(., ARTEFACTNAME = sapply(strsplit(.$short_path, split = "/"), function(x) x[[4]])) %>% 
#   dplyr::mutate(., ARTEFACTNAME = sapply(strsplit(.$ARTEFACTNAME, split = ".jpg"), function(x) x[[1]])) %>% 
#   dplyr::mutate(., long_path = paste0("/home/au656892/Documents/Doktor/2_projects/CLIOARCH_workshop_2_2020/final/1_data_paper/1_data/artefact_images_by_region/",
#                                       .$short_path)) %>% 
#   subset(.,
#          ARTEFACTNAME %in% metadata$ARTEFACTNAME) %>% 
#   select(-short_path)
# 
# head(all_AR_images_single)
# 
# for(i_EVENT in unique(metadata$Event)){
#   
#   current_event_subset <- 
#     subset(metadata, Event == i_EVENT)
#   
#   for(i_REGION in unique(current_event_subset$Macro_region_code)){
#     
#     new_path <- file.path("1_data", "Outlines", "armatures_to_keep_new", i_EVENT, i_REGION)
#     dir.create(new_path, recursive = T)
#     
#     current_region_current_event_subset <- 
#       subset(current_event_subset, Macro_region_code == i_REGION)$ARTEFACTNAME
#     
#     old_long_path <- all_AR_images_single[which(all_AR_images_single$ARTEFACTNAME %in% current_region_current_event_subset), "long_path"]
#     
#     file.copy(from = old_long_path,
#               to = new_path)
#     
#   }
#   
# }
# 
# 
# #############################################################################################
# #############################################################################################

# set.seed(1)
# subsample_stratified_n1_perSite <-
# splitstackshape::stratified(outlines_centered_scaled_C14$fac,
#                             group = "Site", # column "site" is site+layer; however: some of the unique site/layer combinations have the same 14C dates associated with them.
#                             size = 1)

cleaned_subset <-
  Momocs::filter(outlines_centered_scaled_C14,
                 ARTEFACTNAME %in% gsub(x = list.files(path = file.path("1_data",
                                                                        "Outlines",
                                                                        "armatures_to_keep_new_fr",
                                                                        "all")),
                                        pattern = ".jpg",
                                        replacement = ""))
  

set.seed(1)
subsample_stratified_n1_perSite <-
splitstackshape::stratified(cleaned_subset$fac, #outlines_centered_scaled_C14$fac,
                            group = c("TaxUnit", "Region",
                                      # "TaxUnit_unique", 
                                      "median_SPD_age_calBP"
                                      ), # column "site" is site+layer
                            size = 1)

outlines_centered_scaled_subs <-
  Momocs::filter(outlines_centered_scaled_C14,
                 ARTEFACTNAME %in% subsample_stratified_n1_perSite$ARTEFACTNAME)

outlines_centered_scaled_subs$fac$Region <- factor(outlines_centered_scaled_subs$fac$Region)
# outlines_centered_scaled_subs$fac$Expert_editor <- factor(outlines_centered_scaled_subs$fac$Expert_editor)

Momocs::panel(outlines_centered_scaled_subs,
              fac = "Region")

table(outlines_centered_scaled_subs$fac$Region)

filtered_C14_data_calibrated_for_outlines$Site_ID_split[!(filtered_C14_data_calibrated_for_outlines$Site_ID_split %in% outlines_centered_scaled_subs$fac$Site_ID_split)]




for(i in outlines_centered_scaled_subs$fac$ARTEFACTNAME){
  
  file.copy(from = paste0("/home/au656892/Documents/Doktor/2_projects/stone_tool_evolution_article_2022/1_data/Outlines/armatures_to_keep_new_fr/all/",
                          i,
                          ".jpg"),
            to = paste0("/home/au656892/Documents/Doktor/2_projects/stone_tool_evolution_article_2022/1_data/Outlines/final_subset_outlines/",
                        i,
                        ".jpg"))
}

final_subset_outlines <- 
  outlineR::get_outlines(outpath = file.path("1_data",
                                             "Outlines",
                                             "final_subset_outlines")) %>% 
  outlineR::combine_outlines()

final_subset_outlines <-
  Momocs::Out(final_subset_outlines$coo,
              fac = outlines_centered_scaled_subs$fac[order(match(outlines_centered_scaled_subs$fac$ARTEFACTNAME,names(final_subset_outlines$coo))),])

png(filename = file.path("3_output", "final_subset_outlines_panel.png"),
    width = 3000, height = 3000, units = "px", bg = "white")
Momocs::panel(final_subset_outlines,
              cols = "black")
dev.off()

Momocs::filter(final_subset_outlines,
               Site_ID_split == "Nadap") %>% 
  Momocs::panel()

final_subset_outlines_centered <- Momocs::coo_center(final_subset_outlines)
final_subset_outlines_centered_scaled <- Momocs::coo_scale(final_subset_outlines_centered)



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

saveRDS(final_subset_outlines_centered_scaled_PCA,
        file.path("1_data", 
                  "Outlines",
                  paste0("TAXA", nrow(final_subset_outlines_centered_scaled_PCA$x), 
                         "_PCs", ncol(final_subset_outlines_centered_scaled_PCA$x), 
                         "_final_subset_outlines_centered_scaled_seed1_PCA.RDS")))




hclust(d = dist(final_subset_outlines_centered_scaled_PCA$x),
       method = "ward.D2") %>% 
  plot()


#############################################################################################
#############################################################################################

## check which PC axis represents what part of the shapespace

Momocs::PCcontrib(final_subset_outlines_centered_scaled_PCA,
                  nax = 1:5,
                  sd.r = c(-2.5,-2,-1,0,1,2,2.5))
# 
Momocs::scree_plot(final_subset_outlines_centered_scaled_PCA)
# 
Momocs::plot_PCA(final_subset_outlines_centered_scaled_PCA,
                 f = "Region",
                 # axes = c(2,3),
                 morphospace = T,
                 labelpoints = F)

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


plot_of_selected_artefacts_and_ages_ageUncertainties <- 
  ggplot2::ggplot(data = taxa_file, 
                  aes(y = reorder(taxon, -max), 
                      x=max,
                      xmin = oneSigma_rangeMax, 
                      xmax = oneSigma_rangeMin,
                      color = max)) + 
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
                  aes(y = reorder(Site_ID_split, -median_SPD_age_calBP), 
                      x=median_SPD_age_calBP,
                      xmin = oneSigma_rangeMax, 
                      xmax = oneSigma_rangeMin,
                      color = median_SPD_age_calBP)) + 
  ggplot2::geom_pointrange() +
  scale_x_reverse() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Age calBP") +
  ylab("Sites")

ggsave(plot_of_selected_artefacts_and_ages_ageUncertainties_SITES,
       filename = file.path("3_output", "plot_of_selected_artefacts_and_ages_ageUncertainties_SITES.png"),
       width = 40, height = 25, units = "cm", device = "png")




####################
### distribution map
####################

world <- rgeos::gBuffer(rworldmap::getMap(resolution = "high"), byid=TRUE, width=0)

# create an extent which spans only the distribution of our samples
# potentially, the extents have to be manually in-/decreased
data_extent <- as(raster::extent(min(filtered_C14_data_calibrated$Long, #minimum longitude
                                     na.rm = T)-3, 
                                 max(filtered_C14_data_calibrated$Long, #maximum longitude
                                     na.rm = T)+3, 
                                 min(filtered_C14_data_calibrated$Lat, #minimum latitude
                                     na.rm = T)-3, 
                                 max(filtered_C14_data_calibrated$Lat, #maximum latidude
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
  mutate(Event = case_when(median_SPD_age_calBP >= 14600 ~ "GS-2",
                           median_SPD_age_calBP < 14600 & median_SPD_age_calBP >= 12900~ "GI-1",
                           median_SPD_age_calBP < 12900 & median_SPD_age_calBP >= 11700~ "GS-1",
                           median_SPD_age_calBP < 11700 ~ "Holocene")) %>% 
  mutate(Event = factor(Event, 
                        levels = c("GS-2", "GI-1", "GS-1", "Holocene"))) %>% 
  # select(., Event, Long, Lat, median_SPD_age_calBP, Region, 
  #        # TaxUnit_unique, 
  #        TaxUnit) %>% 
  unique()

map_of_selected_artefacts_and_ages_ClimateEvents <- 
  base_map +
  ggrepel::geom_label_repel(data = data_w_events,
                           aes(x = Long, y = Lat,
                               label = paste0(Site_ID_split, "\n(",median_SPD_age_calBP," calBP)")),
                           # alpha = 0.7,
                           size = 4,
                           min.segment.length = 0) +
  geom_jitter(data = data_w_events,
              aes(x = Long, y = Lat,
                  fill = -median_SPD_age_calBP),
              # alpha = 0.7,
              shape = 21,
              size = 4,
              nudge_x = 0.6,
              nudge_y = 0.3) +
  facet_wrap(~Event) +
  theme(legend.position = "none")


ggsave(map_of_selected_artefacts_and_ages_ClimateEvents,
       filename = file.path("3_output", "map_of_selected_artefacts_and_ages_ClimateEvents.png"),
       width = 50, height = 40, units = "cm", device = "png")
