library(magrittr)
library(ggplot2)
#############################################################################################
#############################################################################################

data_w_events_raw <- 
  readr::read_csv(file = file.path("1_data", "data_w_events.csv"))

site_names_ID_vs_clean <- 
  readr::read_csv(file = file.path("1_data", "site_names_ID_vs_clean.csv"))

data_w_events_raw <- 
  dplyr::left_join(data_w_events_raw,
                   site_names_ID_vs_clean,
                   by = "Site")

data_w_events_raw$Event <- factor(data_w_events_raw$Event, 
                              levels = c("GS-2", "GI-1", "GS-1", "Holocene"))

data_w_events <- 
  data_w_events_raw %>% 
  dplyr::arrange(., 
                 Event, desc(Lat), Long) %>% 
  dplyr::mutate(row_ID_by_Event_Lat_Long = 1:nrow(.))


final_subset_outlines_centered_scaled_PCA <- readRDS(file = file.path("1_data",
                                                                      "Outlines",
                                                                      "final_subset_outlines_centered_scaled_seed1_PCA.RDS"))
final_subset_outlines_centered_scaled_PCA$fac <- 
  dplyr::left_join(final_subset_outlines_centered_scaled_PCA$fac,
                   dplyr::select(data_w_events, ARTEFACTNAME, Event),
                   by="ARTEFACTNAME")
final_subset_outlines_centered_scaled_PCA$fac <-
  dplyr::left_join(final_subset_outlines_centered_scaled_PCA$fac,
                   site_names_ID_vs_clean,
                   by = "Site")

taxa_file <- readr::read_tsv(file.path("1_data",
                                       "final_subset_outlines_centered_scaled_FAD_LAD_C14_oneSigmaMinMax.tsv"))



#############################################################################################
#############################################################################################


data_w_events %>% 
  ggplot2::ggplot(data = .,
                  ggplot2::aes(y = forcats::fct_reorder(TaxUnit, -median_SPD_age_calBP), 
                               x = median_SPD_age_calBP)) +
  ggplot2::geom_rect(ggplot2::aes(ymin = levels(forcats::fct_reorder(TaxUnit, -median_SPD_age_calBP))[1], 
                                  ymax = levels(forcats::fct_reorder(TaxUnit, -median_SPD_age_calBP))[length(levels(forcats::fct_reorder(TaxUnit, -median_SPD_age_calBP)))], 
                                  xmin = 16000, xmax = 14600), 
                     color="transparent", 
                     fill="grey90", 
                     alpha=0.03) +
  ggplot2::geom_rect(ggplot2::aes(ymin = levels(forcats::fct_reorder(TaxUnit, -median_SPD_age_calBP))[1], 
                                  ymax = levels(forcats::fct_reorder(TaxUnit, -median_SPD_age_calBP))[length(levels(forcats::fct_reorder(TaxUnit, -median_SPD_age_calBP)))], 
                                  xmin = 12900, xmax = 11700), 
                     color="transparent", 
                     fill="grey90", 
                     alpha=0.03) +
  ggplot2::geom_line(linewidth = 3,
                     color = "grey40") +
  ggplot2::geom_point(ggplot2::aes(y = forcats::fct_reorder(TaxUnit, -median_SPD_age_calBP),
                                   x = median_SPD_age_calBP,
                                   fill = Region),
                      color = "white",
                      shape = 21,
                      size = 3) +
  ggplot2::scale_x_reverse() +
  ggplot2::theme_bw()



#############################################################################################
#############################################################################################

pc_contrib_plot_nax1_9 <- 
  Momocs::PCcontrib(final_subset_outlines_centered_scaled_PCA,
                    nax = 1:9,
                    sd.r = c(-2, -1, 0, 1, 2))$gg
ggplot2::ggsave(plot = pc_contrib_plot_nax1_9,
                filename = file.path("3_output",
                                     "pc_contrib_1-9.png"))
ggplot2::ggsave(plot = pc_contrib_plot_nax1_9,
                filename = file.path("3_output",
                                     "pc_contrib_1-9.eps"))


Momocs::scree(final_subset_outlines_centered_scaled_PCA)

pc_scree_plot <- 
  Momocs::scree_plot(final_subset_outlines_centered_scaled_PCA) +
  ggplot2::theme_bw()
ggplot2::ggsave(plot = pc_scree_plot,
                filename = file.path("3_output",
                                     "screeplot.png"))
ggplot2::ggsave(plot = pc_scree_plot,
                filename = file.path("3_output",
                                     "screeplot.eps"))


cowplot_PCA_scree_contrib_plot <- 
  cowplot::plot_grid(pc_scree_plot,
                     pc_contrib_plot_nax1_9,
                     labels = "AUTO",
                     ncol = 2,
                     rel_widths = c(2, 1))
cowplot_PCA_scree_contrib_plot

ggplot2::ggsave(plot = cowplot_PCA_scree_contrib_plot,
                filename = file.path("3_output",
                                     "screeplot+pccontrib_skizze.png"),
                width = 20, height = 20*2/3, units = "cm",
                bg = "white")
ggplot2::ggsave(plot = cowplot_PCA_scree_contrib_plot,
                filename = file.path("3_output",
                                     "screeplot+pccontrib_skizze.eps"),
                width = 20, height = 20*2/3, units = "cm",
                bg = "white")

#############################################################################################
# PCA
#############################################################################################

pca_data <- as.data.frame(final_subset_outlines_centered_scaled_PCA$x)
pca_data$ARTEFACTNAME <- rownames(pca_data)

pca_data_metaInfo <-
  dplyr::left_join(final_subset_outlines_centered_scaled_PCA$fac,
                   pca_data,
                   by = "ARTEFACTNAME")

a <- 
  ggplot2::ggplot(data = pca_data_metaInfo,
                  ggplot2::aes(x = PC1, y = PC2,
                               fill = as.factor(Event))) +
  geom_point(size = 3,
             # fill = "white",
             shape = 21) +
  coord_fixed(ratio =1) +
  theme_classic() +
  theme(legend.position = "right",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  geom_hline(yintercept=0, linetype="dashed", alpha = 0.5) + 
  geom_vline(xintercept=0, linetype="dashed", alpha = 0.5) +
  ggthemes::scale_fill_colorblind() +
  xlab(paste0("PC1 (", 
              round(final_subset_outlines_centered_scaled_PCA$eig[1]*100, digits = 1), 
              "%)")) +
  ylab(paste0("PC2 (", 
              round(final_subset_outlines_centered_scaled_PCA$eig[2]*100, digits = 1), 
              "%)")) +
  guides(color = FALSE, 
         fill = FALSE) #+
  # geom_polygon(stat = "ellipse", aes(color = as.factor(Event), fill = "white"), alpha = 0)
a

cowplot::plot_grid(a,
                   cowplot::plot_grid(pc_scree_plot,
                                      pc_contrib_plot_nax1_9,
                                      labels = c('B', 'C'),
                                      ncol = 2, byrow = T,
                                      rel_widths = c(2, 1)),
                   labels = c('A', ''),
                   nrow = 2)

abcd <- 
cowplot::plot_grid(a+ theme(aspect.ratio = 1),
                   a+ theme(aspect.ratio = 1)+aes(x=PC3, y=PC4) +
                     xlab(paste0("PC3 (", 
                                 round(final_subset_outlines_centered_scaled_PCA$eig[3]*100, digits = 1), 
                                 "%)")) +
                     ylab(paste0("PC4 (", 
                                 round(final_subset_outlines_centered_scaled_PCA$eig[4]*100, digits = 1), 
                                 "%)")),
                   a+ theme(aspect.ratio = 1)+aes(x=PC5, y=PC6) +
                     xlab(paste0("PC5 (", 
                                 round(final_subset_outlines_centered_scaled_PCA$eig[5]*100, digits = 1), 
                                 "%)")) +
                     ylab(paste0("PC6 (", 
                                 round(final_subset_outlines_centered_scaled_PCA$eig[6]*100, digits = 1), 
                                 "%)")),
                   a+ theme(aspect.ratio = 1)+aes(x=PC7, y=PC8) +
                     xlab(paste0("PC7 (", 
                                 round(final_subset_outlines_centered_scaled_PCA$eig[7]*100, digits = 1), 
                                 "%)")) +
                     ylab(paste0("PC8 (", 
                                 round(final_subset_outlines_centered_scaled_PCA$eig[8]*100, digits = 1), 
                                 "%)")),
                   ncol = 2,
                   labels = "AUTO")
abcd

cowplot::plot_grid(abcd,
                   cowplot::plot_grid(pc_scree_plot,
                                      pc_contrib_plot_nax1_9,
                                      labels = c('E', 'F'),
                                      ncol = 2, byrow = T#,
                                      # rel_widths = c(2, 1)
                                      ),
                   labels = c(''),
                   nrow = 1)
#############################################################################################
#############################################################################################

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
                  aes(y = reorder(Site_nice, -median_SPD_age_calBP), 
                      x=median_SPD_age_calBP,
                      xmin = oneSigma_rangeMax, 
                      xmax = oneSigma_rangeMin,
                      color = -median_SPD_age_calBP)) + 
  ggplot2::geom_pointrange() +
  ggplot2::annotate(xmin = c(16000, 12900), 
                    xmax = c(14600, 11700),
                    ymin = final_subset_outlines_centered_scaled_PCA$fac$Site_nice[1], 
                    ymax = final_subset_outlines_centered_scaled_PCA$fac$Site_nice[87],
                    "rect", 
                    alpha = 0.2, fill = c("grey60", "grey60")) +
  ggplot2::annotate("text", 
                    y = "Unken", 
                    x = c(15300, 13750, 12300, 11100), 
                    label = unique(data_w_events$Event),
                    size = 6) +
  # ggplot2::geom_vline(xintercept = c(14600, 12900, 11700),
  #                     alpha = 0.5,
  #                     linetype = "dotted") +
  scale_x_reverse(limits = c(16000, 
                             round(min(final_subset_outlines_centered_scaled_PCA$fac$oneSigma_rangeMin), -2)), 
                  expand = c(0, 0),
                  minor_breaks = seq(to = 16000,
                                     from = round(min(final_subset_outlines_centered_scaled_PCA$fac$oneSigma_rangeMin), -2),
                                     by = 100),
                  breaks = seq(to = 16000,
                               from = round(min(final_subset_outlines_centered_scaled_PCA$fac$oneSigma_rangeMin), -3),
                               by = 500)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_line(colour="grey95"),
        panel.grid.major = element_line(colour="grey88")) +
  xlab("Age calBP") +
  ylab("Sites") 


plot_of_selected_artefacts_and_ages_ageUncertainties_SITES

ggsave(plot_of_selected_artefacts_and_ages_ageUncertainties_SITES,
       filename = file.path("3_output", "plot_of_selected_artefacts_and_ages_ageUncertainties_SITES.png"),
       width = 23, height = 25, units = "cm", device = "png")
ggsave(plot_of_selected_artefacts_and_ages_ageUncertainties_SITES,
       filename = file.path("3_output", "plot_of_selected_artefacts_and_ages_ageUncertainties_SITES.tiff"),
       width = 23, height = 25, units = "cm", device = "tiff")

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



map_of_selected_artefacts_and_ages_ClimateEvents <-
  base_map +
  ggrepel::geom_text_repel(data = data_w_events,
                           aes(x = Long, y = Lat,
                               # label = paste0(Site_nice, "\n(",median_SPD_age_calBP," calBP)")),
                               label = row_ID_by_Event_Lat_Long), # Site_nice
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
ggsave(map_of_selected_artefacts_and_ages_ClimateEvents,
       filename = file.path("3_output", "map_of_selected_artefacts_and_ages_ClimateEvents.eps"),
       width = 30, height = 30, units = "cm", device = "eps")

cat("", file = file.path("3_output", "data_w_events_Sitename+Number.txt"))
for(i in 1:nrow(data_w_events)){
  cat(paste0(data_w_events[i,"row_ID_by_Event_Lat_Long"], ": ", data_w_events[i,"Site_nice"], ", "),
      append = T,
      file = file.path("3_output", "data_w_events_Sitename+Number.txt"))
  
}

