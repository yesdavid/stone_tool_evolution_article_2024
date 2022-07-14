# install.packages("neotoma")

search_query_csv <- readr::read_csv(file = file.path("1_data",
                                                     "neotoma_search_query_results.csv"))

library(neotoma)

dataset <- neotoma::get_dataset(x = search_query_csv$siteid[1:10])
neotoma::get_site(sitename = search_query_csv$sitename[1])

# data <- neotoma::get_download(x = search_query_csv$siteid[1:10])

dataset <- get_dataset(taxonname='Poaceae', 
                       ageyoung = 10000,
                       ageold = 16000,
                       loc = c(-13, 34, 35, 62))
dataset[1:10]
dataset

for (i in names(dataset)) {
  dataset[[i]]$site.data
}

search_query_csv[1:10,]

data <- get_download(dataset)
   
saveRDS(object = data,
        file = file.path("3_output",
                         "neotoma_pollendata.RDS"))

compiled.sites <- compile_taxa(data, list.name='WS64')


#  Extract the Pseudotsuga curves for the sites:
get.curve <- function(x, taxa) {
  if (taxa %in% colnames(x$counts)) {
    count <- x$counts[,taxa]/rowSums(x$counts, na.rm=TRUE)
  } else {
    count <- rep(0, nrow(x$count))
  }
  data.frame(site = x$dataset$site.data$site.name,
             age = x$sample.meta$age,
             count = count)
}

curves <- do.call(rbind.data.frame,
                  lapply(compiled.sites, get.curve, taxa = 'Poaceae'))

data_with_pollen <- 
dplyr::left_join(search_query_csv,
                 curves,
                 by = c("sitename" = "site")) %>% 
  subset(., age > 9000) %>% 
  subset(., age < 16000)



library(rworldmap)
library(ggplot2)

world <- rgeos::gBuffer(rworldmap::getMap(resolution = "high"), byid=TRUE, width=0)
clipper_europe <- as(raster::extent(-13, 34, 35, 62), "SpatialPolygons")
proj4string(clipper_europe) <- CRS(proj4string(world))
world_clip <- raster::intersect(world, clipper_europe)
world_clip_f <- fortify(world_clip)

plot_map <- 
ggplot() +
  geom_polygon(data = world_clip_f, 
               aes(x = long, y = lat, group = group),
               fill = NA, colour = "grey") +   
  # scale_fill_manual(values = nicolas_colors_without_outliers) +
  # facet_wrap(~cluster,
  #            scales = "fixed", 
  #            labeller = as_labeller(cluster_names)) +
  coord_quickmap() +  
  theme_classic() + 
  xlab("Longitude") +
  ylab("Latitude")  +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "none",
        text = element_text(size=20))


  plot_map + 
    geom_jitter(data = data_with_pollen,
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
#############################################################################################
  # Zone A #
  # above 55N
  zone_a <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude >= 55),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  #Zone A1 #
  # above 55N east of 7E
  zone_a1 <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude >= 55 & longitude <= 4),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  #Zone A2 #
  # above 55N east of 7E
  zone_a2 <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude >= 55 & longitude >= 4 & longitude <= 17),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  #Zone A3 #
  # above 55N east of 7E
  zone_a3 <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude >= 55 & longitude >= 17),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
#############################################################################################
  # Zone B #
  #below 55N and above 47N
  zone_b <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude <= 55 & latitude >= 47.5),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  
  # Zone B1 #
  #below 55N and above 47N
  # east of 7E #
  zone_b1 <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude <= 55 & latitude >= 47.5 & longitude <= 4),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  # Zone B2 #
  #below 55N and above 47N
  # west of 7E and east of 17E#
  zone_b2 <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude <= 55 & latitude >= 47.5 & longitude >= 4 & longitude <= 17),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  
  # Zone B3 #
  #below 55N and above 47N
  # east of 7W #
  zone_b3 <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude <= 55 & latitude >= 47.5 & longitude >= 17),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  #############################################################################################
  
  # Zone C #
  # below 47N
  zone_c <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude <= 47.5),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  # Zone C1 #
  # below 47N
  zone_c1 <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude <= 47.5 & longitude <= 4),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  # Zone C2 #
  # below 47N
  zone_c2 <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude <= 47.5 & longitude >= 4 & longitude <= 17),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  # Zone C3 #
  # below 47N
  zone_c3 <- 
  plot_map + 
    geom_jitter(data = subset(data_with_pollen, latitude <= 47.5 & longitude >= 17),
                aes(x = longitude, y = latitude, fill = age), 
                shape = 21,
                size = 3,
                width = 0.1, height = 0.1) 
  
  
  
  #############################################################################################
  
  # pollen curve Poaceae
  ggplot(data = data_with_pollen,
         aes(x = age, y = count)) +
    geom_smooth(span = 0.4) +
    scale_x_reverse()
  
  ##########################
  # above 55N
  pollen_zone_a <- 
  ggplot(data = subset(data_with_pollen, latitude >= 55),
         aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  pollen_zone_a1 <- 
    ggplot(data = subset(data_with_pollen, latitude >= 55 & longitude <= 4 ),
           aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  pollen_zone_a2 <- 
    ggplot(data = subset(data_with_pollen, latitude >= 55 & longitude >= 4 & longitude <= 17) ,
           aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  pollen_zone_a3 <- 
    ggplot(data = subset(data_with_pollen, latitude >= 55 & longitude >= 17),
           aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  
  
  ##########################
  #below 55N and above 47N
  pollen_zone_b <- 
  ggplot(data = subset(data_with_pollen, latitude <= 55 & latitude >= 47.5),
         aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  pollen_zone_b1 <- 
    ggplot(data = subset(data_with_pollen, latitude <= 55 & latitude >= 47.5 & longitude <= 4),
           aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  pollen_zone_b2 <- 
    ggplot(data = subset(data_with_pollen, latitude <= 55 & latitude >= 47.5 & longitude >= 4 & longitude <= 17 ),
           aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  pollen_zone_b3 <- 
    ggplot(data = subset(data_with_pollen, latitude <= 55 & latitude >= 47.5& longitude >= 17),
           aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  
  
  
  ##########################
  #below 47N
  pollen_zone_c <- 
  ggplot(data = subset(data_with_pollen, latitude <= 47.5),
         aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  pollen_zone_c1 <- 
    ggplot(data = subset(data_with_pollen, latitude <= 47.5 & longitude <= 4),
           aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  pollen_zone_c2 <- 
    ggplot(data = subset(data_with_pollen, latitude <= 47.5  & longitude >= 4 & longitude <= 17 ),
           aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  pollen_zone_c3 <- 
    ggplot(data = subset(data_with_pollen, latitude <= 47.5 & longitude >= 17),
           aes(x = age, y = count)) +
     geom_point() + geom_smooth(span = 0.4) +ylim(0,0.6)+ ylim(0,0.6)+
    scale_x_reverse()
  
  
  
  
  #############################################################################################
  combined_plots_thirds <- 
  cowplot::plot_grid(zone_a, zone_b, zone_c,
                     pollen_zone_a +
                       geom_vline(xintercept = 13006,  # laacher see volcano
                                  color = "green", size = 2) + 
                       geom_vline(xintercept = 14600,  # end of late pleniglacial
                                  color = "red", size = 2) +
                       geom_vline(xintercept = 12900,  # end of bølling allerød complex
                                  color = "red", size = 2) + 
                       geom_vline(xintercept = 11700,  # end of younger dryas complex
                                  color = "red", size = 2), 
                     pollen_zone_b +
                       geom_vline(xintercept = 13006,  # laacher see volcano
                                  color = "green", size = 2) + 
                       geom_vline(xintercept = 14600,  # end of late pleniglacial
                                  color = "red", size = 2) +
                       geom_vline(xintercept = 12900,  # end of bølling allerød complex
                                  color = "red", size = 2) + 
                       geom_vline(xintercept = 11700,  # end of younger dryas complex
                                  color = "red", size = 2), 
                     pollen_zone_c +
                       geom_vline(xintercept = 13006,  # laacher see volcano
                                  color = "green", size = 2) + 
                       geom_vline(xintercept = 14600,  # end of late pleniglacial
                                  color = "red", size = 2) +
                       geom_vline(xintercept = 12900,  # end of bølling allerød complex
                                  color = "red", size = 2) + 
                       geom_vline(xintercept = 11700,  # end of younger dryas complex
                                  color = "red", size = 2) ,
                     metrics_zone_a, metrics_zone_b, metrics_zone_c,
                     ncol = 3,
                     labels = c("A", "", "",
                                "B", "", "",
                                "C", "", ""))
  

  
  combined_plots_ninths <- 
    cowplot::plot_grid(zone_a1, zone_a2,zone_a3,
                       pollen_zone_a1, pollen_zone_a2, pollen_zone_a3, 
                       zone_b1, zone_b2, zone_b3, 
                       pollen_zone_b1, pollen_zone_b2, pollen_zone_b3, 
                       zone_c1, zone_c2, zone_c3,
                       pollen_zone_c1,pollen_zone_c2,pollen_zone_c3,
                       ncol = 3)
  
  
  
  
  
  
