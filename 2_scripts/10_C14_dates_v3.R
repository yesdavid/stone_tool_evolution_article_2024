
library(dplyr)
library(magrittr)
library(ggplot2)
library(rcarbon)



####################
# load the site-data
####################
datasheet_sites <- readr::read_csv(file = file.path("1_data",
                                                    "1511NAC_Database_mod",
                                                    "MASTERTABLE_v8_all_revised_sites_with_outline_siteID.csv"))

# remove sites with problematic dating quality
rank_1_sites <- 
datasheet_sites %>% 
  dplyr::filter(Dating_qual != "Problematic")  %>% 
  dplyr::mutate(Site_fau2 = dplyr::case_when(
    Site_fau %in% "Preseved" ~ TRUE, 
    Site_fau %in% "No" ~ FALSE)) %>% 
  # dplyr::filter(Site_strat == "Stratified") %>%  
  # dplyr::filter(Site_fau2 == TRUE) %>% # if TRUE: faunal remains MUST be present in the site
  dplyr::select("Site_ID", "Key_site","Level_layer_concentration",
                "TaxUnit","TaxUnit_unique",
                "Timeslice","Long","Lat",
                "Dating_BP", "Dating_method","Dating_qual",
                "Site_type","Site_strat","Ass_pos","Ass_coh",
                "Site_fau2","Site_excav","Key_references",
                "Macro_region","Macro_region_code")

# table(rank_1_sites$Macro_region)
# table(rank_1_sites$Site_strat)
# table(rank_1_sites$Site_type)
# table(rank_1_sites$Dating_qual)
# table(rank_1_sites$Dating_method)
# 
# head(rank_1_sites)


# unwrap the Site_ID in the site-data
rank_1_sites$rank_1_sites_unique_row_ID <- 1:nrow(rank_1_sites)

rank_1_sites_list_separate_siteID <- list()
for(i in 1:nrow(rank_1_sites)){
  
  current_row <- rank_1_sites[i,]
  
  split_Site_ID <- strsplit(current_row$Site_ID, ",")[[1]]
  
  if(length(split_Site_ID) > 1){
    
    current_Site_ID_list <- list()
    for (current_Site_ID in 1:length(split_Site_ID)){
      
      current_Site_ID_list[[current_Site_ID]] <- 
        data.frame(Site_ID_split = split_Site_ID[[current_Site_ID]],
                   Site_ID = current_row$Site_ID,
                   rank_1_sites_unique_row_ID = current_row$rank_1_sites_unique_row_ID)
    }
    df <- do.call(rbind.data.frame, current_Site_ID_list)
    
  }  else if(is.na(current_row$Site_ID)) {
    
    df <- 
      data.frame(Site_ID_split = current_row$Key_site,
                 Site_ID = NA,
                 rank_1_sites_unique_row_ID = current_row$rank_1_sites_unique_row_ID)
    
  } else {
    
    df <- 
      data.frame(Site_ID_split = strsplit(current_row$Site_ID, ",")[[1]],
                 Site_ID = current_row$Site_ID,
                 rank_1_sites_unique_row_ID = current_row$rank_1_sites_unique_row_ID)
    
  }
  
  
  rank_1_sites_list_separate_siteID[[i]] <- df
}
rank_1_sites_list_separate_siteID_df <- 
  dplyr::left_join(do.call(rbind.data.frame, rank_1_sites_list_separate_siteID),
                   rank_1_sites,
                   by = "rank_1_sites_unique_row_ID")


###################
# load outlines
###################
outlines_AR_centered_w_metric_measures <- readRDS(file = file.path("1_data",
                                                                   "1511NAC_Database_mod",
                                                                   "outlines",
                                                                   "outlines_AR_centered_w_metric_measures.RDS"))

outlines_AR_centered_w_metric_measures


#####################
# load 14C data sheet
#####################
C14_data <- 
  readr::read_csv(file = file.path("1_data",
                                   "CLIOARCH_14C_datasheet_v5.csv")) %>% 
  unique()
# View(C14_data)
# names(C14_data)

# clean 14C spreadsheet
C14_data_selected <- 
  C14_data %>% 
  dplyr::select(!c("...19")) %>% 
  dplyr::filter(!(`LabCode/Ref` %in% "?") & !is.na(`LabCode/Ref`)) %>% # remove entries where 14C-labcode is unknown/not available
  unique()


C14_data_selected$C14_data_selected_unique_row_ID <- 1:nrow(C14_data_selected)

C14_data_selected_list_separate_siteID <- list()
for(i in 1:nrow(C14_data_selected)){
  
  current_row <- C14_data_selected[i,]
  
  split_Site_ID <- strsplit(current_row$`Site_ID (from main datasheet)`, ",")[[1]]
  
  if(length(split_Site_ID) > 1){
    
    current_Site_ID_list <- list()
    for (current_Site_ID in 1:length(split_Site_ID)){
      
      current_Site_ID_list[[current_Site_ID]] <- 
        data.frame(Site_ID_split = split_Site_ID[[current_Site_ID]],
                   Site_ID = current_row$`Site_ID (from main datasheet)`,
                   C14_data_selected_unique_row_ID = current_row$C14_data_selected_unique_row_ID)
    }
    df <- do.call(rbind.data.frame, current_Site_ID_list)
    
  }  else if(is.na(current_row$`Site_ID (from main datasheet)`)) {
    
    df <- 
      data.frame(Site_ID_split = current_row$`Site (from main datasheet + 14c-palaeolithic)`,
                 Site_ID = NA,
                 C14_data_selected_unique_row_ID = current_row$C14_data_selected_unique_row_ID)
    
  } else {
    
    df <- 
      data.frame(Site_ID_split = strsplit(current_row$`Site_ID (from main datasheet)`, ",")[[1]],
                 Site_ID = current_row$`Site_ID (from main datasheet)`,
                 C14_data_selected_unique_row_ID = current_row$C14_data_selected_unique_row_ID)
    
  }
  
  
  C14_data_selected_list_separate_siteID[[i]] <- df
}
C14_data_selected_list_separate_siteID_df <- 
  dplyr::left_join(do.call(rbind.data.frame, C14_data_selected_list_separate_siteID),
                   C14_data_selected,
                   by = "C14_data_selected_unique_row_ID")



###################################
# join quality sites with 14C data
###################################
rank_1_sites_C14 <-
  dplyr::left_join(rank_1_sites_list_separate_siteID_df %>% 
                     select(-"rank_1_sites_unique_row_ID"),
                   C14_data_selected_list_separate_siteID_df %>% 
                     select(-"C14_data_selected_unique_row_ID"),
                   relationship = "many-to-many",
                   by = c("Site_ID_split"#,
                          #"Site_ID" = "Site_ID (from main datasheet)",
                          # "Key_site" = "Site (from main datasheet)"#,
                          #"Level_layer_concentration" = "Layer/Context (from main datasheet)"#,
                          #"TaxUnit" = "TaxUnit (from main datasheet)"
                          )) %>% 
  dplyr::filter((!is.na(Conv_C14_date) | !is.na(AMS_C14_date) | !is.na(`OSL/TL date`) )) %>% # select only the ones that have a dating
  dplyr::select(-c("Site_ID.x", "Site_ID.y")) %>% 
  unique() %>%
  dplyr::as_tibble()


#######################
list_separate_siteID_df_separate_siteId <- rank_1_sites_C14


# which sites in the meta data match the sites in the outline?
unique_sites_artefacts <-
Momocs::filter(outlines_AR_centered_w_metric_measures,
               Site %in% list_separate_siteID_df_separate_siteId$Site_ID_split)$fac$Site %>% 
  unique()


# meta data subset by the available outlines
filtered_C14_data <- 
  list_separate_siteID_df_separate_siteId %>% 
  filter(Site_ID_split %in% unique_sites_artefacts) 

####################
# calibrate 14C-data
####################
filtered_C14_data$median_SPD_age_calBP <- NA
filtered_C14_data$oneSigma_rangeMin <- NA
filtered_C14_data$oneSigma_rangeMax <- NA
for(i in unique(filtered_C14_data$Site_ID_split)){

    a <-
      dplyr::filter(filtered_C14_data, Site_ID_split == i) 
    
    a_ams <-
    a %>% 
      dplyr::select("AMS_C14_date",
                    `AMS_Pm/Std`) %>% 
      na.omit()
    
    a_conv <-
    a %>% 
      dplyr::select("Conv_C14_date",
                    "Conv_std") %>% 
      na.omit()
    
    caldates <-
      rcarbon::calibrate(x=na.omit(c(a_ams$AMS_C14_date, a_conv$Conv_C14_date)),
                         errors=na.omit(c(a_ams$`AMS_Pm/Std`, a_conv$Conv_std)),
                         calCurves='intcal20')
    
    if(nrow(a) > 1){
      spd_out <- rcarbon::spd(caldates,
                              timeRange=c(19000,9000))
    } else {
      spd_out <- caldates
      spd_out$grid <- spd_out$grids$`1`
    }
    
    filtered_C14_data[which(filtered_C14_data$Site_ID_split == i),"median_SPD_age_calBP"] <- 
      spd_out$grid %>% 
      arrange(., desc(calBP)) %>%
      mutate(cumsum_prdens = cumsum(PrDens)) %>% 
      filter(cumsum_prdens <= 0.5*sum(PrDens)) %>% 
      filter(cumsum_prdens == max(cumsum_prdens)) %>% 
      pull(calBP) %>% 
      max()
    
    onesigma <- 
      spd_out$grid %>% 
      arrange(., desc(PrDens)) %>% 
      mutate(cumsum_prdens = cumsum(PrDens)) %>% 
      filter(cumsum_prdens <= 0.6827*sum(PrDens))
    
    # ggplot(data = onesigma, 
    #        aes(x = calBP, y = PrDens)) + 
    #   geom_col() +
    #   geom_vline(xintercept = dplyr::pull(filtered_C14_data[which(filtered_C14_data$Site_ID_split == i),"median_SPD_age_calBP"]),
    #              color = "red", size = 2)
    
    filtered_C14_data[which(filtered_C14_data$Site_ID_split == i),"oneSigma_rangeMin"] <- min(onesigma$calBP)
    filtered_C14_data[which(filtered_C14_data$Site_ID_split == i),"oneSigma_rangeMax"] <- max(onesigma$calBP)
    
}

# subset dates to our timeframe
filtered_C14_data_calibrated <- 
  filtered_C14_data %>% 
  select(-c(Key_site,
            "OSL/TL date",
            "Std-dev"#,
            # rank_1_sites_C14_unique_row_ID,
            # Site_ID_split,
            #"Site_ID.y"
            )) %>% 
  unique() %>% 
  filter(., median_SPD_age_calBP <= 16000 &
           median_SPD_age_calBP >= 11000) 

hist(filtered_C14_data_calibrated$median_SPD_age_calBP)

# save
readr::write_csv(filtered_C14_data_calibrated,
                 file.path("3_output", 
                           "C14",
                           "filtered_C14_data_calibrated.csv"))

# select(filtered_C14_data_calibrated,
#        Site_ID_split, Dating_BP, median_SPD_age_calBP, Timeslice, `LabCode/Ref`, Conv_C14_date, AMS_C14_date) %>%
#   unique() %>%
#   # View()
#   readr::write_csv(.,
#                    file.path("3_output", 
#                              "C14",
#                              "filtered_C14_data_calibrated_COMPARISON.csv"))


###########
# visualise
###########
filtered_C14_data_calibrated %>% 
  ggplot(aes(y = reorder(Site_ID_split, -median_SPD_age_calBP)))+
  geom_point(aes(x = median_SPD_age_calBP,
                 color = -median_SPD_age_calBP), 
             size = 2) +
  scale_x_reverse() +
  theme(legend.position = "none")

timerange_per_siteID_fig <- 
  filtered_C14_data_calibrated %>% 
  ggplot(aes(y = reorder(Site_ID_split, -median_SPD_age_calBP)))+
  geom_pointrange(aes(y = reorder(Site_ID_split, -median_SPD_age_calBP), 
                      x = median_SPD_age_calBP,
                      xmin = oneSigma_rangeMax, 
                      xmax = oneSigma_rangeMin,
                      color = -median_SPD_age_calBP)) +
  scale_x_reverse() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Years calBP") +
  ylab("Site ID")

ggsave(timerange_per_siteID_fig,
       path = file.path("3_output",
                        "C14"),
       filename = "timerange_per_siteID_fig.png",
       width = 20,
       height = 30,
       units = "cm")

# filtered_C14_data_calibrated %>% 
#   ggplot(aes(y = reorder(Macro_region, Lat)))+
#   geom_point(aes(x = median_SPD_age_calBP,
#                  color = -median_SPD_age_calBP), 
#              size = 2) +
#   scale_x_reverse() +
#   theme(legend.position = "none") 

filtered_C14_data_calibrated_for_outlines <-
filtered_C14_data_calibrated %>% 
   select(Site_ID_split, Long, Lat,
          median_SPD_age_calBP,
          oneSigma_rangeMin,
          oneSigma_rangeMax
     # -c(Site_ID,
  #           # rank_1_sites_C14_unique_row_ID,
  #           Dating_method,Dating_qual,Timeslice,TaxUnit,
  #           "LabCode/Ref",
  #           "Additional info on layer/context not originally in the datasheet",
  #           "Conv_C14_date","Conv_std",
  #           "AMS_C14_date","AMS_Pm/Std",
  #           "Dated material","Dated species",
  #           "Informaton on pretreatment","Source (Database/Literature)",
  #           "Org. Lit. Reference (if given)","Other comments/observation"
     ) %>% 
  unique() 

readr::write_csv(filtered_C14_data_calibrated_for_outlines,
                 file.path("3_output", 
                           "C14",
                           "filtered_C14_data_calibrated_for_outlines.csv"))

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
               fill = NA, 
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

#################################
# assign dates to climatic events
#################################
data_w_events <- 
  filtered_C14_data_calibrated %>% 
  mutate(Event = case_when(median_SPD_age_calBP >= 14600 ~ "GS-2",
                           median_SPD_age_calBP < 14600 & median_SPD_age_calBP >= 12900~ "GI-1",
                           median_SPD_age_calBP < 12900 & median_SPD_age_calBP >= 11700~ "GS-1",
                           median_SPD_age_calBP < 11700 ~ "Holocene")) %>% 
  mutate(Event = factor(Event, 
                        levels = c("GS-2", "GI-1", "GS-1", "Holocene"))) %>% 
  select(., Event, Long, Lat, median_SPD_age_calBP, Macro_region, TaxUnit_unique, TaxUnit) %>% 
  unique()

# plot sites + dates by events
dates_sites_per_event_fig <- 
base_map +
  geom_jitter(data = data_w_events,
              aes(x = Long, y = Lat,
                  color = median_SPD_age_calBP),
              # alpha = 0.7,
              size = 2) +
  ggrepel::geom_text_repel(data = data_w_events,
                           aes(x = Long, y = Lat,
                               label = median_SPD_age_calBP),
                           # alpha = 0.7,
                           size = 4) +
  facet_wrap(~Event) +
  theme(legend.position = "none")

dates_sites_per_event_fig


ggsave(dates_sites_per_event_fig,
       path = file.path("3_output",
                        "C14"),
       filename = "dates_sites_per_event_fig.png",
       width = 30,
       height = 30,
       units = "cm")








