
library(dplyr)
library(magrittr)
library(ggplot2)
library(rcarbon)





datasheet_sites <- readr::read_csv(file = file.path("1_data",
                                                    "MASTERTABLE_v6_all_revised_sites_with_outline_siteID.csv"))

rank_1_sites <- 
datasheet_sites %>% 
  dplyr::filter(Dating_qual != "Problematic")  %>% 
  dplyr::mutate(Site_fau2 = dplyr::case_when(
    Site_fau %in% "Preseved" ~ TRUE, 
    Site_fau %in% "No" ~ FALSE)) %>% 
  dplyr::filter(Site_fau2 == TRUE) %>% 
  dplyr::select("Site_ID", "Key_site","Level_layer_concentration",
                "TaxUnit","TaxUnit_unique",
                "Timeslice","Long","Lat",
                "Dating_method","Dating_qual",
                "Site_type","Site_strat","Ass_pos","Ass_coh",
                "Site_fau2","Site_excav","Key_references",
                "Expert_editor","Macro_region","Macro_region_code","Unique_row_ID")

table(rank_1_sites$Macro_region)
table(rank_1_sites$Site_strat)
table(rank_1_sites$Site_type)
table(rank_1_sites$Dating_qual)
table(rank_1_sites$Dating_method)

head(rank_1_sites)



# load outlines
outlines_AR_centered_w_metric_measures <- readRDS(file = file.path("..",
                                                                       "CLIOARCH_workshop_2_2020",
                                                                       "post_workshop_2021", 
                                                                       "3_output", 
                                                                       "outlines_AR_centered_w_metric_measures.RDS"))

outlines_AR_centered_w_metric_measures <- outlines_AR_centered_w_metric_measures %>%
  Momocs::filter(!Region %in% "CASIA")




# re-calibrate


# load 14C data sheet
C14_v3 <- 
  readr::read_csv(file = file.path("1_data",
                                   "CLIOARCH_14C_datasheet_v4.csv")) %>% 
  unique()

# View(C14_v3)
# names(C14_v3)

# clean 14C spreadsheet
C14_v3_selected <- 
C14_v3 %>% 
  dplyr::select(!c("Site_ID (from main datasheet)","...19")) %>% 
  filter(!(`LabCode/Ref` %in% "?") & !is.na(`LabCode/Ref`))

# join quality sites with 14C data
rank_1_sites_C14 <-
  dplyr::left_join(rank_1_sites,
                   C14_v3_selected,
                   by = c("Key_site" = "Site (from main datasheet)",
                          "Level_layer_concentration" = "Layer/Context (from main datasheet)",
                          "TaxUnit" = "TaxUnit (from main datasheet)")) %>% 
  dplyr::filter(!is.na(Conv_C14_date) | !is.na(AMS_C14_date) | !is.na(`OSL/TL date`) ) %>% # select only the ones that have a dating
  unique() 




rank_1_sites_C14$rank_1_sites_C14_unique_row_ID <- 1:nrow(rank_1_sites_C14)

list_separate_siteID <- list()
for(i in 1:nrow(rank_1_sites_C14)){
  
  current_row <- rank_1_sites_C14[i,]
  
  split_Site_ID <- strsplit(current_row$Site_ID, ",")[[1]]
  
  if(length(split_Site_ID) > 1){
    
    current_Site_ID_list <- list()
    for (current_Site_ID in 1:length(split_Site_ID)){
      
      current_Site_ID_list[[current_Site_ID]] <- 
        data.frame(Site_ID_split = split_Site_ID[[current_Site_ID]],
                   Site_ID = current_row$Site_ID,
                   rank_1_sites_C14_unique_row_ID = current_row$rank_1_sites_C14_unique_row_ID)
    }
    df <- do.call(rbind.data.frame, current_Site_ID_list)
  
    } else{
      
      df <- 
        data.frame(Site_ID_split = strsplit(current_row$Site_ID, ",")[[1]],
                   Site_ID = current_row$Site_ID,
                   rank_1_sites_C14_unique_row_ID = current_row$rank_1_sites_C14_unique_row_ID)
      
  }

  
  list_separate_siteID[[i]] <- df
}
list_separate_siteID_df <- do.call(rbind.data.frame, list_separate_siteID)



# merge
list_separate_siteID_df_separate_siteId <- dplyr::left_join(rank_1_sites_C14,
                                                            list_separate_siteID_df,
                                                            by = c("rank_1_sites_C14_unique_row_ID"))



Momocs::filter(list_separate_siteID_df_separate_siteId ,
               !(Site_ID_split %in% outlines_AR_centered_w_metric_measures$Site))$Site_ID_split %>% 
  unique()

unique_sites_artefacts <- 
Momocs::filter(outlines_AR_centered_w_metric_measures,
               Site %in% list_separate_siteID_df_separate_siteId$Site_ID_split)$fac$Site %>% 
  unique()


filtered_C14_v4 <- 
  list_separate_siteID_df_separate_siteId %>% 
  filter(Site_ID_split %in% unique_sites_artefacts) 


filtered_C14_v4$median_SPD_age_calBP <- NA
for(i in unique(filtered_C14_v4$Site_ID.x)){

    a <-
      dplyr::filter(filtered_C14_v4, Site_ID.x == i)
    
    caldates <-
      rcarbon::calibrate(x=na.omit(c(a$AMS_C14_date, a$Conv_C14_date)),
                         errors=na.omit(c(a$`AMS_Pm/Std`, a$Conv_std)),
                         calCurves='intcal20')
    
    if(nrow(a) > 1){
      spd_out <- spd(caldates,
                     timeRange=c(19000,9000))
    } else {
      spd_out <- caldates
      spd_out$grid <- spd_out$grids$`1`
    }
    
    filtered_C14_v4[which(filtered_C14_v4$Site_ID.x == i),"median_SPD_age_calBP"] <- 
      spd_out$grid %>% 
      arrange(., desc(calBP)) %>%
      mutate(cumsum_prdens = cumsum(PrDens)) %>% 
      filter(cumsum_prdens <= 0.5*sum(PrDens)) %>% 
      filter(cumsum_prdens == max(cumsum_prdens)) %>% 
      pull(calBP) %>% 
      max()
    
}

filtered_C14_v4_calibrated <- 
  filtered_C14_v4 %>% 
  select(-c(Key_site,
            "OSL/TL date",
            "Std-dev",
            rank_1_sites_C14_unique_row_ID,
            # Site_ID_split,
            "Site_ID.y")) %>% 
  unique() %>% 
  filter(., median_SPD_age_calBP <= 16000 &
           median_SPD_age_calBP >= 11000) 

readr::write_csv(filtered_C14_v4_calibrated,
                 file.path("3_output", 
                           "C14",
                           "filtered_C14_v4_calibrated.csv"))

filtered_C14_v4_calibrated %>% 
  ggplot(aes(y = Site_ID.x))+
  geom_point(aes(x = median_SPD_age_calBP), 
             color = "green", 
             size = 2) +
  scale_x_reverse() +
  theme(legend.position = "none") 

filtered_C14_v4_calibrated_for_outlines <- 
filtered_C14_v4_calibrated %>% 
  select(-c(Site_ID.x,
            Unique_row_ID,
            Dating_method,Dating_qual,Timeslice,TaxUnit,
            "LabCode/Ref",
            "Additional info on layer/context not originally in the datasheet",
            "Conv_C14_date","Conv_std",
            "AMS_C14_date","AMS_Pm/Std",
            "Dated material","Dated species",
            "Informaton on pretreatment","Source (Database/Literature)",
            "Org. Lit. Reference (if given)","Other comments/observation")) %>% 
  unique()

readr::write_csv(filtered_C14_v4_calibrated_for_outlines,
                 file.path("3_output", 
                           "C14",
                           "filtered_C14_v4_calibrated_for_outlines.csv"))


  ### distribution map
  world <- rgeos::gBuffer(rworldmap::getMap(resolution = "high"), byid=TRUE, width=0)
  
  # create an extent which spans only the distribution of our samples
  # potentially, the extents have to be manually in-/decreased
  data_extent <- as(raster::extent(min(filtered_C14_v4_calibrated$Long, #minimum longitude
                                       na.rm = T)-3, 
                                   max(filtered_C14_v4_calibrated$Long, #maximum longitude
                                       na.rm = T)+3, 
                                   min(filtered_C14_v4_calibrated$Lat, #minimum latitude
                                       na.rm = T)-3, 
                                   max(filtered_C14_v4_calibrated$Lat, #maximum latidude
                                       na.rm = T)+3), # order: xmin, xmax, ymin, ymax
                    "SpatialPolygons")
  
  sp::proj4string(data_extent) <- sp::CRS(sp::proj4string(world)) # set the coordinate reference system of the data to be the same as the world map.
  
  world_clip <- raster::intersect(world, data_extent) # select only those parts of the world map within our bounding box/extent 
  
  world_clip_f <- fortify(world_clip) # transforms it into a data frame
  
  # plot
  ggplot() +
    geom_polygon(data = world_clip_f, 
                 aes(x = long, 
                     y = lat, 
                     group = group),
                 fill = NA, 
                 colour = "grey") +
    geom_jitter(data = filtered_C14_v4_calibrated,
               aes(x = Long,
                   y = Lat,
                   color = median_SPD_age_calBP,
                   size = scales::rescale(median_SPD_age_calBP, to = c(0,1)))) +
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
  
  


# ########################
# 
# datasheet_sites <- readr::read_csv(file = "./post_workshop_2021/001_datasheet/MASTERTABLE_v6_all_revised_sites_with_outline_siteID.csv")
# datasheet_sites$site_counter <- 1:nrow(datasheet_sites)
# 
# 
# list_separate_datasheet_sites <- list()
# for(i in 1:nrow(datasheet_sites)){
#   
#   current_row <- datasheet_sites[i,]
#   
#   list_separate_datasheet_sites[[i]] <- 
#     data.frame(Site_ID_separate = strsplit(current_row$Site_ID, ",")[[1]],
#                site_counter = current_row$site_counter
#     )
#   
# }
# list_separate_datasheet_sites_df <- do.call(rbind.data.frame, list_separate_datasheet_sites)
# 
# datasheet_sites_separate_df <- dplyr::left_join(list_separate_datasheet_sites_df,
#                                                 datasheet_sites,
#                                                 by = "site_counter")
# 
# joined_datasheet_sites_C14_v3 <-
# dplyr::left_join(
# 
#   C14_v3[,which(colnames(C14_v3) %in% c("C14_ID",
#                                          "Site (from main datasheet)",
#                                          # "Site_ID (from main datasheet)",
#                                          "Layer/Context (from main datasheet)"))],
#   datasheet_sites[,
#                   which(colnames(datasheet_sites) %in% c("Key_site",
#                                          "Level_layer_concentration",
#                                          "Site_ID",
#                                          "site_counter",
#                                          "TaxUnit_unique"))],
#   by = c("Site (from main datasheet)" = "Key_site",
#          # "Site_ID (from main datasheet)" = "Site_ID",
#          "Layer/Context (from main datasheet)" = "Level_layer_concentration"
#          )
#                  )
# 
# length(unique(joined_datasheet_sites_C14_v3$site_counter))
# nrow(datasheet_sites)
# # 
# # View(joined_datasheet_sites_C14_v3)
# # 
# # subset(joined_datasheet_sites_C14_v3, is.na(site_counter)) %>% 
# #   # dplyr::select(#"Site_ID (from main datasheet)",
# #   #               # "Layer/Context (from main datasheet)",
# #   #               "Site (from main datasheet)") %>% 
# #   unique() %>% 
# #   View()
# # 
# # # which outline sites not in c14 sites?
# unique(outlines_AR_centered_w_metric_measures$fac$Site)[!(which(unique(outlines_AR_centered_w_metric_measures$fac$Site) %in%
#                                                                   unique(joined_datasheet_sites_C14_v3$Site_ID)))]
# # 
# # # which c14 sites not in outline sites?
# unique(joined_datasheet_sites_C14_v3$Site_ID)[!(which(unique(joined_datasheet_sites_C14_v3$Site_ID) %in% unique(outlines_AR_centered_w_metric_measures$fac$Site)))]
# 
# 
# 
# 
# joined_datasheet_sites_C14_v3 <- 
#   dplyr::left_join(
#     C14_v3_calibrated %>% 
#       dplyr::select(-`TaxUnit (from main datasheet)`),
#     datasheet_sites_separate_df %>% 
#       dplyr::select(c(Key_site, Level_layer_concentration, TaxUnit_unique, 
#                       Long, Lat, Site_ID_separate,
#                       Site_type, Site_strat, Ass_pos, Ass_coh, Site_fau, Site_excav, Key_references,
#                       Expert_editor, Macro_region, Macro_region_code, Unique_row_ID, Site_ID, quality_score, quality_rank, site_counter)),
#     by = c("Site_ID" = "Site_ID_separate")
#   )
# 
# 
# 
# 
# blubb <- 
# joined_datasheet_sites_C14_v3 %>% 
# subset(!is.na(site_counter) &
#          median_calBP >= 9000 & median_calBP <= 16000) %>% 
#   dplyr::group_by(`Site (from main datasheet)`, `Site_ID (from main datasheet)`, `Layer/Context (from main datasheet)`) %>% 
#   dplyr::arrange(desc(median_calBP)) %>% 
#   dplyr::mutate(FirstDate=dplyr::first(median_calBP),LastDate=dplyr::last(median_calBP),
#          duration = LastDate-FirstDate)
# 
# # View(blubb)
# 
# # blubb %>% 
# #   dplyr::select(`Site (from main datasheet)`, `Site_ID (from main datasheet)`, `Layer/Context (from main datasheet)`,
# #          # dating_method,
# #          # TaxUnit_unique, Long, Lat, Site_type, Site_strat,
# #          # median_calBP,
# #          site_counter,
# #          Macro_region_code,
# #          FirstDate, LastDate, duration) %>% 
# #   dplyr::arrange(desc(FirstDate)) %>% 
# #   unique() %>% 
# #   ggplot() +
# #   geom_segment(aes(x = FirstDate, xend = LastDate, y = reorder(site_counter, FirstDate), yend = reorder(site_counter, FirstDate))) +
# #   scale_x_reverse() +
# #   facet_wrap(~ Macro_region_code,
# #              scales = "free_y")
#   
# # # geom_ridges
# # blubb %>% 
# #   select(`Site (from main datasheet)`, `Site_ID (from main datasheet)`, `Layer/Context (from main datasheet)`,
# #          # dating_method,
# #          # TaxUnit_unique, Long, Lat, Site_type, Site_strat,
# #          # median_calBP,
# #          site_counter,
# #          Macro_region_code,
# #          FirstDate, LastDate, duration) %>% 
# #   arrange(desc(FirstDate)) %>% 
# #   unique() %>% 
# # ggplot(., 
# #        aes(x = FirstDate, 
# #            y = Macro_region_code)) + 
# #   ggridges::geom_density_ridges2() +
# #   scale_x_reverse() 
# 
# ######################################################
# # combine FAD (first appearance date) with outlines. 
# # confidence interval? the dates can probably not be exactly the same for each outline due to how the FBDR process works (Rachel Warnock)
# # give only the oldest date as terminus _post quem_. we don't know when they went "extinct" which is why we cannot give a LAD (last appearance date) (Cathrine Klein)
# 
# unique(outlines_AR_centered_w_metric_measures$fac$Site) %in% blubb$`Site_ID (from main datasheet)`
# 
# # which outline sites are not listed in the c14 sites?
# unique(outlines_AR_centered_w_metric_measures$fac$Site)[(which(!(unique(outlines_AR_centered_w_metric_measures$fac$Site) %in% 
#                                                                    unique(joined_datasheet_sites_C14_v3$Site_ID))))]
# unique(outlines_AR_centered_w_metric_measures$fac$Site)
# 
# # which c14 sites are not listed in the outline sites?
# unique(joined_datasheet_sites_C14_v3$Site_ID)[(which(!(unique(joined_datasheet_sites_C14_v3$Site_ID) %in% 
#                                                                    unique(outlines_AR_centered_w_metric_measures$fac$Site))))]
# # which c14 sites are not listed in the outline sites?
# unique(joined_datasheet_sites_C14_v3$Site_ID)[(which(!(unique(joined_datasheet_sites_C14_v3$Site_ID) %in% 
#                                                          unique(outlines_AR_centered_w_metric_measures$fac$Site))))]
# # ###
# # joined_datasheet_sites_C14_v3_separated_SiteIDs <-
# # list_separate_siteID_df %>%
# #   dplyr::left_join(joined_datasheet_sites_C14_v3 %>%
# #                      dplyr::select(-Site_ID),
# #                    by = "C14_ID")
# # joined_datasheet_sites_C14_v3_separated_SiteIDs <-
# #   list_separate_siteID_df %>%
# #   dplyr::left_join(blubb %>%
# #                      dplyr::select(-Site_ID),
# #                    by = "C14_ID")
# # 
# # unique(joined_datasheet_sites_C14_v3_separated_SiteIDs$Site_ID)[(which(!(unique(joined_datasheet_sites_C14_v3_separated_SiteIDs$Site_ID) %in% 
# #                                                          unique(outlines_AR_centered_w_metric_measures$fac$Site))))]
# # unique(blubb$`Site (from main datasheet)`)[(which(!(unique(blubb$`Site (from main datasheet)`) %in% 
# #                                  unique(outlines_AR_centered_w_metric_measures$fac$Site))))]
# # 
# # unique(outlines_AR_centered_w_metric_measures$fac$Site)[(which(!(unique(outlines_AR_centered_w_metric_measures$fac$Site) %in% 
# #                                                                            unique(joined_datasheet_sites_C14_v3_separated_SiteIDs$Site_ID))))]
# 
# outlines_AR_centered_w_metric_measures_w_c14dates_w_siteData <- 
# 
#   dplyr::left_join(outlines_AR_centered_w_metric_measures$fac,
#                    
#                    blubb %>% 
#                      dplyr::select(Site_ID,C14_ID,"Site (from main datasheet)", "LabCode/Ref",
#                             median_calBP,#sd, 
#                             dating_method, from_calBP, to_calBP,
#                             "Layer/Context (from main datasheet)","Additional info on layer/context not originally in the datasheet",
#                             "Dated species","Informaton on pretreatment","Source (Database/Literature)","Org. Lit. Reference (if given)"),
#                    
#                    by = c("Site" ="Site_ID")) # should be Site (from main datasheet)
# 
# # for each outline, find the oldest date
# unique_artefactnames <- unique(outlines_AR_centered_w_metric_measures_w_c14dates_w_siteData$ARTEFACTNAME)
# outlines_AR_taxa_c14_list <- list()
# for(i in unique_artefactnames){
#   current_artefact <- subset(outlines_AR_centered_w_metric_measures_w_c14dates_w_siteData, 
#                              ARTEFACTNAME == i)
#   
#   current_artefact <- unique(current_artefact[which(current_artefact$median_calBP == max(current_artefact$median_calBP)),]) %>% 
#     dplyr::select(c(from_calBP, to_calBP, median_calBP))
#   
#   current_artefact_max <- subset(current_artefact, median_calBP == max(current_artefact$median_calBP))
#   
#   if(nrow(current_artefact) > 0 & all(current_artefact$median_calBP > 9000 & current_artefact$median_calBP < 16100)){
#     
#     
#     outlines_AR_taxa_c14_list[[i]] <- 
#     unique(  data.frame(taxon = i,
#                         max = current_artefact_max$median_calBP, # max == oldest date; terminus post quem
#                         # max = runif(1, as.numeric(current_artefact$to_calBP), as.numeric(current_artefact$from_calBP) ), #randomized medial_calBP
#                         min = 0,# min == youngest date; extinction date unknown, therefore 0.
#                         
#                         oneSigma_rangeMax = current_artefact_max$from_calBP,
#                         oneSigma_rangeMin = current_artefact_max$to_calBP))
#     
#     
#   }
#   
# 
# }
# 
# # outlines_AR_taxa_c14_df[,1][which(!(outlines_AR_taxa_c14_df[,1] %in% na.omit(do.call(rbind.data.frame, outlines_AR_taxa_c14_list))[,1]))]
# 
# outlines_AR_taxa_c14_df <- unique(na.omit(do.call(rbind.data.frame, outlines_AR_taxa_c14_list)))
# taxa_file <- outlines_AR_taxa_c14_df[which(is.finite(outlines_AR_taxa_c14_df$max) == T),]
# 
# 
# min(taxa_file$max)
# max(taxa_file$max)
# hist(outlines_AR_taxa_c14_df$max)
# 
# # taxa_file$min <- taxa_file$min - 
# # taxa_file$max <- round(taxa_file$max - min(taxa_file$max), digits = 1)
# # taxa_file[which(unique(taxa_file$max) == T),]
# readr::write_tsv(taxa_file,
#                  path = file.path("post_workshop_2021", "001_data_14c", "outlines_AR_taxa_14c_9k16k_v2.tsv"))
# readr::write_tsv(taxa_file,
#                  path = file.path("post_workshop_2021", "4_revbayes_HPC", "data", "outlines_AR_taxa_14c_9k16k_v2.tsv"))
# 
# 
# # 14c sampling bias per region!
# outlines_AR_centered_w_metric_measures_w_c14dates_w_siteData
# 
# base_map +
#   geom_jitter(#data = subset(outlines_AR_centered_w_metric_measures_w_c14dates_w_siteData, median_calBP <= 12000 & median_calBP >= 10000),
#               data = outlines_AR_centered_w_metric_measures_w_c14dates_w_siteData,
#               aes(x = Long,
#                   y = Lat),
#               color = "black",
#               shape = 4,
#               alpha = 0.3,
#               width = 0.1,
#               height = 0.1) +
#   geom_density_2d_filled(data=outlines_AR_centered_w_metric_measures_w_c14dates_w_siteData,
#                          mapping=aes(x=Long,
#                                      y=Lat),
#                          alpha=0.3,
#                          contour_var = "count",
#                          n = 40,
#                          adjust = 0.75,
#                          show.legend = F) +
#   facet_wrap(Timeslice_range_from~Timeslice_range_to)






