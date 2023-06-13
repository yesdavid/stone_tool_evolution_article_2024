
#############################################################################################
#############################################################################################

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

# saveRDS(outlines_centered_scaled_C14,
#         file.path("1_data",
#                   "Outlines",
#                   "outlines_toKeep_centered_scaled_C14metadata.RDS"))

# #############################################################################################
# #############################################################################################
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

Momocs::panel(outlines_centered_scaled_subs,
              fac = "Region")

table(outlines_centered_scaled_subs$fac$Region)

filtered_C14_data_calibrated_for_outlines$Site_ID_split[!(filtered_C14_data_calibrated_for_outlines$Site_ID_split %in% outlines_centered_scaled_subs$fac$Site_ID_split)]

# name of all subsampled artefact outlines to keep!
outlines_centered_scaled_subs$fac$ARTEFACTNAME



# for(i in outlines_centered_scaled_subs$fac$ARTEFACTNAME){
#   
#   file.copy(from = paste0("1_data/Outlines/armatures_to_keep_new_fr/all/",
#                           i,
#                           ".jpg"),
#             to = paste0("1_data/Outlines/final_subset_outlines/",
#                         i,
#                         ".jpg"))
# }