library(magrittr)
##################################
# load outlines + PCA data
outlines_centered_scaled_subset_PCA <- readRDS(file = file.path("1_data",
                                                                "Outlines",
                                                                "TAXA88_PCs88_final_subset_outlines_centered_scaled_seed1_PCA.RDS"))

pcs <- outlines_centered_scaled_subset_PCA$x

number_of_artefacts_used <- nrow(pcs)
number_of_pc_axes_used <- ncol(pcs) 
##################################

##################################
chainlength_in_millions <- 30
printgen <- 300000
##################################

##################################
# load taxa file
taxa_file_raw <- readr::read_tsv(file.path("1_data",
                                       "final_subset_outlines_centered_scaled_FAD_LAD_C14_oneSigmaMinMax.tsv"))
number_of_taxa <- nrow(taxa_file_raw)

age_scaler <- 1000

root_age <- 30000 # origin age BP


# for skyline: set birthRateChangeTimes to sensible points in time (i.e. NGRIP event stratigraphy)
RATECHANGETIMES_raw <- c(14600, 12900, 11700)
N_O_RATECHANGETIMES <- length(RATECHANGETIMES_raw)+1

#### modify the ages!

# FBD_constant_rate.xml + FBD_skyline.xml
## one single fixed age for each artefact, no age uncertainty. single fixed ages have to be rescaled so that the youngest age = 0.

# FBD_constant_rate_age_uncertainty.xml + FBD_skyline_age_uncertainty.xml
## one single age for each artefact (derived from the maximum of the Summed Probability Distribution [SPD] of the 14C-dates) plus an age range corresponding to the one-sigma range of the SPD.
## here, the youngest age is set to 0 and all one-sigma range ages are pruned so they cannot be younger than 0.
### For FBD_skyline_age_uncertainty.xml with the youngest age fixed to zero, the zero age specimen won’t have a prior or operator assigned to it!!!

# FBD_constant_rate_age_uncertainty_woffset.xml
## one single age for each artefact (derived from the maximum of the Summed Probability Distribution [SPD] of the 14C-dates) plus an age range corresponding to the one-sigma range of the SPD.
## here, the youngest one-sigma range age is set to 0. all other ages are scaled appropriately. 
## the age offset is now (?) the maximum SPD age of the age whose youngest one-sigma range age is set to 0 (????)

##################################
# NO OFFSET
## scale so that youngest age = 0
### attention! the youngest age may only be = 0 if A) the ages are fixed, or B) there is no age uncertainty associated with this age=0 specimen!
### for fully extinct trees, the minimum age has to be 0+min(oneSigma_rangeMin)!
taxa_file <- taxa_file_raw

age_offset <- min(taxa_file$max)
age_offset_beast <- age_offset/age_scaler # scale age to range 0-10

taxa_file$max <- round((taxa_file$max - age_offset)/age_scaler, digits = 3)
taxa_file$oneSigma_rangeMax <- round((taxa_file$oneSigma_rangeMax - age_offset)/age_scaler, digits = 3)
taxa_file$oneSigma_rangeMin <- round((taxa_file$oneSigma_rangeMin - age_offset)/age_scaler, digits = 3)
taxa_file$oneSigma_rangeMin[which(taxa_file$oneSigma_rangeMin < 0)] <- 0 # minimum sigma range cannot be less than 0!

subset_taxa <- taxa_file

plot_age_uncertainty_youngest_date_zero <- 
  ggplot2::ggplot(data = dplyr::arrange(subset_taxa, max) %>% dplyr::top_n(n=-10), 
                  ggplot2::aes(y = reorder(taxon, -max), 
                               x=max,
                               xmin = oneSigma_rangeMax, 
                               xmax = oneSigma_rangeMin,
                               color = max)) + 
  ggplot2::geom_errorbar(width = 0.5,
                         lty = "dashed") +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_reverse() +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::xlab("Age calBP") +
  ggplot2::ylab("Artefactnames") +
  ggplot2::geom_vline(xintercept = 0,
                      lty = "dashed",
                      color = "red") +
  ggplot2::ggtitle("Starting ages with age uncertainty AND all age ranges\n(age uncertainties cannot be <0)")
##################################
# with OFFSET
# 'The "offset" parameter should be set to the starting age of the youngest fossil.' https://taming-the-beast.org/tutorials/FBD-tutorial/FBD-tutorial.pdf

# starting values for fossil ages = on the “real” timescale
# offset = the youngest fossil age in “real” time
# age range = the 2-sigma range in “real” time

taxa_file_WOFFSET <- taxa_file_raw

taxa_file_WOFFSET$max <- round((taxa_file_WOFFSET$max)/age_scaler, digits = 3)
taxa_file_WOFFSET$oneSigma_rangeMax <- round((taxa_file_WOFFSET$oneSigma_rangeMax)/age_scaler, digits = 3)
taxa_file_WOFFSET$oneSigma_rangeMin <- round((taxa_file_WOFFSET$oneSigma_rangeMin)/age_scaler, digits = 3)

age_offset_WOFFSET <- min(taxa_file_WOFFSET$max)

subset_taxa_WOFFSET <- taxa_file_WOFFSET

plot_age_uncertainty_WOFFSET <- 
  ggplot2::ggplot(data = dplyr::arrange(subset_taxa_WOFFSET, max) %>% dplyr::top_n(n=-10), 
                  ggplot2::aes(y = reorder(taxon, -max), 
                               x= max,
                               xmin = oneSigma_rangeMax, 
                               xmax = oneSigma_rangeMin,
                               color = max)) + 
  ggplot2::geom_errorbar(width = 0.5,
                         lty = "dashed") +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_reverse() +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::xlab("Age calBP") +
  ggplot2::ylab("Artefactnames") +
  ggplot2::geom_vline(xintercept = 0,
                      lty = "dashed",
                      color = "red") +
  ggplot2::ggtitle("Starting ages with age uncertainty\n(for offset-script)")
##################################

##################################
# Youngest age fixed to 0 for SKYLINE script
subset_taxa_skyline <- subset_taxa
subset_taxa_skyline[which(subset_taxa_skyline$max == 0),]
subset_taxa_skyline[which(subset_taxa_skyline$max == 0),"oneSigma_rangeMax"] <- 0

plot_age_uncertainty_youngest_date_fixed <- 
  ggplot2::ggplot(data = dplyr::arrange(subset_taxa_skyline, max) %>% dplyr::top_n(n=-10), 
                  ggplot2::aes(y = reorder(taxon, -max), 
                               x=max,
                               xmin = oneSigma_rangeMax, 
                               xmax = oneSigma_rangeMin,
                               color = max)) + 
  ggplot2::geom_errorbar(width = 0.5,
                         lty = "dashed") +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_reverse() +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::xlab("Age calBP") +
  ggplot2::ylab("Artefactnames") +
  ggplot2::geom_vline(xintercept = 0,
                      lty = "dashed",
                      color = "red") +
  ggplot2::ggtitle("Starting ages with youngest age fixed\n(for skyline model with age uncertainty)")
##################################

##################################
cowplot::plot_grid(plot_age_uncertainty_youngest_date_zero,
                   plot_age_uncertainty_youngest_date_fixed,
                   plot_age_uncertainty_WOFFSET,
                   ncol = 1)

##################################

sh_run <- readr::read_file(file.path("2_scripts",
                                     "new_xmls",
                                     "run_BLANK.sh"))
sh_resume <- readr::read_file(file.path("2_scripts",
                                     "new_xmls",
                                     "resume_BLANK.sh"))
xml_script_names <- 
  list.files(file.path("2_scripts",
                       "new_xmls"),
             pattern = "_BLANK.xml")


number_of_independent_runs <- 4

for(RUN_PLACEHOLDER in 1:number_of_independent_runs){
  
  
  
  CURRENT_FOLDER_PLACEHOLDER <- file.path(paste0("TAXA_", number_of_artefacts_used, "_PCs_", number_of_pc_axes_used),
                                          paste0("independent_run_",RUN_PLACEHOLDER))
  
  for(current_script in xml_script_names){
    xml <- readr::read_file(file.path("2_scripts",
                                      "new_xmls",
                                      current_script))
    
    ##################################
    ## CHAINLENGTH_PLACEHOLDER
    xml_1 <- gsub(pattern = "CHAINLENGTH_PLACEHOLDER",
                  replacement = format(chainlength_in_millions*10^6, scientific = FALSE),
                  x = xml)
    # PRINTGEN_PLACEHOLDER
    xml_1 <- gsub(pattern = "PRINTGEN_PLACEHOLDER",
                  replacement = format(printgen, scientific = FALSE),
                  x = xml_1)
    ##################################
    
    ##################################
    ## ROOT_AGE_PLACEHOLDER
    xml_1 <- gsub(pattern = "ROOT_AGE_PLACEHOLDER",
                  replacement = round((root_age - age_offset)/age_scaler, digits = 3),
                  x = xml_1)
    ## ROOT_AGE_WOFFSET_PLACEHOLDER
    xml_1 <- gsub(pattern = "ROOT_AGE_WOFFSET_PLACEHOLDER",
                  replacement = round((root_age)/age_scaler, digits = 3),
                  x = xml_1)
    ##################################
    
    ##################################
    ## SKYLINE RATECHANGETIMES
    xml_1 <- gsub(pattern = "RATECHANGETIMES_PLACEHOLDER",
                  replacement = paste(c(round((RATECHANGETIMES_raw - age_offset)/age_scaler, digits = 3),0), collapse = " "), #last entry must always be zero
                  x = xml_1)
    
    xml_1 <- gsub(pattern = "BirthDeath_RATES_PLACEHOLDER",
                  replacement = paste(rep("1.0", times = N_O_RATECHANGETIMES), collapse = " "), #last entry must always be zero
                  x = xml_1)
    xml_1 <- gsub(pattern = "Sampling_RATES_PLACEHOLDER",
                  replacement = paste(rep("0.1", times = N_O_RATECHANGETIMES), collapse = " "), #last entry must always be zero
                  x = xml_1)
    
    N_O_RATECHANGETIMES
    ##################################
    
    ##################################
    ## DNA_SEQUENCE_PLACEHOLDER
    # <sequence id="1" taxon="TS1_UMag_IBA_LapadosCoelhos4_Gameiroetal2020_AR__pseudo_no_3" totalcount="4" value="N"/>
    DNA_SEQUENCE_PLACEHOLDER_list <- list()
    for(i in 1:nrow(subset_taxa)){
      DNA_SEQUENCE_PLACEHOLDER_list[[i]] <- 
        paste0("<sequence id=\"", i, "\" taxon=\"", subset_taxa[i,]$taxon, "\" totalcount=\"4\" value=\"N\"/>\n")
    }
    # cat(unlist(DNA_SEQUENCE_PLACEHOLDER_list))
    xml_1 <- gsub(pattern = "DNA_SEQUENCE_PLACEHOLDER",
                  replacement = paste0((unlist(DNA_SEQUENCE_PLACEHOLDER_list)), collapse = ""),
                  x = xml_1)
    ##################################
    
    ##################################
    ## TAXON_ID_PLACEHOLDER
    # <taxon id="TS1_UMag_IBA_LapadosCoelhos4_Gameiroetal2020_AR__pseudo_no_3" spec="Taxon"/>
    TAXON_ID_PLACEHOLDER_list <- list()
    for(i in 1:nrow(subset_taxa)){
      TAXON_ID_PLACEHOLDER_list[[i]] <- 
        paste0("<taxon id=\"", subset_taxa[i,]$taxon,"\" spec=\"Taxon\"/>\n")
    }
    # cat(unlist(TAXON_ID_PLACEHOLDER_list))
    xml_1 <- gsub(pattern = "TAXON_ID_PLACEHOLDER",
                  replacement = paste0((unlist(TAXON_ID_PLACEHOLDER_list)), collapse = ""),
                  x = xml_1)
    ##################################
    
    ##################################
    ## FOSSILSET_PLACEHOLDER
    # TS1_UMag_IBA_LapadosCoelhos4_Gameiroetal2020_AR__pseudo_no_3=2.78,
    FOSSILSET_PLACEHOLDER_list <- list()
    for(i in 1:nrow(subset_taxa)){
      taxon_name <- subset_taxa[i,]$taxon
      if(i != nrow(subset_taxa)){ # add comma at the end of the line as long as it's not the last entry
        FOSSILSET_PLACEHOLDER_list[[i]] <- 
          paste0(taxon_name,
                 # as.character(subset_taxa[i,"taxon"]), 
                 "=", 
                 # as.numeric(subset_taxa[i,"max"]), 
                 round(subset(subset_taxa, taxon == taxon_name)$max, digits = 3),
                 ",\n", 
                 collapse ="")
      } else {
        FOSSILSET_PLACEHOLDER_list[[i]] <- 
          paste0(taxon_name,
                 # as.character(subset_taxa[i,"taxon"]), 
                 "=", 
                 # as.numeric(subset_taxa[i,"max"]),  #no comma at the end of the line, no new line
                 round(subset(subset_taxa, taxon == taxon_name)$max, digits = 3),
                 collapse ="")
      }
    }
    xml_1 <- gsub(pattern = "FOSSILSET_PLACEHOLDER",
                  replacement = paste0((unlist(FOSSILSET_PLACEHOLDER_list)), collapse = ""),
                  x = xml_1)
    
    ## FOSSILSET_WOFFSET_PLACEHOLDER
    FOSSILSET_WOFFSET_PLACEHOLDER_list <- list()
    for(i in 1:nrow(subset_taxa_WOFFSET)){
      taxon_name <- subset_taxa_WOFFSET[i,]$taxon
      if(i != nrow(subset_taxa_WOFFSET)){ # add comma at the end of the line as long as it's not the last entry
        FOSSILSET_WOFFSET_PLACEHOLDER_list[[i]] <- 
          paste0(taxon_name,
                 # as.character(subset_taxa_WOFFSET[i,"taxon"]), 
                 "=", 
                 # as.numeric(subset_taxa_WOFFSET[i,"max"]), 
                 round(subset(subset_taxa_WOFFSET, taxon == taxon_name)$max, digits = 3),
                 ",\n", 
                 collapse ="")
      } else {
        FOSSILSET_WOFFSET_PLACEHOLDER_list[[i]] <- 
          paste0(taxon_name,
                 # as.character(subset_taxa_WOFFSET[i,"taxon"]), 
                 "=", 
                 # as.numeric(subset_taxa_WOFFSET[i,"max"]),  #no comma at the end of the line, no new line
                 round(subset(subset_taxa_WOFFSET, taxon == taxon_name)$max, digits = 3),
                 collapse ="")
      }
    }
    xml_1 <- gsub(pattern = "FOSSILSET_WOFFSET_PLACEHOLDER",
                  replacement = paste0((unlist(FOSSILSET_WOFFSET_PLACEHOLDER_list)), collapse = ""),
                  x = xml_1)
    
    ##################################
    
    ##################################
    ## MINORDIMENSION_PLACEHOLDER <- number of PC-axes
    #number_of_pc_axes_used
    xml_1 <- gsub(pattern = "MINORDIMENSION_PLACEHOLDER",
                  replacement = number_of_pc_axes_used,
                  x = xml_1)
    ##################################
    
    ##################################
    ## MINORDIMENSION_PLACEHOLDER <- number of PC-axes
    xml_1 <- gsub(pattern = "MINORDIMENSION_PLACEHOLDER",
                  replacement = number_of_pc_axes_used,
                  x = xml_1)
    ## ONETRAITDATA_PLACEHOLDER
    xml_1 <- gsub(pattern = "ONETRAITDATA_PLACEHOLDER",
                  replacement = paste0(rownames(pcs), collapse = " "),
                  x = xml_1)
    ## TRAITDATA_PLACEHOLDER
    PCaxis_traits <- c()
    for(i in 1:nrow(subset_taxa)){
      PCaxis_traits <- c(PCaxis_traits, pcs[i,])
    }
    names(PCaxis_traits) <- NULL
    xml_1 <- gsub(pattern = "TRAITDATA_PLACEHOLDER",
                  replacement = paste0(PCaxis_traits, collapse = " "),
                  x = xml_1)
    ##################################
    
    ##################################
    ## FOSSILSET_PLACEHOLDER
    # TS1_UMag_IBA_LapadosCoelhos4_Gameiroetal2020_AR__pseudo_no_3=2.78,
    FOSSILSET_PLACEHOLDER_list <- list()
    for(i in 1:nrow(subset_taxa)){
      taxon_name <- subset_taxa[i,]$taxon
      if(i != nrow(subset_taxa)){ # add comma at the end of the line as long as it's not the last entry
        FOSSILSET_PLACEHOLDER_list[[i]] <- 
          paste0(taxon_name,
                 "=", 
                 round(subset(subset_taxa, taxon == taxon_name)$max, digits = 3),
                 ",\n", 
                 collapse ="")
      } else {
        FOSSILSET_PLACEHOLDER_list[[i]] <- 
          paste0(taxon_name,
                 "=", 
                 round(subset(subset_taxa, taxon == taxon_name)$max, digits = 3),
                 collapse ="")
      }
    }
    xml_1 <- gsub(pattern = "FOSSILSET_PLACEHOLDER",
                  replacement = paste0((unlist(FOSSILSET_PLACEHOLDER_list)), collapse = ""),
                  x = xml_1)
    ##################################
    
    ##################################
    # RATECATASSIGN_PLACEHOLDER 
    # for constant rate; has to be 0 times number of TAXA
    xml_1 <- gsub(pattern = "RATECATASSIGN_PLACEHOLDER",
                  replacement = paste(rep(0, times= nrow(subset_taxa)), collapse = " "),
                  x = xml_1)
    ##################################
    
    
    if(current_script == "FBD_skyline_age_uncertainty_BLANK.xml" ){
      subset_taxa_VARIABLE <- subset_taxa %>% 
        dplyr::filter(!(taxon %in% dplyr::pull(subset_taxa_skyline[which(subset_taxa_skyline$max == 0),"taxon"]))) # remove all information for the taxon which has been fixed to age 0, in order for it to not get assigned any prior, operator, age range or logger
    } else {
      subset_taxa_VARIABLE <- subset_taxa
    }
    
    ##################################
    # PRIOR_FOSSILAGES_PLACEHOLDER fossil ages
    ## <distribution id="TS4_Sauv_IBM_Filador4_Romanetal2021_AR_na_1_pseudo_no_14.prior" spec="beast.math.distributions.MRCAPrior" tipsonly="true" tree="@TheTree">
    ##  <taxonset id="TS4_Sauv_IBM_Filador4_Romanetal2021_AR_na_1_pseudo_no_14fossilAgeSampl" spec="TaxonSet">
    ##    <taxon idref="TS4_Sauv_IBM_Filador4_Romanetal2021_AR_na_1_pseudo_no_14"/>
    ##  </taxonset>    
    ##  <Uniform id="Uniform.35" lower="0" name="distr" upper="1.179"/>
    ## </distribution>
    PRIOR_FOSSIL_AGES_PLACEHOLDER_list <- list()
    for(i in 1:nrow(subset_taxa_VARIABLE)){
      taxon_name <- subset_taxa_VARIABLE[i,]$taxon
      PRIOR_FOSSIL_AGES_PLACEHOLDER_list[[i]] <- 
        paste0("<distribution id=\"", taxon_name, ".prior\" spec=\"beast.math.distributions.MRCAPrior\" tipsonly=\"true\" tree=\"@TheTree\">
                <taxonset id=\"", taxon_name,"fossilAgeSampl", "\" spec=\"TaxonSet\">\n
                    <taxon idref=\"", taxon_name, "\"/>\n
                </taxonset>\n
                <Uniform id=\"Uniform.",i, "\" lower=\"", round(subset(subset_taxa_VARIABLE, taxon == taxon_name)$oneSigma_rangeMin, digits = 3), #find lower value for 14C date
               "\" name=\"distr\" upper=\"", round(subset(subset_taxa_VARIABLE, taxon == taxon_name)$oneSigma_rangeMax, digits = 3), #find upper value for 14C date
               "\"/>\n 
            </distribution>\n"
        )
    }
    xml_1 <- 
      gsub(pattern = "PRIOR_FOSSILAGES_PLACEHOLDER",
           replacement = paste0((unlist(PRIOR_FOSSIL_AGES_PLACEHOLDER_list)), collapse = ""),
           x = xml_1)
    
    #PRIOR_FOSSILAGES_WOFFSET_PLACEHOLDER
    PRIOR_FOSSILAGES_WOFFSET_PLACEHOLDER_list <- list()
    for(i in 1:nrow(subset_taxa_WOFFSET)){
      taxon_name <- subset_taxa_WOFFSET[i,]$taxon
      PRIOR_FOSSILAGES_WOFFSET_PLACEHOLDER_list[[i]] <-
        paste0("<distribution id=\"", taxon_name, ".prior\" spec=\"beast.math.distributions.MRCAPrior\" tipsonly=\"true\" tree=\"@TheTree\">
                <taxonset id=\"", taxon_name,"fossilAgeSampl", "\" spec=\"TaxonSet\">\n
                    <taxon idref=\"", taxon_name, "\"/>\n
                </taxonset>\n
                <Uniform id=\"Uniform.",i, "\" lower=\"", round(subset(subset_taxa_WOFFSET, taxon == taxon_name)$oneSigma_rangeMin, digits = 3), #find lower value for 14C date
               "\" name=\"distr\" upper=\"", round(subset(subset_taxa_WOFFSET, taxon == taxon_name)$oneSigma_rangeMax, digits = 3), #find upper value for 14C date
               "\"/>\n
            </distribution>\n"
        )
    }
    xml_1 <-
      gsub(pattern = "PRIOR_FOSSILAGES_WOFFSET_PLACEHOLDER",
           replacement = paste0((unlist(PRIOR_FOSSILAGES_WOFFSET_PLACEHOLDER_list)), collapse = ""),
           x = xml_1)
    ##################################
    
    ##################################
    # MOVE_FOSSILAGES_PLACEHOLDER moves on fossil ages
    MOVE_FOSSILAGES_PLACEHOLDER_list <- list()
    MOVE_FOSSILAGES_WOFFSET_PLACEHOLDER_list <- list()
    for(i in 1:nrow(subset_taxa_VARIABLE)){
      taxon_name <- subset_taxa_VARIABLE[i,]$taxon
      
      MOVE_FOSSILAGES_PLACEHOLDER_list[[i]] <- 
        paste0("<operator id=\"tipDatesSampler.", taxon_name, "fossilAgeSampl",
               "\" spec=\"SampledNodeDateRandomWalker\" taxonset=\"@", taxon_name, "fossilAgeSampl",
               "\" tree=\"@TheTree\" weight=\"3.0\" windowSize=\"0.003\"/>\n")
      
      MOVE_FOSSILAGES_WOFFSET_PLACEHOLDER_list[[i]] <- 
        paste0("<operator id=\"tipDatesSampler.", taxon_name, "fossilAgeSampl",
               "\" spec=\"SampledNodeDateRandomWalker\" taxonset=\"@", taxon_name, "fossilAgeSampl",
               "\" tree=\"@TheTree\" treeWOffset=\"@treeWOffset\" weight=\"3.0\" windowSize=\"0.003\"/>\n")
    }
    xml_1 <- 
      gsub(pattern = "MOVE_FOSSILAGES_PLACEHOLDER",
           replacement = paste0(unlist(MOVE_FOSSILAGES_PLACEHOLDER_list), collapse = ""),
           x = xml_1)
    xml_1 <- 
      gsub(pattern = "MOVE_FOSSILAGES_WOFFSET_PLACEHOLDER",
           replacement = paste0(unlist(MOVE_FOSSILAGES_WOFFSET_PLACEHOLDER_list), collapse = ""),
           x = xml_1)
    ##################################
    
    ##################################
    # TREEOFFSET_PLACEHOLDER age offset
    xml_1 <-
      gsub(pattern = "TREEOFFSET_PLACEHOLDER",
           replacement = age_offset_WOFFSET,
           x = xml_1)
    ##################################
    
    ##################################
    # LOG_FOSSILAGES_PLACEHOLDER log fossil ages for each taxon
    ## <log idref="TS4_Sauv_IBM_Filador4_Romanetal2021_AR_na_1_pseudo_no_14.prior"/>
    PRIOR_FOSSIL_AGES_PLACEHOLDER_logger_list <- list()
    for(i in 1:nrow(subset_taxa_VARIABLE)){
      taxon_name <- subset_taxa_VARIABLE[i,]$taxon
      PRIOR_FOSSIL_AGES_PLACEHOLDER_logger_list[[i]] <- 
        paste0("<log idref=\"", taxon_name, ".prior\"/>\n")
    }
    xml_1 <- 
      gsub(pattern = "LOG_FOSSILAGES_PLACEHOLDER",
           replacement = paste0((unlist(PRIOR_FOSSIL_AGES_PLACEHOLDER_logger_list)), collapse = ""),
           x = xml_1)
    ##################################
    
    CURRENT_SCRIPT_PLACEHOLDER <- gsub(x = current_script, 
                                       pattern = "_BLANK.xml",
                                       replacement = "")
    
    # modify output path in xml
    xml_1 <- 
      gsub(pattern = "CURRENT_FOLDER_PLACEHOLDER",
           replacement = CURRENT_FOLDER_PLACEHOLDER,
           x = xml_1)
    
    
    # modify .sh files
    ## sh run
    sh_run_1 <- 
      gsub(pattern = "RUN_PLACEHOLDER",
           replacement = RUN_PLACEHOLDER,
           x = sh_run)
    sh_run_1 <- 
      gsub(pattern = "CURRENT_FOLDER_PLACEHOLDER",
           replacement = CURRENT_FOLDER_PLACEHOLDER,
           x = sh_run_1)
    sh_run_1 <- 
      gsub(pattern = "CURRENT_SCRIPT_PLACEHOLDER",
           replacement = CURRENT_SCRIPT_PLACEHOLDER,
           x = sh_run_1)
    ## sh resume
    sh_resume_1 <- 
      gsub(pattern = "RUN_PLACEHOLDER",
           replacement = RUN_PLACEHOLDER,
           x = sh_resume)
    sh_resume_1 <- 
      gsub(pattern = "CURRENT_FOLDER_PLACEHOLDER",
           replacement = CURRENT_FOLDER_PLACEHOLDER,
           x = sh_resume_1)
    sh_resume_1 <- 
      gsub(pattern = "CURRENT_SCRIPT_PLACEHOLDER",
           replacement = CURRENT_SCRIPT_PLACEHOLDER,
           x = sh_resume_1)
    
    
    ##################################
    # save script
    out_path <- file.path("2_scripts",
                          "new_xmls",
                          CURRENT_FOLDER_PLACEHOLDER)
    dir.create(file.path(out_path, "output"),
               recursive = T)
    
    readr::write_file(x = xml_1,
                      file = file.path(out_path,
                                       gsub(x = current_script, 
                                            pattern = "_BLANK",
                                            replacement = "")))
    readr::write_file(x = sh_run_1,
                      file = file.path(out_path,
                                       paste0("run_", gsub(x = current_script, 
                                                           pattern = "_BLANK.xml",
                                                           replacement = ".sh"))))
    readr::write_file(x = sh_resume_1,
                      file = file.path(out_path,
                                       paste0("resume_", gsub(x = current_script, 
                                                              pattern = "_BLANK.xml",
                                                              replacement = ".sh"))))
  }
  
  
}






