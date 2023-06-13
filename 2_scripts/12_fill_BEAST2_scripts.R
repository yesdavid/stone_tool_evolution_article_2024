
##################################
# load outlines + PCA data
outlines_centered_scaled_subset_PCA <- readRDS(file = file.path("1_data",
                                                                "Outlines",
                                                                "TAXA88_PCs88_final_subset_outlines_centered_scaled_seed1_PCA.RDS"))

pcs <- outlines_centered_scaled_subset_PCA$x
number_of_pc_axes_used <- ncol(pcs) 
##################################

##################################
chainlength_in_millions <- 30
printgen <- 30000
##################################

##################################
# load taxa file
taxa_file_raw <- readr::read_tsv(file.path("1_data",
                                       "final_subset_outlines_centered_scaled_FAD_LAD_C14_oneSigmaMinMax.tsv"))
number_of_taxa <- nrow(taxa_file_raw)

age_scaler <- 1000


root_age <- 30000 # origin age BP

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
##################################
# with OFFSET
# 'The "offset" parameter should be set to the starting age of the youngest fossil.' https://taming-the-beast.org/tutorials/FBD-tutorial/FBD-tutorial.pdf
# Now, it's set to the youngest _possible_ age.

taxa_file_WOFFSET <- taxa_file_raw

# ATTENTION! the youngest age does not necessarily have to be the one whose age-range ranges to the youngest date overall!
set_offset_WOFFSET <- min(taxa_file_WOFFSET$oneSigma_rangeMin)
taxa_file_WOFFSET$max <- round((taxa_file_WOFFSET$max - set_offset_WOFFSET)/age_scaler, digits = 3)
taxa_file_WOFFSET$oneSigma_rangeMax <- round((taxa_file_WOFFSET$oneSigma_rangeMax - set_offset_WOFFSET)/age_scaler, digits = 3)
taxa_file_WOFFSET$oneSigma_rangeMin <- round((taxa_file_WOFFSET$oneSigma_rangeMin - set_offset_WOFFSET)/age_scaler, digits = 3)

age_offset_WOFFSET <- subset(subset_taxa_WOFFSET, oneSigma_rangeMin == 0)$max

subset_taxa_WOFFSET <- taxa_file_WOFFSET


##################################
##################################
##################################

xml_script_names <- 
  list.files(file.path("2_scripts",
                       "new_xmls"),
             pattern = "_BLANK.xml")

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
                replacement = round((root_age - set_offset_WOFFSET)/age_scaler, digits = 3),
                x = xml_1)
  ##################################
  
  ##################################
  ## DNA_SEQUENCE_PLACEHOLDER
  # <sequence id="1" taxon="TS1_UMag_IBA_LapadosCoelhos4_Gameiroetal2020_AR__pseudo_no_3" totalcount="4" value="N"/>
  DNA_SEQUENCE_PLACEHOLDER_list <- list()
  for(i in 1:nrow(subset_taxa)){
    DNA_SEQUENCE_PLACEHOLDER_list[[i]] <- 
      paste0("<sequence id=\"", i, "\" taxon=\"", subset_taxa[i,]$taxon, "\" totalcount=\"4\" value=\"N\"/>\n")
  }
  cat(unlist(DNA_SEQUENCE_PLACEHOLDER_list))
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
  cat(unlist(TAXON_ID_PLACEHOLDER_list))
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
  # for constant rate; has to be 0 times n.o. TAXA
  xml_1 <- gsub(pattern = "RATECATASSIGN_PLACEHOLDER",
                replacement = paste(rep(0, times= nrow(subset_taxa)), collapse = " "),
                x = xml_1)
  ##################################
  
  ##################################
  # PRIOR_FOSSILAGES_PLACEHOLDER fossil ages
  ## <distribution id="TS4_Sauv_IBM_Filador4_Romanetal2021_AR_na_1_pseudo_no_14.prior" spec="beast.math.distributions.MRCAPrior" tipsonly="true" tree="@TheTree">
  ##  <taxonset id="TS4_Sauv_IBM_Filador4_Romanetal2021_AR_na_1_pseudo_no_14fossilAgeSampl" spec="TaxonSet">
  ##    <taxon idref="TS4_Sauv_IBM_Filador4_Romanetal2021_AR_na_1_pseudo_no_14"/>
  ##  </taxonset>    
  ##  <Uniform id="Uniform.35" lower="0" name="distr" upper="1.179"/>
  ## </distribution>
  PRIOR_FOSSIL_AGES_PLACEHOLDER_list <- list()
  for(i in 1:nrow(subset_taxa)){
    taxon_name <- subset_taxa[i,]$taxon
    PRIOR_FOSSIL_AGES_PLACEHOLDER_list[[i]] <- 
      paste0("<distribution id=\"", taxon_name, ".prior\" spec=\"beast.math.distributions.MRCAPrior\" tipsonly=\"true\" tree=\"@TheTree\">
                <taxonset id=\"", taxon_name,"fossilAgeSampl", "\" spec=\"TaxonSet\">\n
                    <taxon idref=\"", taxon_name, "\"/>\n
                </taxonset>\n
                <Uniform id=\"Uniform.",i, "\" lower=\"", round(subset(subset_taxa, taxon == taxon_name)$oneSigma_rangeMin, digits = 3), #find lower value for 14C date
             "\" name=\"distr\" upper=\"", round(subset(subset_taxa, taxon == taxon_name)$oneSigma_rangeMax, digits = 3), #find upper value for 14C date
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
  for(i in 1:nrow(subset_taxa)){
    taxon_name <- subset_taxa[i,]$taxon
    
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
  for(i in 1:nrow(subset_taxa)){
    taxon_name <- subset_taxa[i,]$taxon
    PRIOR_FOSSIL_AGES_PLACEHOLDER_logger_list[[i]] <- 
      paste0("<log idref=\"", taxon_name, ".prior\"/>\n")
  }
  xml_1 <- 
    gsub(pattern = "LOG_FOSSILAGES_PLACEHOLDER",
         replacement = paste0((unlist(PRIOR_FOSSIL_AGES_PLACEHOLDER_logger_list)), collapse = ""),
         x = xml_1)
  ##################################
  
  ##################################
  # save script
  readr::write_file(x = xml_1,
                    file = file.path("2_scripts",
                                     "new_xmls",
                                     "current_data",
                                     gsub(x = current_script, 
                                          pattern = "_BLANK",
                                          replacement = ""))
  )
}







