# xml helper function



xml_helper_function <- # has to be loaded first
  function(fossil_age_uncertainty,
           fully_extinct,
           skyline_BDMM,
           underPrior,
           chainlength_in_millions,
           printgen,
           walltime_spec,
           blank_file_path){
    
    ####### fossil ages
    
    ## scale so that youngest age = 0
    ### attention! the youngest age may only be = 0 if A) the ages are fixed, or B) there is no age uncertainty associated with this age=0 specimen!
    ### for fully extinct trees, the minimum age has to be 0+oneSigma_rangeMin !
    if(fully_extinct == F){
      age_offset <- min(taxa_file_raw$max)
      
    } else if(fully_extinct == T){
      
      age_offset <- min(taxa_file_raw$oneSigma_rangeMin)
      
    }
    
    taxa_file <- taxa_file_raw
    
    taxa_file$max <- taxa_file$max - age_offset
    taxa_file$oneSigma_rangeMax <- taxa_file$oneSigma_rangeMax - age_offset
    taxa_file$oneSigma_rangeMin <- taxa_file$oneSigma_rangeMin - age_offset
    
    taxa_file_subset <-
      subset(taxa_file,
             taxon %in% outlines_AR_subset_PCA$fac$ARTEFACTNAME) %>% 
      dplyr::distinct()
    
    # scale age to range 0-10
    age_offset_beast <- age_offset/1000
    print(age_offset_beast)
    taxa_file_subset$max <- round(taxa_file_subset$max/1000, digits = 3)
    taxa_file_subset$oneSigma_rangeMax <- round(taxa_file_subset$oneSigma_rangeMax/1000, digits = 3)
    taxa_file_subset$oneSigma_rangeMin <- round(taxa_file_subset$oneSigma_rangeMin/1000, digits = 3)
    
    number_of_taxa <- nrow(taxa_file_subset)
    subset_taxa <- taxa_file_subset
    
    
    
    
    # input file will be chosen accordingly
    if(underPrior == TRUE){
      if(fossil_age_uncertainty == FALSE & skyline_BDMM == FALSE){
        blank_file_name <- "BMPruneLikelihood_calval_1_FBDbds_underPrior"  # fbds with without age uncertainty
      } else if(fossil_age_uncertainty == FALSE & skyline_BDMM == TRUE){
        blank_file_name <- "BMPruneLikelihood_calval_1_FBDbds_BDMMprime_underPrior_blank" # BDMM-prime skyline without age uncertainty
      } else if(fossil_age_uncertainty == TRUE & skyline_BDMM == FALSE){
        blank_file_name <- "BMPruneLikelihood_calval_1_ageUncertainty_underPrior" # fbds with with age uncertainty
      } else if(fossil_age_uncertainty == TRUE & skyline_BDMM == TRUE){
        blank_file_name <- "BMPruneLikelihood_calval_1_BDMMprime_ageUncertainty_underPrior"  # BDMM-prime skyline with age uncertainty
      }
    } else if(underPrior == FALSE){
      if(fossil_age_uncertainty == FALSE & skyline_BDMM == FALSE){
        blank_file_name <- "BMPruneLikelihood_calval_1_FBDbds" # fbds with without age uncertainty
      } else if(fossil_age_uncertainty == FALSE & skyline_BDMM == TRUE){
        blank_file_name <- "BMPruneLikelihood_calval_1_FBDbds_BDMMprime_blank" # BDMM-prime skyline without age uncertainty
      } else if(fossil_age_uncertainty == TRUE & skyline_BDMM == FALSE){
        blank_file_name <- "BMPruneLikelihood_calval_1_ageUncertainty" # fbds with with age uncertainty
      } else if(fossil_age_uncertainty == TRUE & skyline_BDMM == TRUE){
        blank_file_name <- "BMPruneLikelihood_calval_1_BDMMprime_ageUncertainty" # BDMM-prime skyline with age uncertainty
      }
    }
    
    # input file will be read
    xml <- readr::read_file(file.path(blank_file_path,
                                      paste0(blank_file_name,
                                             "_blank.xml")))
    
    # output names will be adjusted if fullExtinct
    if(fully_extinct == TRUE){
      blank_file_name <- paste0(blank_file_name, "_fullyExtinct")
    }
    
    
    
    
    ## CHAINLENGTH_PLACEHOLDER
    xml_1 <- gsub(pattern = "CHAINLENGTH_PLACEHOLDER",
                  replacement = format(chainlength_in_millions*10^6, scientific = FALSE),
                  x = xml)
    
    # PRINTGEN_PLACEHOLDER
    xml_1 <- gsub(pattern = "PRINTGEN_PLACEHOLDER",
                  replacement = format(printgen, scientific = FALSE),
                  x = xml_1)
    
    
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
    
    
    
    ## TYPETRAITSET_PLACEHOLDER
    # value="CHICKEN_HEBEI_326_2005=NOT_SET,CHICKEN_HONGKONG_915_1997=NOT_SET,...."
    TYPETRAITSET_PLACEHOLDER_list <- list()
    for(i in 1:nrow(subset_taxa)){
      
      TYPETRAITSET_PLACEHOLDER_list[[i]] <- 
        paste0(subset_taxa[i,]$taxon,"=NOT_SET")
      
    }
    cat(unlist(TYPETRAITSET_PLACEHOLDER_list))
    
    xml_1 <- gsub(pattern = "TYPETRAITSET_PLACEHOLDER",
                  replacement = paste0((unlist(TYPETRAITSET_PLACEHOLDER_list)), collapse = ","),
                  x = xml_1)
    
    
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
    
    
    ## MINORDIMENSION_PLACEHOLDER <- number of PC-axes
    number_of_pc_axes_used
    
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
    # PCaxis_traits
    
    names(PCaxis_traits) <- NULL
    # PCaxis_traits
    
    xml_1 <- gsub(pattern = "TRAITDATA_PLACEHOLDER",
                  replacement = paste0(PCaxis_traits, collapse = " "),
                  x = xml_1)
    
    ###############################
    # fossil age uncertainty https://taming-the-beast.org/tutorials/FBD-tutorial/FBD-tutorial.pdf section 3.2.10
    ###############################
    
    
    if(fossil_age_uncertainty == TRUE){
      
      if(fully_extinct == F){
        subset_taxa <- subset(subset_taxa, !(taxon == subset(subset_taxa, max == 0)$taxon))
      }
      
      # prior fossil ages
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
      # cat(unlist(PRIOR_FOSSIL_AGES_PLACEHOLDER_list))
      xml_1 <- 
        gsub(pattern = "<!--PRIOR_FOSSIL_AGES_PLACEHOLDER-->",
             replacement = paste0((unlist(PRIOR_FOSSIL_AGES_PLACEHOLDER_list)), collapse = ""),
             x = xml_1)
      
      
      
      ## operator
      PRIOR_FOSSIL_AGES_PLACEHOLDER_operator_list <- list()
      for(i in 1:nrow(subset_taxa)){
        
        taxon_name <- subset_taxa[i,]$taxon
        
        if(fully_extinct == FALSE){
          PRIOR_FOSSIL_AGES_PLACEHOLDER_operator_list[[i]] <- 
            paste0("<operator id=\"tipDatesSampler.", taxon_name, "fossilAgeSampl",
                   "\" spec=\"SampledNodeDateRandomWalker\" taxonset=\"@", taxon_name, "fossilAgeSampl",
                   "\" tree=\"@TheTree\" weight=\"1.0\" windowSize=\"0.1\"/>\n")
        } else if(fully_extinct == TRUE) { # if fully_extinct == TRUE, insert TREEWOFFSET_PLACEHOLDER_operator: treeWOffset=\"@treeWOffset\"
          PRIOR_FOSSIL_AGES_PLACEHOLDER_operator_list[[i]] <- 
            paste0("<operator id=\"tipDatesSampler.", taxon_name, "fossilAgeSampl",
                   "\" spec=\"SampledNodeDateRandomWalker\" taxonset=\"@", taxon_name, "fossilAgeSampl",
                   "\" tree=\"@TheTree\" treeWOffset=\"@treeWOffset\" weight=\"1.0\" windowSize=\"0.1\"/>\n")
        }
        
        
      }
      # cat(unlist(PRIOR_FOSSIL_AGES_PLACEHOLDER_operator_list))
      xml_1 <- 
        gsub(pattern = "<!--PRIOR_FOSSIL_AGES_PLACEHOLDER_operator-->",
             replacement = paste0((unlist(PRIOR_FOSSIL_AGES_PLACEHOLDER_operator_list)), collapse = ""),
             x = xml_1)
      
      
      ## logger
      PRIOR_FOSSIL_AGES_PLACEHOLDER_logger_list <- list()
      for(i in 1:nrow(subset_taxa)){
        
        taxon_name <- subset_taxa[i,]$taxon
        
        PRIOR_FOSSIL_AGES_PLACEHOLDER_logger_list[[i]] <- 
          paste0("<log idref=\"", taxon_name, ".prior\"/>\n")
        
      }
      # cat(unlist(PRIOR_FOSSIL_AGES_PLACEHOLDER_logger_list))
      xml_1 <- 
        gsub(pattern = "<!--PRIOR_FOSSIL_AGES_PLACEHOLDER_logger-->",
             replacement = paste0((unlist(PRIOR_FOSSIL_AGES_PLACEHOLDER_logger_list)), collapse = ""),
             x = xml_1)
      
      
    }
    
    
    ###############################
    # fully extinct clades https://taming-the-beast.org/tutorials/FBD-tutorial/FBD-tutorial.pdf section 5
    ###############################
    # delete "comment" to "activate" fully extinct in script
    
    if(fully_extinct == TRUE){
      xml_1 <- 
        gsub(pattern = "<!--TREEWOFFSET_PLACEHOLDER\n",
             replacement = "",
             x = xml_1)
      xml_1 <- 
        gsub(pattern = "TREEWOFFSET_PLACEHOLDER-->",
             replacement = "",
             x = xml_1)
      
      # add logger 
      xml_1 <- 
        gsub(pattern = "<!--TREEWOFFSET_PLACEHOLDER_logger-->",
             replacement = "<log id=\"offset\" spec=\"beast.evolution.tree.OffsetLogger\" treeWOffset=\"@treeWOffset\"/>",
             x = xml_1)
      
    }
    
    
    
    ###############################
    # finishing
    ###############################
    # create folder for this constellation
    current_output_path<- file.path(getwd(),"post_workshop_2021","2_script","BEAST2_contraband",
                                    paste0("TAXA",number_of_taxa,
                                           "_PCs",number_of_pc_axes_used))
    
    current_output_path_output <- file.path(current_output_path, 
                                            "output")
    
    dir.create(path = current_output_path_output,
               recursive = T)
    
    
    
    # output names and paths
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.log",
                  replacement = file.path(".", 
                                          paste0("TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used), 
                                          "output",
                                          paste0("out_",
                                                 blank_file_name,
                                                 "_TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used, 
                                                 "_nIter", chainlength_in_millions, "m",
                                                 ".log")
                  ),
                  x = xml_1)
    
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.trees",
                  replacement = file.path(".", 
                                          paste0("TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used), 
                                          "output",
                                          paste0("out_",
                                                 blank_file_name,
                                                 "_TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used, 
                                                 "_nIter", chainlength_in_millions, "m",
                                                 ".trees")
                  ),
                  x = xml_1)
    
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.typed.trees",
                  replacement = file.path(".", 
                                          paste0("TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used), 
                                          "output",
                                          paste0("out_",
                                                 blank_file_name,
                                                 "_TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used, 
                                                 "_nIter", chainlength_in_millions, "m",
                                                 ".typed.trees")
                  ),
                  x = xml_1)
    
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.typed.node.trees",
                  replacement = file.path(".", 
                                          paste0("TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used), 
                                          "output",
                                          paste0("out_",
                                                 blank_file_name,
                                                 "_TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used, 
                                                 "_nIter", chainlength_in_millions, "m",
                                                 ".typed.node.trees")
                  ),
                  x = xml_1)
    
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.traj",
                  replacement = file.path(".", 
                                          paste0("TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used), 
                                          "output",
                                          paste0("out_",
                                                 blank_file_name,
                                                 "_TAXA",number_of_taxa,
                                                 "_PCs",number_of_pc_axes_used, 
                                                 "_nIter", chainlength_in_millions, "m",
                                                 ".traj")
                  ),
                  x = xml_1)
    
    
    
    
    
    ###############################
    # save script
    readr::write_file(x = xml_1,
                      file = file.path(current_output_path,
                                       paste0("script_",
                                              blank_file_name,
                                              "_TAXA", number_of_taxa, 
                                              "_PCs",number_of_pc_axes_used,
                                              "_nIter", chainlength_in_millions, "m",
                                              ".xml"))
    )
    
    
    
    
    ###############################
    # modify .sh file for HPC
    
    job_sh <- readr::read_file(file.path(blank_file_path,
                                         "run_contraband_blank.sh"))
    
    # WALLTIME_PLACEHOLDER 00:00:00 (hh:mm:ss)
    job_sh_1 <- gsub(pattern = "WALLTIME_PLACEHOLDER",
                     replacement = walltime_spec,
                     x = job_sh)
    
    # JOBNAME_PLACEHOLDER
    job_sh_1 <- gsub(pattern = "JOBNAME_PLACEHOLDER",
                     replacement = paste0(blank_file_name,
                                          "_TAXA", number_of_taxa, 
                                          "_PCs",number_of_pc_axes_used,
                                          "_nIter", chainlength_in_millions, "m"),
                     x = job_sh_1)
    
    # SCRIPTPATH_PLACEHOLDER
    job_sh_1 <- gsub(pattern = "SCRIPTPATH_PLACEHOLDER",
                     replacement = file.path(".",
                                             paste0("TAXA",number_of_taxa,
                                                    "_PCs",number_of_pc_axes_used),
                                             paste0("script_",
                                                    blank_file_name,
                                                    "_TAXA", number_of_taxa, 
                                                    "_PCs",number_of_pc_axes_used,
                                                    "_nIter", chainlength_in_millions, "m",
                                                    ".xml")),
                     x = job_sh_1)
    
    
    
    
    # write job.sh
    readr::write_file(x = job_sh_1,
                      file = file.path(blank_file_path,
                                       paste0("job_",
                                              blank_file_name,
                                              "_TAXA", number_of_taxa, 
                                              "_PCs",number_of_pc_axes_used,
                                              "_nIter", chainlength_in_millions, "m",
                                              ".sh")
                      )
    )
    
    
    
    
  }
