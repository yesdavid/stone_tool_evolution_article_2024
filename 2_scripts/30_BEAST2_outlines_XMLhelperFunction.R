# xml helper function



xml_helper_function <- # has to be loaded first
  function(taxa_file, # age has to be in column called "max"
           number_of_pc_axes_used, # number of axes chosen
           pcs, # pca scores
           root_age, 
           clockmodel,
           fossil_age_uncertainty,
           fully_extinct,
           skyline_BDMM,
               timebins, # this helper function does not work for timebins <2. Has to be adjusted manually.
               estimate_changeTimes, # logical, if false, provide the following parameters:
               changeTimes, # the date(s) when the timebins change; has to be of length(timebins-1); has to be in the same format as the raw dates provided in taxa_file_raw
               birthParameter,
               deathParameter, 
               samplingParameter, 
               removalParameter,
           substitution_tree,
           BDS_ExponentialMean,
           SteppingStone,
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
      age_offset <- min(taxa_file$max)
      
    } else if(fully_extinct == T){
      
      age_offset <- min(taxa_file$oneSigma_rangeMin)
      
    }
    
    # taxa_file <- taxa_file_raw
    
    taxa_file$max <- taxa_file$max - age_offset
    taxa_file$oneSigma_rangeMax <- taxa_file$oneSigma_rangeMax - age_offset
    taxa_file$oneSigma_rangeMin <- taxa_file$oneSigma_rangeMin - age_offset
    
    # taxa_file_subset <-
    #   subset(taxa_file,
    #          taxon %in% outlines_AR_subset_PCA$fac$ARTEFACTNAME) %>% 
    #   dplyr::distinct()
    
    # scale age to range 0-10
    age_scaler <- 1000
    age_offset_beast <- age_offset/age_scaler
    print(age_offset_beast)
    taxa_file$max <- round(taxa_file$max/age_scaler, digits = 3)
    taxa_file$oneSigma_rangeMax <- round(taxa_file$oneSigma_rangeMax/age_scaler, digits = 3)
    taxa_file$oneSigma_rangeMin <- round(taxa_file$oneSigma_rangeMin/age_scaler, digits = 3)
    
    number_of_taxa <- nrow(taxa_file)
    subset_taxa <- taxa_file
    
    
    
    
    # input file will be chosen accordingly
      if(fossil_age_uncertainty == FALSE & skyline_BDMM == FALSE){
        blank_file_name <- "BMPruneLikelihood_clockModel_FBDbds" # fbds with without age uncertainty
      } else if(fossil_age_uncertainty == FALSE & skyline_BDMM == TRUE){
        blank_file_name <- "BMPruneLikelihood_clockModel_FBDbds_BDMMprime" # BDMM-prime skyline without age uncertainty
      } else if(fossil_age_uncertainty == TRUE & skyline_BDMM == FALSE){
        blank_file_name <- "BMPruneLikelihood_clockModel_ageUncertainty" # fbds with with age uncertainty
      } else if(fossil_age_uncertainty == TRUE & skyline_BDMM == TRUE){
        blank_file_name <- "BMPruneLikelihood_clockModel_ageUncertainty_BDMMprime" # BDMM-prime skyline with age uncertainty
      }
    
    # input file will be read
    xml <- readr::read_file(file.path(blank_file_path,
                                      paste0(blank_file_name,
                                             "_blank.xml")))
    
    # add beast age offset to output name
    blank_file_name <- paste0(blank_file_name, "_Offset",age_offset_beast)
    
    # output names will be adjusted if fullExtinct
    if(fully_extinct == TRUE){
      blank_file_name <- paste0(blank_file_name, "_fullyExtinct")
    }
    
    
    
    
    ## CHAINLENGTH_PLACEHOLDER
    xml_1 <- gsub(pattern = "CHAINLENGTH_PLACEHOLDER",
                  replacement = format(chainlength_in_millions*10^6, scientific = FALSE),
                  x = xml)
    
    
    ## ROOT_AGE_PLACEHOLDER
    xml_1 <- gsub(pattern = "ROOT_AGE_PLACEHOLDER",
                  replacement = round((root_age - age_offset)/age_scaler, digits = 3),
                  x = xml_1)
    
    
    # PRINTGEN_PLACEHOLDER
    xml_1 <- gsub(pattern = "PRINTGEN_PLACEHOLDER",
                  replacement = format(printgen, scientific = FALSE),
                  x = xml_1)
    
    
    # RHO_PLACEHOLDER
    if(fossil_age_uncertainty == TRUE & fully_extinct == FALSE){
      xml_1 <- gsub(pattern = "RHO_PLACEHOLDER",
                    replacement = 0.00000001,
                    x = xml_1)
    } else if(fossil_age_uncertainty == TRUE & fully_extinct == TRUE){
      xml_1 <- gsub(pattern = "RHO_PLACEHOLDER",
                    replacement = 1,
                    x = xml_1)
    }
    
    
    # sampleFromPrior="true"
    if(underPrior == FALSE){
      xml_1 <- gsub(pattern = " sampleFromPrior=\"true\"",
                    replacement = "",
                    x = xml_1)
    }
    
    
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
    # CLOCKMODEL_PLACEHOLDER
    ############################### 
    
    ############################### strict clock
    if(clockmodel == "strict") {
      blank_file_name <- gsub(x = blank_file_name,
                              pattern = "clockModel",
                              replacement = paste0(clockmodel, "Clock"))
      # <!-- the morphological clock model used here is equivalent to a strict clock -->
      #   <!-- the rates input here refers to the branch rates -->
      xml_1 <- gsub(pattern = "CLOCKMODEL_PLACEHOLDER_model",
                    replacement = "<branchRateModel id=\"rateCatClock\" spec=\"contraband.clock.RateCategoryClockModel\" nCat=\"1\">\n
                                    <rates id=\"rateValues\" spec=\"parameter.RealParameter\">1</rates>\n
                                    <rateCatAssign id=\"rateAssignments\" spec=\"parameter.IntegerParameter\" lower=\"0\" upper=\"1\">0</rateCatAssign>\n
                                    <tree idref=\"TheTree\"/>\n
                                   </branchRateModel>\n",
                    x = xml_1)
      
      # logger
      xml_1 <- gsub(pattern = "<!--CLOCKMODEL_PLACEHOLDER_loggers-->",
                    replacement = 
                      "<log idref=\"rateValues\"/>",
                    x = xml_1)
      
      # CLOCKMODEL_PLACEHOLDER_outputTreeLog
      xml_1 <- gsub(pattern = "CLOCKMODEL_PLACEHOLDER_outputTreeLog",
                    replacement = "<log id=\"TreeWithMetaDataLogger\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@TheTree\"/>",
                    x = xml_1)
      
      if(substitution_tree == TRUE){
        # TREELOG_SUBSTITUTION_PLACEHOLDER
        xml_1 <- gsub(pattern = "<!--TREELOG_SUBSTITUTION_PLACEHOLDER",
                      replacement = "",
                      x = xml_1)
        xml_1 <- gsub(pattern = "TREELOG_SUBSTITUTION_PLACEHOLDER-->",
                      replacement = "",
                      x = xml_1)
        xml_1 <- gsub(pattern = "CLOCKMODEL_PLACEHOLDER_substiOutputTreeLog",
                      replacement = "<log id=\"TreeWithMetaDataLoggerSubstitution\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" substitutions=\"true\" tree=\"@TheTree\"/>",
                      x = xml_1)
      }

      
      ############################### nCat clock
    } else if (clockmodel == "nCat") {
      blank_file_name <- gsub(x = blank_file_name,
                              pattern = "clockModel",
                              replacement = paste0(clockmodel, "Clock"))
      #state node
      xml_1 <- gsub(pattern = "<!--CLOCKMODEL_PLACEHOLDER_stateNode-->",
                    replacement = 
                      "<stateNode idref=\"rateAssignments\"/>\n
                      <stateNode idref=\"rateValues\"/>",
                    x = xml_1)
      
      # from https://github.com/fkmendes/contraband/blob/master/examples/testing/OUMVNLikelihoodOneTrait_FBDTree_RateCatClock.xml lines 103-108
      xml_1 <- gsub(pattern = "CLOCKMODEL_PLACEHOLDER_model",
                    replacement = 
                      paste0("<branchRateModel id=\"rateCatClock\" spec=\"contraband.clock.RateCategoryClockModel\" nCat=\"2\">\n
                      <rates id=\"rateValues\" spec=\"parameter.RealParameter\" lower=\"-Infinity\" upper=\"Infinity\">0.1 0.5</rates>\n
                      <rateCatAssign id=\"rateAssignments\" spec=\"parameter.IntegerParameter\" lower=\"0\" upper=\"1\">", paste(rep(0, times= nrow(pcs)), collapse = " "),"</rateCatAssign>\n
                      <tree idref=\"TheTree\"/>\n
                      </branchRateModel>\n",
                             collapse = " "),
                    x = xml_1)

      # operator
      xml_1 <- gsub(pattern = "<!--CLOCKMODEL_PLACEHOLDER_operators-->",
                    replacement = 
      "<operator id=\"RateAssignmentWalker\" spec=\"IntRandomWalkOperator\" parameter=\"@rateAssignments\" windowSize=\"1\" weight=\"10.0\"/>\n
        <operator id=\"RateValueScaler\" spec=\"ScaleOperator\" parameter=\"@rateValues\" scaleFactor=\"0.75\" weight=\"3.0\"/>",
      x = xml_1)
      
      
      
      # logger
      xml_1 <- gsub(pattern = "<!--CLOCKMODEL_PLACEHOLDER_loggers-->",
                    replacement = 
                      "<log idref=\"rateValues\"/>\n
                      <log idref=\"rateAssignments\"/>",
                    x = xml_1)

      
      # CLOCKMODEL_PLACEHOLDER_outputTreeLog
      xml_1 <- gsub(pattern = "CLOCKMODEL_PLACEHOLDER_outputTreeLog",
                    replacement = 
                    "<log id=\"TreeWithMetaDataLogger\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@TheTree\" branchratemodel=\"@rateCatClock\" sort=\"false\"/>",
                    x = xml_1)
      
      
      
      if(substitution_tree == TRUE){
        # TREELOG_SUBSTITUTION_PLACEHOLDER
        xml_1 <- gsub(pattern = "<!--TREELOG_SUBSTITUTION_PLACEHOLDER",
                      replacement = "",
                      x = xml_1)
        xml_1 <- gsub(pattern = "TREELOG_SUBSTITUTION_PLACEHOLDER-->",
                      replacement = "",
                      x = xml_1)
        xml_1 <- gsub(pattern = "CLOCKMODEL_PLACEHOLDER_substiOutputTreeLog",
                      replacement = 
                        "<log id=\"TreeWithMetaDataLoggerSubstitution\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" substitutions=\"true\" tree=\"@TheTree\" branchratemodel=\"@rateCatClock\" sort=\"false\"/>",
                      x = xml_1)
      }
      
      
      ############################### relaxed ln clock
    } else if (clockmodel == "relaxed") {
      # from https://github.com/fkmendes/contraband/blob/master/examples/testing/Carnivora_Morph_BDSS.xml
      # needs beast2 package "orc"
      
      blank_file_name <- gsub(x = blank_file_name,
                              pattern = "clockModel",
                              replacement = paste0(clockmodel, "Clock"))
      
      #state node
      xml_1 <- gsub(pattern = "<!--CLOCKMODEL_PLACEHOLDER_stateNode-->",
                    replacement = 
                      "<stateNode id=\"branchRate\" spec=\"parameter.RealParameter\" value=\"1.0\"/>\n
                      <stateNode id=\"ucldMean\" spec=\"parameter.RealParameter\" value=\"0.0005\"/>\n
                      <stateNode id=\"ucldStdev\" spec=\"parameter.RealParameter\" lower=\"0.0\" value=\"0.5\"/>",
                    x = xml_1)
      
      # priors
      xml_1 <- gsub(pattern = "<!--CLOCKMODEL_PLACEHOLDER_priors-->",
                    replacement = 
                      "<distribution id=\"rateValuesPrior\" spec=\"beast.math.distributions.Prior\" x=\"@branchRate\">\n
                      <distr id=\"LogNormal.rateValues\" spec=\"beast.math.distributions.LogNormalDistributionModel\" S=\"@ucldStdev\" M=\"1.0\" meanInRealSpace=\"true\"/>\n
                      </distribution>\n
                      \n
                      <prior id=\"ucldStdevPrior\"  name=\"distribution\" x=\"@ucldStdev\">\n
                      <distr id=\"LogNormal.ucldstdev\" spec=\"beast.math.distributions.LogNormalDistributionModel\" S=\"0.3417393\" M=\"0.9403022\" meanInRealSpace=\"true\"/>\n
                      </prior>\n
                      \n
                      <distribution id=\"ucldMeanPrior\" spec=\"beast.math.distributions.Prior\" x=\"@ucldMean\">\n
                      <Gamma id=\"Gamma.ucldmean\" mode=\"ShapeRate\" name=\"distr\">\n
                      <parameter id=\"RealParameter.3\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"alpha\">2.0</parameter>\n
                      <parameter id=\"RealParameter.4\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"beta\">5.0</parameter>\n
                      </Gamma>\n
                      </distribution>\n",
                    x = xml_1)
      
      # actual clockmodel
      xml_1 <- gsub(pattern = "CLOCKMODEL_PLACEHOLDER_model",
                    replacement = "<branchRateModel id=\"RelaxedClock\" spec=\"beast.evolution.branchratemodel.UCRelaxedClockModel\" clock.rate=\"@ucldMean\" rates=\"@branchRate\" tree=\"@TheTree\">\n
                                    <LogNormal id=\"FastLogNormalDistributionModel\" S=\"@ucldStdev\" meanInRealSpace=\"true\" name=\"distr\">\n
                                    <parameter id=\"RealParameter.0\" spec=\"parameter.RealParameter\" estimate=\"false\" name=\"M\">1.0</parameter>\n
                                    </LogNormal>\n
                                   </branchRateModel>\n",
                    x = xml_1)

      # operator
      xml_1 <- gsub(pattern = "<!--CLOCKMODEL_PLACEHOLDER_operators-->",
                    replacement = 
                      
                      "<operator id=\"RateValueScaler\" spec=\"ScaleOperator\" parameter=\"@branchRate\" scaleFactor=\"0.75\" weight=\"8.0\"/>\n
                      <operator id=\"RateValueRandomWalker\" spec=\"RealRandomWalkOperator\" parameter=\"@branchRate\" windowSize=\"1.0\" weight=\"8.0\"/>\n
                      <operator id=\"InternalnodesOperator\" spec=\"consoperators.InConstantDistanceOperator\" clockModel=\"@RelaxedClock\"\n
                      twindowSize=\"1.0\"  tree=\"@TheTree\" rates=\"@branchRate\"  weight=\"10.0\"/>\n
                      <operator id=\"FastRootOperator1\" spec=\"consoperators.SimpleDistance\" clockModel=\"@RelaxedClock\" rates=\"@branchRate\" tree=\"@TheTree\" twindowSize=\"1.0\" weight=\"1.0\"/>\n
                      <operator id=\"FastRootOperator2\" spec=\"consoperators.SmallPulley\" clockModel=\"@RelaxedClock\" dwindowSize=\"1.0\" rates=\"@branchRate\" tree=\"@TheTree\" weight=\"1.0\"/>\n
                      \n
                      <operator id=\"UcldStdevScaler\" spec=\"consoperators.UcldScalerOperator\"\n
                      rates=\"@branchRate\" stdev=\"@ucldStdev\" distr=\"@LogNormal.rateValues\" scaleFactor=\"0.5\" weight=\"3.0\"/>\n
                      \n
                      <operator id=\"ucldMeancaler\" spec=\"ScaleOperator\" parameter=\"@ucldMean\" scaleFactor=\"0.75\" weight=\"20.0\"/>\n
                      <operator id=\"relaxedUpDownOperator.c:anolis\" spec=\"UpDownOperator\" scaleFactor=\"0.75\" weight=\"20.0\">\n
                      <up idref=\"ucldMean\"/>\n
                      <down idref=\"TheTree\"/>\n
                      </operator>\n",
                    x = xml_1)
      
      # logger
        xml_1 <- gsub(pattern = "<!--CLOCKMODEL_PLACEHOLDER_loggers-->",
                      replacement = 
                        "<log id=\"ratesStat\" spec=\"beast.evolution.branchratemodel.RateStatistic\" branchratemodel=\"@RelaxedClock\" tree=\"@TheTree\"/>\n
                        <log idref=\"ucldMean\"/>\n
                        <log idref=\"ucldStdev\"/>\n",
                      x = xml_1)
      
        
        # CLOCKMODEL_PLACEHOLDER_outputTreeLog
        xml_1 <- gsub(pattern = "CLOCKMODEL_PLACEHOLDER_outputTreeLog",
                      replacement = "<log id=\"TreeWithMetaDataLogger\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@TheTree\" branchratemodel=\"@RelaxedClock\" sort=\"false\"/>",
                      x = xml_1)
        
        
        
        if(substitution_tree == TRUE){
          # TREELOG_SUBSTITUTION_PLACEHOLDER
          xml_1 <- gsub(pattern = "<!--TREELOG_SUBSTITUTION_PLACEHOLDER",
                        replacement = "",
                        x = xml_1)
          xml_1 <- gsub(pattern = "TREELOG_SUBSTITUTION_PLACEHOLDER-->",
                        replacement = "",
                        x = xml_1)
          xml_1 <- gsub(pattern = "CLOCKMODEL_PLACEHOLDER_substiOutputTreeLog",
                        replacement = "<log id=\"TreeWithMetaDataLoggerSubstitution\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" substitutions=\"true\" tree=\"@TheTree\" branchratemodel=\"@RelaxedClock\" sort=\"false\"/>",
                        x = xml_1)
        }
        
      
    }
    
    
    
    
    ###############################
    # BDS_EXPONENTIAL_MEAN_PLACEHOLDER BDS_ExponentialMean
    ###############################
    xml_1 <- 
      gsub(pattern = "BDS_EXPONENTIAL_MEAN_PLACEHOLDER",
           replacement = BDS_ExponentialMean,
           x = xml_1)
    
    
    # output names will be adjusted 
    blank_file_name <- paste0(blank_file_name, "_BDSExp", BDS_ExponentialMean)
    
    
    
    ###############################
    # BDMM-prime Skyline model setup
    ###############################
    
    if(skyline_BDMM == T){
      
      if(estimate_changeTimes == TRUE){
        
        blank_file_name <- paste0(blank_file_name, "_ChTmsEst")
        
        xml_1 <- 
          gsub(pattern = "<!--IF_changeTimeEstimate_TRUE",
               replacement = "",
               x = xml_1)
        xml_1 <- 
          gsub(pattern = "IF_changeTimeEstimate_TRUE-->",
               replacement = "",
               x = xml_1)
        
        
      } else if(estimate_changeTimes == FALSE){
        xml_1 <- 
          gsub(pattern = "<!--IF_changeTimeEstimate_FALSE",
               replacement = "",
               x = xml_1)
        xml_1 <- 
          gsub(pattern = "IF_changeTimeEstimate_FALSE-->",
               replacement = "",
               x = xml_1)
      }
      
      
      
      # N_TIMEBINS_PLACEHOLDER
      xml_1 <- 
        gsub(pattern = "N_TIMEBINS_PLACEHOLDER",
             replacement = timebins,
             x = xml_1)
      
      # output names will be adjusted 
      blank_file_name <- paste0(blank_file_name, "_tBins", timebins)
                                
                                
      
      # CHANGE_TIMES_PLACEHOLDER of length(timebins-1)
      xml_1 <- 
        gsub(pattern = "CHANGE_TIMES_PLACEHOLDER",
             replacement = paste0(round((changeTimes - age_offset)/age_scaler, digits = 3), collapse = " "),
             x = xml_1)
      
      # BIRTH_TIMEBINS_PLACEHOLDER 1.0 1.0
      xml_1 <- 
        gsub(pattern = "BIRTH_TIMEBINS_PLACEHOLDER",
             replacement = paste0(rep(birthParameter, timebins), collapse = " "),
             x = xml_1)
      
      # DEATH_TIMEBINS_PLACEHOLDER 1.0 1.0
      xml_1 <- 
        gsub(pattern = "DEATH_TIMEBINS_PLACEHOLDER",
             replacement = paste0(rep(deathParameter, timebins), collapse = " "),
             x = xml_1)
      
      # SAMPLING_TIMEBINS_PLACEHOLDER 0.1 0.1
      xml_1 <- 
        gsub(pattern = "SAMPLING_TIMEBINS_PLACEHOLDER",
             replacement = paste0(rep(samplingParameter, timebins), collapse = " "),
             x = xml_1)
      
      # REMOVALPROB_TIMEBINS_PLACEHOLDER 0.0 0.0
      xml_1 <- 
        gsub(pattern = "REMOVALPROB_TIMEBINS_PLACEHOLDER",
             replacement = paste0(rep(removalParameter, timebins), collapse = " "),
             x = xml_1)
      
      
      # if change times are estimate, the number of dimensions has to be 1 less for the changeTimesCanonical parameters than for the non-changetimes canonical parameters.
      if(estimate_changeTimes == TRUE){
        timebins_estChTi <- timebins-1
        xml_1 <- 
          gsub(pattern = "N_TIMEBINS_estChTi_PLACEHOLDER",
               replacement = timebins_estChTi,
               x = xml_1)
        
        xml_1 <- 
          gsub(pattern = "BIRTH_TIMEBINS_estChTi_PLACEHOLDER",
               replacement = paste0(rep(birthParameter, timebins_estChTi), collapse = " "),
               x = xml_1)
        
        xml_1 <- 
          gsub(pattern = "DEATH_TIMEBINS_estChTi_PLACEHOLDER",
               replacement = paste0(rep(deathParameter, timebins_estChTi), collapse = " "),
               x = xml_1)
        
        xml_1 <- 
          gsub(pattern = "SAMPLING_TIMEBINS_estChTi_PLACEHOLDER",
               replacement = paste0(rep(samplingParameter, timebins_estChTi), collapse = " "),
               x = xml_1)
        
        xml_1 <- 
          gsub(pattern = "REMOVALPROB_TIMEBINS_estChTi_PLACEHOLDER",
               replacement = paste0(rep(removalParameter, timebins_estChTi), collapse = " "),
               x = xml_1)
      }
      
      
    }
    
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
    current_output_path<- file.path(getwd(),"2_scripts","BEAST2_contraband",
                                    paste0("TAXA",number_of_taxa,
                                           "_PCs",number_of_pc_axes_used))
    
    current_output_path_output <- file.path(current_output_path, 
                                            "output")
    
    dir.create(path = current_output_path_output,
               recursive = T)
    
    current_folder_name <- paste0("TAXA",number_of_taxa,
                                  "_PCs",number_of_pc_axes_used)
    
    current_analysis_name <- paste0(blank_file_name,
                                    "_TAXA",number_of_taxa,
                                    "_PCs",number_of_pc_axes_used, 
                                    "_nIter", chainlength_in_millions, "m")
    
    if(underPrior == TRUE) {
      current_analysis_name <- paste0(current_analysis_name, "_underPrior")
    }
    if(SteppingStone == TRUE) {
      current_analysis_name <- paste0(current_analysis_name, "_SteppingStone")
    }  
    if(substitution_tree == TRUE) {
      current_analysis_name <- paste0(current_analysis_name, "_substitutionTree")
    }
    
    
    # output names and paths
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.log",
                  replacement = file.path(".", 
                                          current_folder_name, 
                                          "output",
                                          paste0("out_",
                                                 current_analysis_name,
                                                 ".log")
                  ),
                  x = xml_1)
    
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.trees",
                  replacement = file.path(".", 
                                          current_folder_name, 
                                          "output",
                                          paste0("out_",
                                                 current_analysis_name,
                                                 ".trees")
                  ),
                  x = xml_1)
    
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.typed.trees",
                  replacement = file.path(".", 
                                          current_folder_name, 
                                          "output",
                                          paste0("out_",
                                                 current_analysis_name,
                                                 ".typed.trees")
                  ),
                  x = xml_1)
    
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.typed.node.trees",
                  replacement = file.path(".", 
                                          current_folder_name, 
                                          "output",
                                          paste0("out_",
                                                 current_analysis_name,
                                                 ".typed.node.trees")
                  ),
                  x = xml_1)
    
    xml_1 <- gsub(pattern = "BMPruneLikelihood_calval_1.traj",
                  replacement = file.path(".", 
                                          current_folder_name, 
                                          "output",
                                          paste0("out_",
                                                 current_analysis_name,
                                                 ".traj")
                  ),
                  x = xml_1)
    
    
    ###############################
    # SteppingStone
    ###############################
    
    if(SteppingStone == T){
      
      xml_1 <- gsub(pattern = "<!--STEPPINGSTONE_PATHSAMPLING_TRUE_PLACEHOLDER",
                    replacement = "",
                    x = xml_1)
      xml_1 <- gsub(pattern = "STEPPINGSTONE_PATHSAMPLING_TRUE_PLACEHOLDER-->",
                    replacement = "",
                    x = xml_1)
      
      xml_1 <- gsub(pattern = "STEPPINGSTONE_PATHSAMPLING_TRUE_PATH_PLACEHOLDER",
                    replacement = current_output_path_output,
                    x = xml_1)
      
    } else {
      
      xml_1 <- gsub(pattern = "<!--STEPPINGSTONE_PATHSAMPLING_FALSE_PLACEHOLDER",
                    replacement = "",
                    x = xml_1)
      xml_1 <- gsub(pattern = "STEPPINGSTONE_PATHSAMPLING_FALSE_PLACEHOLDER-->",
                    replacement = "",
                    x = xml_1)
      
    }
    
    
    
    
    
    ###############################
    # save script
    readr::write_file(x = xml_1,
                      file = file.path(current_output_path,
                                       paste0("script_",
                                              current_analysis_name,
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
                     replacement = current_analysis_name,
                     x = job_sh_1)
    
    # SCRIPTPATH_PLACEHOLDER
    job_sh_1 <- gsub(pattern = "SCRIPTPATH_PLACEHOLDER",
                     replacement = file.path(".",
                                             current_folder_name,
                                             paste0("script_",
                                                    current_analysis_name,
                                                    ".xml")),
                     x = job_sh_1)
    
    
    
    
    # write job.sh
    readr::write_file(x = job_sh_1,
                      file = file.path(blank_file_path,
                                       paste0("job_",
                                              current_analysis_name,
                                              ".sh")
                      )
    )
    
    
    
    
  }
