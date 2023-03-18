library(Momocs)
library(splitstackshape)
library(readr)
library(magrittr)

# load outlines + PCA data
outlines_centered_scaled_subset_PCA <- readRDS(file = file.path("1_data",
                                                                "Outlines",
                                                                "TAXA35_PCs35_outlines_centered_scaled_subset_seed1_PCA.RDS"))

# load taxa file
# taxa_file_raw <- readr::read_tsv(file.path("1_data", "outlines_centered_scaled_subset_FAD_LAD_C14.tsv"))
taxa_file_raw <- readr::read_tsv(file.path("1_data","outlines_centered_scaled_subset_FAD_LAD_C14_oneSigmaMinMax.tsv"))

ggplot2::ggplot(data = taxa_file_raw, 
                aes(y = reorder(taxon, -max), 
                    x=max,
                    xmin = oneSigma_rangeMax, 
                    xmax = oneSigma_rangeMin)) + 
  ggplot2::geom_pointrange() 

# # determine the number of PC axes to use
# # select number of PCs which account for XX.X% of variation
minimum_no_of_pcs_outlines_AR <- Momocs::scree_min(outlines_centered_scaled_subset_PCA,
                                                   prop = 0.999)
minimum_no_of_pcs_outlines_AR

Momocs::scree_plot(outlines_centered_scaled_subset_PCA,
                   nax = 1:15)
Momocs::PCcontrib(outlines_centered_scaled_subset_PCA,
                  nax = c(1:7),
                  sd.r = c(-3, -2, -1, 0, 1, 2, 3))

# A) number of axes that account for 99.9% of the total variation
# number_of_pc_axes_used <- minimum_no_of_pcs_outlines_AR # or other arbitrary number
# B) all axes available
number_of_pc_axes_used <- ncol(outlines_centered_scaled_subset_PCA$x) 
# C) or any other arbitrary number
# number_of_pc_axes_used <- 

pcs <- outlines_centered_scaled_subset_PCA$x[,1:number_of_pc_axes_used]



# modify BEAST xml template
# load xml_helper_function
source(file.path("2_scripts",
                 "30_BEAST2_outlines_XMLhelperFunction.R"))

# xml file set up
# OBS! rho = 0, or rho = 0.00000001?
xml_helper_function(taxa_file = taxa_file_raw, # age has to be in column called "max"
                    number_of_pc_axes_used = number_of_pc_axes_used, # number of axes chosen
                    pcs = pcs, # pca axes
                    root_age = 40000,
                    clockmodel = "nCat", # "strict", "relaxed", "nCat", "nCat3", nCat4. "relaxed" or "nCat" works so far only for fossil_age_uncertainty = F,fully_extinct = F,skyline_BDMM = F, ### "relaxed" needs BEAST2 package "ORC"
                    fossil_age_uncertainty = T,
                    fully_extinct = F,
                    skyline_BDMM = T,
                        timebins = 4, # this helper function does not work for timebins <2. Has to be adjusted manually.
                        estimate_changeTimes = T, # logical, if false, provide the following parameters:
                        changeTimes = c(14600, 12900, 11700),  # 13006,  #13,006+-9 calBP is the year of the Laacher See eruption (2021) https://www.nature.com/articles/s41586-021-03608-x   # the date(s) when the timebins change; has to be of length(timebins-1); has to be in the same format as the raw dates provided in taxa_file_raw
                        birthParameter = "1.0",
                        deathParameter = "1.0", 
                        samplingParameter = "0.1", 
                        removalParameter = "0.0",
                    substitution_tree = F,
                    BDS_ExponentialMean = "1.0",
                    SteppingStone = F,
                    underPrior = F,
                    printgen = 20000, # print ever _printgen_ iteration; set it to: chainlength_in_millions/printgen = 10000
                    chainlength_in_millions = 200,
                    walltime_spec = "23:55:00",
                    blank_file_path <- file.path(getwd(), "2_scripts","BEAST2_contraband") # path to folder where the blank .xml files are
)


############################################################
############################################################
# for pathSampling/steppingStone:
############################################################
############################################################
# 1. take the above created *_SteppingStone.xml and cut&paste it into a dedicated pathSampling folder, i.e. BEAST2_contraband/TAXA42_PCs42/pathSampling_80steps/nCat2_SS/
# 
# 2. run the .xml file, it will create a step0 folder + beast.xml script, and a run.sh script in the output-folder. These files lie where contraband.sh lies. cut&paste them into the dedicated pathSampling folder
# 
# 3. at the bottom of the beast.xml file created in folder step0, remove the following lines
# <logger id="Logger" spec="Logger" fileName="likelihood.log" logEvery="50000">
#   
#   <log idref="likelihood"/>
#   
# </logger>
# and the lines for the trees/ .tree file. (not needed for this)

# 4. modify the output path for the log file to this:
# fileName="./TAXA42_PCs42/pathSampling_100steps/nCat2_SS/step0/likelihood.log", where "nCat2_SS" is the folder name of the current run, specified in step 1, and "pathSampling_100steps" the folder from step 1.

# 5. insert the beta argument like so:
# <run id="PathSamplingStep" spec="beast.inference.PathSamplingStep" chainLength="50000000" beta="1" preBurnin="1000000">

############################################################
# 6. run this script:
############################################################
nSteps = 100
taxa_and_pcs <- "TAXA42_PCs42"
model_folder <- "nCat2_SS"
path_to_folder = file.path(getwd(), "2_scripts","BEAST2_contraband",taxa_and_pcs, paste0("pathSampling_", nSteps,"steps"), model_folder)
############################################################

# read step0 beast.xml
xml <- readr::read_file(file.path(path_to_folder,
                                  "step0",
                                  "beast.xml"))
# read step0 run.sh
run.sh <- readr::read_file(file.path(path_to_folder,
                                     "run0.sh"))
############################################################

for(i in nSteps:1){
  
  if(i > 1){
    current_beta_value <- toupper(qbeta(i/nSteps, 0.3, 1.0))
  } else if (i == 1){
    current_beta_value <- 0
  }
  
  current_step <- paste0("step",nSteps-i)
  
  current_out_path <- file.path(path_to_folder,
                                current_step)
  dir.create(path = current_out_path,
             recursive = T)
  
  # print(paste(current_step, current_beta_value))
  
  xml_current <- 
    gsub(pattern = "beta=\"1\"",
         replacement = paste0("beta=\"",current_beta_value,"\""),
         x = xml)
  
  xml_current <- 
    gsub(pattern = "step0",
         replacement = current_step,
         x = xml_current)
  
  run.sh_current <- 
    gsub(pattern = "step0",
         replacement = current_step,
         x = run.sh)
  
  
  
  # write xml
  readr::write_file(x = xml_current,
                    file = file.path(current_out_path,
                                     "beast.xml"))
  # write sh
  readr::write_file(x = run.sh_current,
                    file = file.path(path_to_folder,
                                     paste0(model_folder,"_run", nSteps-i, ".sh")))
  Sys.sleep(0.1)
}
############################################################

# 1. use 48-128 steps
# 1.1. compare stepping-stone and pathsampling in revbayes. If output the same, then enough steps (https://revbayes.github.io/tutorials/model_selection_bayes_factors/bf_intro.html). 
# You can use the .log files from the BEAST2 - applauncher:pathSampleAnalyser 
#   (however, the theta column has to be renamed to "power" and all columns other than "Step","power","likelihood" have to be deleted). 
# For the same model, the values of ps.marginal() and ss.marginal() have to be _quite_ close (<= 0.1 range).
# 
# For nCat2 clock model:
#   ss = steppingStoneSampler(file="output_nCat2_pathSamplingAnalyser.out", powerColumnName="power", likelihoodColumnName="likelihood")
# ss.marginal()
# ps = pathSampler(file="output_nCat2_pathSamplingAnalyser.out", powerColumnName="power", likelihoodColumnName="likelihood")
# ps.marginal()
# 
# For strickt clock model:
#   ss = steppingStoneSampler(file="output_strictClock_pathSamplingAnalyser.out", powerColumnName="power", likelihoodColumnName="likelihood")
# ss.marginal()
# ps = pathSampler(file="output_strictClock_pathSamplingAnalyser.out", powerColumnName="power", likelihoodColumnName="likelihood")
# ps.marginal()
# 
# 2. high burn-in
# 3. lower chain length (1/10th of current)

############################################################
############################################################
