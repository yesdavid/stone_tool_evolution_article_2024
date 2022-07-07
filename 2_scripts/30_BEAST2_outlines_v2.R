library(Momocs)
library(splitstackshape)
library(readr)
library(magrittr)

outlines_AR_centered_w_metric_measures_raw <- readRDS(file = file.path("..",
                                                                       "CLIOARCH_workshop_2_2020",
                                                                       "post_workshop_2021", 
                                                                       "3_output", 
                                                                       "outlines_AR_centered_w_metric_measures.RDS"))

# load taxa file
taxa_file_raw <- readr::read_tsv(file.path("..",
                                           "CLIOARCH_workshop_2_2020", 
                                           "post_workshop_2021", 
                                           "001_data_14c", 
                                           "outlines_AR_taxa_14c_9k16k_v2.tsv"))

# subset outlines 
## select outlines for which 14C date is available
outlines_AR_with_dates <-
  Momocs::filter(outlines_AR_centered_w_metric_measures_raw, ARTEFACTNAME %in% 
                   taxa_file_raw$taxon)

## stratified sub-sample
set.seed(1234)
stratified_subset_of_names <- 
  splitstackshape::stratified(outlines_AR_with_dates$fac,
                              group = c("Region", "Timeslice_range_from", "TaxUnit"),
                              size = 1)$ARTEFACTNAME
outlines_AR_subset <- 
  Momocs::filter(outlines_AR_with_dates, ARTEFACTNAME %in% stratified_subset_of_names)

### check shapes of stratified sub-sample
Momocs::panel(outlines_AR_subset,
              fac = "TaxUnit")

## pca of selected outlines
current_outlines_centered_w_metric_measures <- outlines_AR_subset

current_outlines_centered_w_metric_measures_scaled <- 
  Momocs::coo_scale(current_outlines_centered_w_metric_measures)

### harmonic calibration. Estimates the number of harmonics required for the Fourier methods implemented in Momocs. This is the only step in this section that produces data we need in the subsequent step.
current_outlines_centered_w_metric_measures_scaled_harmonics <- 
  Momocs::calibrate_harmonicpower_efourier(current_outlines_centered_w_metric_measures_scaled, 
                                           plot = F)
# current_outlines_centered_w_metric_measures_scaled_harmonics$minh

### efourier
current_outlines_centered_w_metric_measures_scaled_efourier <- 
  Momocs::efourier(current_outlines_centered_w_metric_measures_scaled,
                   nb.h = as.matrix(current_outlines_centered_w_metric_measures_scaled_harmonics[["minh"]])[[4,1]], # chooses number of harmonics for 99.9%
                   norm = F) 
### PCA
outlines_AR_subset_PCA <- 
  Momocs::PCA(current_outlines_centered_w_metric_measures_scaled_efourier) # PCA on Coe objects, using prcomp.

plot(outlines_AR_subset_PCA$x[,1], outlines_AR_subset_PCA$x[,2])



# determine the number of PC axes to use
# select number of PCs which account for XX.X% of variation
minimum_no_of_pcs_outlines_AR <- Momocs::scree_min(outlines_AR_subset_PCA,
                                                   prop = 0.990) 
minimum_no_of_pcs_outlines_AR 

Momocs::scree_plot(outlines_AR_subset_PCA,
                   nax = 1:15)
Momocs::PCcontrib(outlines_AR_subset_PCA,
                  nax = c(1:15))

number_of_pc_axes_used <- minimum_no_of_pcs_outlines_AR # or other arbitrary number
pcs <- outlines_AR_subset_PCA$x[,1:number_of_pc_axes_used]





# modify BEAST xml template
# load xml_helper_function
source(file.path("2_scripts",
                 "30_BEAST2_outlines_XMLhelperFunction.R"))

# xml file set up
xml_helper_function(fossil_age_uncertainty = F,
                    fully_extinct = F,
                    skyline_BDMM = F,
                        timebins = 2, # this helper function does not work for timebins <2. Has to be adjusted manually.
                        changeTimes, # the date(s) when the timebins change; has to be of length(timebins-1); has to be in the same format as the raw dates provided in taxa_file_raw
                        birthParameter = "1.0",
                        deathParameter = "1.0", 
                        samplingParameter = "0.1", 
                        removalParameter = "0.0",
                    BDS_ExponentialMean = "1.0",
                    underPrior = F,
                    printgen = 100000, # print ever _printgen_ iteration; set it to: chainlength_in_millions/printgen = 10000
                    chainlength_in_millions = 1000,
                    walltime_spec = "24:00:00",
                    blank_file_path <- file.path(getwd(), "2_scripts","BEAST2_contraband") # path to folder where the blank .xml files are
)




