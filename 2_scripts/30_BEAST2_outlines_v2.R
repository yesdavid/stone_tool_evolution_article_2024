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
set.seed(12)
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

# remove outliers?



# determine the number of PC axes to use
# select number of PCs which account for XX.X% of variation
minimum_no_of_pcs_outlines_AR <- Momocs::scree_min(outlines_AR_subset_PCA,
                                                   prop = 0.999) 
minimum_no_of_pcs_outlines_AR 

Momocs::scree_plot(outlines_AR_subset_PCA,
                   nax = 1:15)
Momocs::PCcontrib(outlines_AR_subset_PCA,
                  nax = c(1:15))

number_of_pc_axes_used <- minimum_no_of_pcs_outlines_AR # or other arbitrary number
pcs <- outlines_AR_subset_PCA$x[,1:number_of_pc_axes_used]

# subset the age information for the taxa to the selected outline data
taxa_file_subset <-
  subset(taxa_file_raw,
         taxon %in% outlines_AR_subset_PCA$fac$ARTEFACTNAME) %>% 
  dplyr::distinct()


# modify BEAST xml template
# load xml_helper_function
source(file.path("2_scripts",
                 "30_BEAST2_outlines_XMLhelperFunction.R"))

# xml file set up
xml_helper_function(taxa_file = taxa_file_subset, # age has to be in column called "max"
                    number_of_pc_axes_used = number_of_pc_axes_used, # number of axes chosen
                    pcs = pcs, # pca axes
                    root_age = 40000,
                    clockmodel = "nCat", # "strict", "relaxed", "nCat". "relaxed" or "nCat" works so far only for fossil_age_uncertainty = F,fully_extinct = F,skyline_BDMM = F, ### "relaxed" needs BEAST2 package "ORC"
                    fossil_age_uncertainty = F,
                    fully_extinct = F,
                    skyline_BDMM = F,
                        timebins = 4, # this helper function does not work for timebins <2. Has to be adjusted manually.
                        estimate_changeTimes = F, # logical, if false, provide the following parameters:
                        changeTimes = c(14600, 12900, 11700),  # 13006,  #13,006+-9 calBP is the year of the Laacher See eruption (2021) https://www.nature.com/articles/s41586-021-03608-x   # the date(s) when the timebins change; has to be of length(timebins-1); has to be in the same format as the raw dates provided in taxa_file_raw
                        birthParameter = "1.0",
                        deathParameter = "1.0", 
                        samplingParameter = "0.1", 
                        removalParameter = "0.0",
                    substitution_tree = T,
                    BDS_ExponentialMean = "1.0",
                    SteppingStone = F,
                    underPrior = F,
                    printgen = 100000, # print ever _printgen_ iteration; set it to: chainlength_in_millions/printgen = 10000
                    chainlength_in_millions = 1000,
                    walltime_spec = "24:00:00",
                    blank_file_path <- file.path(getwd(), "2_scripts","BEAST2_contraband") # path to folder where the blank .xml files are
)


#############
# geiger::fitContinuous()
mccTre <- treeio::read.beast("/home/au656892/Documents/Doktor/2_projects/stone_tool_evolution_article_2022/2_scripts/BEAST2_contraband/TAXA71_PCs12/output/out_BMPruneLikelihood_relaxedClock_FBDbds_Offset9.754_BDSExp1.0_TAXA71_PCs12_nIter1000m.NoNegHeightmcc.tre")


ggtree(mccTre,
       aes(color=(rate_median)),
       size = 1,
       # mrsd ="11-01-01"
) +
  # geom_range(range='length_0.95_HPD', color='blue', alpha=.53, size=2) +
  # geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
  # geom_tiplab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
  geom_nodelab(aes(x=branch, label=round((rate_median), 4)), vjust=-.5, size=3) +
  geom_tiplab(aes(x=branch, label=round((rate_median), 4)), vjust=-.5, size=3) +
  # geom_nodelab(aes(x=branch, label=ln(rate_median)), vjust=-.5, size=3) +
  # geom_tiplab(aes(x=branch, label=ln(rate_median)), vjust=-.5, size=3) +
  scale_color_continuous(low="blue", high="red",
                         name = "Rate of\nevolution") +
  theme(legend.position=c(.1, .8),
        plot.title = element_text(hjust = 0.5)) +
  geom_rootedge(rootedge = 0.2) +
  theme_tree2()+
  geom_tiplab(align = T,
              linesize=.5,
              offset = 0.75,
              color = "black") +
  # ggtitle("Bayesian MAP timescaled tree w/ offset, ngen = 3 Mio.") +
  geom_tippoint(pch = 16) +
  geom_nodepoint(pch = 15) +
  geom_treescale()

ggplot(data = mccTre@data,
       aes(x = height_median+9.754,
           y = rate_median)) +
  geom_vline(xintercept = 13,  # laacher see volcano
             color = "green", size = 2) +
  geom_label(aes(x = 13, y = 0.0017, label = "Laacher See Eruption")) +
  geom_vline(xintercept = 14.6,  # end of late pleniglacial
             color = "red", size = 2) +
  geom_label(aes(x = 14.60, y = 0.0017, label = "end of late pleniglacial")) +
  geom_vline(xintercept = 12.9,  # end of bølling allerød complex
             color = "red", size = 2) +
  geom_label(aes(x = 12.7, y = 0.0017, label = "end of bølling allerød complex")) +
  geom_vline(xintercept = 11.7,  # end of younger dryas complex
             color = "red", size = 2) +
  geom_label(aes(x = 11.7, y = 0.0017, label = "end of younger dryas complex")) +
  geom_point() +
  xlim(16,10)


outlines_AR_subset_lengthCM <- (log(outlines_AR_subset$fac$width_cm))
names(outlines_AR_subset_lengthCM) <- (outlines_AR_subset$fac$ARTEFACTNAME)

outlines_AR_subset_lengthCM_naomit <- (na.omit(outlines_AR_subset_lengthCM))
outlines_AR_subset_lengthCM_naomit_v <- as.vector(outlines_AR_subset_lengthCM_naomit)
names(outlines_AR_subset_lengthCM_naomit_v) <- names(outlines_AR_subset_lengthCM_naomit)

# library(geiger)
BM <-
  geiger::fitContinuous(ape::drop.tip(mccTre@phylo,
                                      names(which(is.na(outlines_AR_subset$fac$width_cm)))),
                        dat = outlines_AR_subset_lengthCM_naomit_v,
                        # niter = 10000,
                        ncores = 7,
                        model = "BM")
EB <- 
  geiger::fitContinuous(ape::drop.tip(mccTre@phylo,
                                      names(which(is.na(outlines_AR_subset$fac$width_cm)))),
                        dat = outlines_AR_subset_lengthCM_naomit_v,
                        # niter = 10000,
                        ncores = 7,
                        model = "EB")
# OU doesnt work
OU <-
  geiger::fitContinuous(ape::drop.tip(mccTre@phylo,
                                      names(which(is.na(outlines_AR_subset$fac$width_cm)))),
                        dat = outlines_AR_subset_lengthCM_naomit_v,
                        # niter = 10000,
                        ncores = 7,
                        model = "OU")
aic.vals<-setNames(c(BM$opt$aicc,
                     OU$opt$aicc,
                     EB$opt$aicc),
                   c("BM",
                     "OU",
                     "EB"))
aic.vals

#############
outlines_AR_subset_region <- factor(outlines_AR_subset$fac$Region)
names(outlines_AR_subset_region) <- outlines_AR_subset$fac$ARTEFACTNAME

fitDiscrete(mccTre@phylo, outlines_AR_subset_region)
fitDiscrete(mccTre@phylo, outlines_AR_subset_region, lambda=TRUE)
fitDiscrete(mccTre@phylo, outlines_AR_subset_region, delta=TRUE)
fitDiscrete(mccTre@phylo, outlines_AR_subset_region, kappa=TRUE)

fitDiscrete(mccTre@phylo, outlines_AR_subset_region, linearchange=TRUE)
fitDiscrete(mccTre@phylo, outlines_AR_subset_region, exponentialchange=TRUE)
fitDiscrete(mccTre@phylo, outlines_AR_subset_region, tworate=TRUE)

#############
# # skyline plot
# devtools::install_github("laduplessis/bdskytools")
# library(bdskytools)
# fname <- "/home/au656892/Documents/Doktor/2_projects/stone_tool_evolution_article_2022/2_scripts/BEAST2_contraband/TAXA71_PCs12/output/out_BMPruneLikelihood_relaxedClock_FBDbds_BDMMprime_Offset9.754_BDSExp1.0_tBins4_TAXA71_PCs12_nIter1000m.log"   
# lf    <- bdskytools::readLogfile(fname, burnin=0.1)
# 
# Re_sky    <- bdskytools::getSkylineSubset(lf, "reproductiveNumber")
# Re_hpd    <- bdskytools::getMatrixHPD(Re_sky)
# delta_hpd <- bdskytools::getHPD(lf$becomeUninfectiousRate)
# 
# bdskytools::plotSkyline(1:10, Re_hpd, type='step', ylab="R")
#############

outlines_AR_with_dates$fac %>% 
  select(ARTEFACTNAME, length_cm, width_cm, area_cm2, perimeter_cm, calliper_cm) %>% 
  dplyr::left_join(., taxa_file_raw,
                   by = c("ARTEFACTNAME" = "taxon")) %>% 
  ggplot(aes(x = max, 
             y = length_cm*-9)) +
  # geom_point() +
  scale_x_reverse() +
  theme_bw() +
  
  geom_vline(xintercept = c(15000,11000)) +
  geom_smooth(method = "loess",
              span = 0.4) +
  xlab("years calBP") +
  # ylab("Artefact length (cm)") +
  geom_line(data = NGRIP_subset,
            aes(x = BP1950, y = degree_celsius),
            color = "black") +
  # geom_line(data = GRIP_subset,
  #           aes(x = BP1950, y = degree_celsius),
  #           color = "green") +
  # geom_line(data = GISP2_subset,
  #           aes(x = BP1950, y = degree_celsius),
  #           color = "blue") +
  geom_vline(xintercept = 13006,  # laacher see volcano
             color = "green", size = 2) + 
  geom_label(aes(x = 13156, y = -11, label = "Laacher See Eruption")) +
  
  geom_vline(xintercept = 14600,  # end of late pleniglacial
             color = "red", size = 2) +
  geom_label(aes(x = 14600, y = -12, label = "end of late pleniglacial")) +
  
  geom_vline(xintercept = 12900,  # end of bølling allerød complex
             color = "red", size = 2) + 
  geom_label(aes(x = 12700, y = -13, label = "end of bølling allerød complex")) +
  
  geom_vline(xintercept = 11700,  # end of younger dryas complex
             color = "red", size = 2) +
  geom_label(aes(x = 11700, y = -14, label = "end of younger dryas complex")) +
  ggtitle("NGRIP") + 
  scale_y_continuous(
    "NGRIP Temperature (C)", 
    sec.axis = sec_axis(~ . /-9, name = "artefact length (cm)")
  )



outlines_to_plot_their_metrics <- 
outlines_AR_with_dates$fac %>% 
  select(ARTEFACTNAME, length_cm, width_cm, area_cm2, perimeter_cm, calliper_cm,
         Lat, Long) %>% 
  dplyr::left_join(., taxa_file_raw,
                   by = c("ARTEFACTNAME" = "taxon")) 

whole_dataset_metrics_over_time <- 
  ggplot(outlines_to_plot_their_metrics,
         aes(x = max, 
             y = length_cm)) +
  geom_smooth(span = 0.4) +
  scale_x_reverse() +
  theme_bw() 

metrics_zone_a <- 
  ggplot(subset(outlines_to_plot_their_metrics, Lat >= 55),
         aes(x = max, 
             y = length_cm)) +
  geom_point() +
  geom_smooth() +
  xlim(16000,9000) +
  ylim(0,10) +
  theme_bw()

metrics_zone_b <- 
  ggplot(subset(outlines_to_plot_their_metrics, Lat <= 55 & Lat >= 47.5),
         aes(x = max, 
             y = length_cm)) +
  geom_point() +
  geom_smooth(span = 0.4) +
  xlim(16000,9000) +
  ylim(0,10) +
  theme_bw()

metrics_zone_c <- 
  ggplot(subset(outlines_to_plot_their_metrics,  Lat <= 47.5),
         aes(x = max, 
             y = length_cm)) +
  geom_point() +
  geom_smooth(span = 0.4) +
  xlim(16000,9000) +
  ylim(0,10) +
  theme_bw()



