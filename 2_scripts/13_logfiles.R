
library(coda)
library(dplyr)
library(magrittr)


burnin_prop <- 0.2 # e.g., 0.5 == 50% burn-in

output_folder_paths <- 
  list.dirs(file.path(
    "2_scripts",
    "new_xmls",
    "BEAST2_contraband"),
    recursive = F,
    full.names = F)

taxa_pc_combination_list <- list()
# pb = txtProgressBar(min = 0, max = length(output_folder_paths), initial = 0) 
for(i in output_folder_paths){
  cat(i, "\n")
  current <- file.path("2_scripts",
                       "new_xmls",
                       "BEAST2_contraband",
                       i)
  
  output_IR1 <-
    list.files(file.path(current,
                         "independent_run_1",
                         "output"),
               pattern = "*.log",
               full.names = T)
  
  # output_IR2 <-
  #   list.files(file.path(current,
  #                        "independent_run_2",
  #                        "output"),
  #              pattern = "*.log",
  #              full.names = T)
  
  
  current_model_list <- list()
  for (ii in output_IR1) {
    
    current_model <- 
      gsub(strsplit(ii, split = "/")[[1]][length(strsplit(ii, split = "/")[[1]])],
           pattern = ".log", replacement = "")
    print(current_model)
    
    logfile <- read.table(ii, header = T) %>%
      dplyr::slice_head(prop = burnin_prop) # % burn-in of slice_head (i.e., the first % rows)
    
    variables_log <- colnames(logfile)
    
    curent_model_variables_log_list <- list()
    for(iii in 2:length(variables_log)){
      curent_model_variables_log_list[[iii]] <- 
        data.frame(TAXA_PC = i,
                   Taxa = strsplit(i, split = "_")[[1]][2],
                   traits = strsplit(i, split = "_")[[1]][4],
                   current_model = current_model,
                   variable = variables_log[iii],
                   mean = mean(coda::as.mcmc(logfile[,variables_log[iii]])),
                   median = median(coda::as.mcmc(logfile[,variables_log[iii]])),
                   HPDinterval = coda::HPDinterval(coda::as.mcmc(logfile[,variables_log[iii]])),
                   ESS = if(length(unique(coda::as.mcmc(logfile[,variables_log[iii]])))>1){
                     coda::effectiveSize(coda::as.mcmc(logfile[,variables_log[iii]]))
                   } else {
                     NA
                   })
      # coda::traceplot(coda::as.mcmc(logfile[,variables_log[iii]]))
    }
    current_model_list[[current_model]] <- 
      do.call(rbind.data.frame, curent_model_variables_log_list)

  } 
  
  taxa_pc_combination_list[[i]] <- 
    do.call(rbind.data.frame, current_model_list)
  
  # setTxtProgressBar(pb,i)
  
}
# close(pb)

all_combinations_results_table <- 
  do.call(rbind.data.frame, taxa_pc_combination_list)

# readr::write_csv(all_combinations_results_table,
#                  file = file.path("3_output",
#                                   "all_combinations_results_table.csv"))


a <-
  all_combinations_results_table %>% 
  # tibble::as.tibble() %>%
  subset(variable %in% c("prior", "likelihood", "posterior")) %>%
  dplyr::filter(ESS < 200) %>% 
  dplyr::select(TAXA_PC, current_model) %>% 
  unique()

sbatch_resume_file <- 
  file.path("2_scripts", "new_xmls", paste0("resume_ESSunder200_", Sys.Date(),".txt"))

sbatch_resume_fun <- 
  function(x){
    cat("sbatch ", x$TAXA_PC, "/independent_run_1/resume_", x$current_model, ".sh", "\n", 
        "sbatch ", x$TAXA_PC, "/independent_run_2/resume_", x$current_model, ".sh", "\n",
        sep = "",
        file = sbatch_resume_file,
        append = T)
  }

file.create(sbatch_resume_file)
for(i in 1:nrow(a)){
  sbatch_resume_fun(a[i,])
}


# Rachel
ESS_check = function(log, var = "posterior") {
  library(coda)
  effectiveSize(as.mcmc(log[var])) >= 200
}

# logf = list.files(path = file.path("2_scripts", "new_xmls", "BEAST2_contraband"), pattern = "*.log", recursive = T, full.names = T)
logf = unique(readLines(file.path("2_scripts", "new_xmls", "resume_ESSunder200_paths.txt")))

cat("total number of log files = ", length(logf), "\n",file = "ESS_log.txt")

for(i in logf){
  
  cat(i, "\n")
  
  if(file.exists(i)){
    log = read.table(i, header = T, stringsAsFactors = F, comment.char = "#")
  } else{
    cat(i, "does not exist!\n")
    next
  }	
  
  #if(!ESS_check(log, var = "prior")) cat(i, "LOW ESS prior!\n", file = "ESS.log", append =TRUE)
  #if(!ESS_check(log, var = "likelihood")) cat(i, "LOW ESS likelihood!\n", file = "ESS.log", append = TRUE)
  #if(!ESS_check(log, var = "posterior")) cat(i, "LOW ESS posterior!\n", file = "ESS.log", append = TRUE)
  
  log <- log[-c(1:(length(log$Sample) * 0.2)),]
  
  if( ((!ESS_check(log, var = "prior")) || (!ESS_check(log, var = "likelihood")) || (!ESS_check(log, var = "posterior"))) ) {
    
    # report convergence issues
    cat(i, "\n", file = "ESS_log2.txt", append = TRUE)
  } else {
    # split string
    j = strsplit(i, "[.]")[[1]][1]

    # move files
    system(paste0("mv ", j, ".* done/"))
    # mkdir -p done/TAXA_87_PCs_9/independent_run_2/output/ && cp TAXA_87_PCs_9/independent_run_2/output/ULNC_FBD_skyline_age_uncertainty.* done/TAXA_87_PCs_9/independent_run_2/output/

  }	
}

readLines("ESS_log2.txt")



###############
taxa_fun <- function(x){as.integer(strsplit(x, split = "_")[[1]][2])}
trait_fun <- function(x){as.integer(strsplit(x, split = "_")[[1]][4])}

##############################
# nCat2
nCat2_rate_results <- 
all_combinations_results_table %>% 
  # dplyr::filter(current_model == "FBD_constant_rates") %>%
  dplyr::filter(variable %in% c("rateValues.1", "rateValues.2")) %>% 
  dplyr::select("TAXA_PC", "variable", "mean", "median", "current_model") %>% 
  dplyr::mutate(current_model = factor(current_model, levels = c("FBD_constant_rates",
                                                                 "FBD_constant_rate_age_uncertainty",
                                                                 "FBD_constant_rate_age_uncertainty_woffset",
                                                                 "FBD_skyline",
                                                                 "FBD_skyline_age_uncertainty")),
                taxa = sapply(TAXA_PC, FUN = taxa_fun),
                traits = sapply(TAXA_PC, FUN = trait_fun)) %>% 
  dplyr::arrange(., taxa, traits) %>% 
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa)),
                traits = factor(traits, levels = unique(traits))) %>% 
  reshape(., idvar = c("TAXA_PC", "current_model"), timevar = "variable", direction = "wide") %>% 
  dplyr::mutate(difference_between_means = mean.rateValues.1 - mean.rateValues.2,
                difference_between_medians = median.rateValues.1 - median.rateValues.2)

nCat2_rate_median_plot <- 
  nCat2_rate_results %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = taxa.rateValues.2, 
                                  y = traits.rateValues.2, 
                                  fill = difference_between_medians),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "medianRate1-medianRate2",
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 2,
                      dir = "v")
nCat2_rate_mean_plot <- 
  nCat2_rate_results %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = taxa.rateValues.2, 
                                  y = traits.rateValues.2, 
                                  fill = difference_between_means),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "meanRate1-meanRate2",
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 2,
                      dir = "v")



##############################
# ULNC
ULNC_rate_results <- 
  all_combinations_results_table %>% 
  # dplyr::filter(current_model == "FBD_constant_rates") %>%
  dplyr::filter(variable %in% c("rate.mean", "rate.variance")) %>% 
  dplyr::select("TAXA_PC", "variable", "mean", "median", "current_model") %>% 
  dplyr::mutate(current_model = factor(current_model, levels = c("ULNC_FBD_constant_rates",
                                                                 "ULNC_FBD_constant_rate_age_uncertainty",
                                                                 "ULNC_FBD_constant_rate_age_uncertainty_woffset",
                                                                 "ULNC_FBD_skyline",
                                                                 "ULNC_FBD_skyline_age_uncertainty")),
                taxa = sapply(TAXA_PC, FUN = taxa_fun),
                traits = sapply(TAXA_PC, FUN = trait_fun)) %>% 
  dplyr::arrange(., taxa, traits) %>% 
  dplyr::mutate(taxa = factor(taxa, levels = unique(taxa)),
                traits = factor(traits, levels = unique(traits)))

ULNC_rate_plot <- 
  ULNC_rate_results %>% 
  tibble::as.tibble() %>% 
  subset(., variable == "rate.mean") %>% 
  ggplot2::ggplot(data = .) +
  ggplot2::geom_tile(ggplot2::aes(x = taxa, 
                                  y = traits, 
                                  fill = log(median)),
                     color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "red") +
  ggplot2::labs(fill = "log(median(rate.mean))",
                x = "Number of taxa",
                y = "Number of traits") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(~current_model,
                      ncol = 2,
                      dir = "v")





cowplot::plot_grid(nCat2_rate_mean_plot,
                   nCat2_rate_median_plot,
                   ULNC_rate_plot,
                   labels = "AUTO",
                   ncol = 1,
                   byrow = T)
