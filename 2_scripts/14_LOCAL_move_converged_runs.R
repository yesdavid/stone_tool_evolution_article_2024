# check whether prior, likelihood, and posterior have converged/ESS>>200. If yes, move to "done", if no, create file to run "resume.sh"

if (!require("coda", character.only=T, quietly=T)) {
  install.packages("coda", repos="http://cran.us.r-project.org")
  library("coda", character.only=T)
} else {
  library("coda", character.only=T)
}


ESS_check = function(log, var = "posterior") {
  effectiveSize(as.mcmc(log[var])) >= 200
}


# setwd("/home/au656892/Documents/Doktor/2_projects/stone_tool_evolution_article_2022/2_scripts/new_xmls/BEAST2_contraband")

logf = list.files(path = file.path("2_scripts", "new_xmls", "BEAST2_contraband"), pattern = "\\.log", recursive = T, full.names = T)
# logf = unique(readLines(file.path("2_scripts", "new_xmls", "resume_ESSunder200_paths.txt")))

cat("total number of log files = ", length(logf), "\n", file = paste0("ESS_log_", Sys.Date(),".txt"))
cat("#!/bin/bash -l\n\n",
    file = paste0("ESS_logPathsResume_", Sys.Date(),".txt"))


for(i in logf){
  
  cat(i, "\n")
  
  log = read.table(i, header = T, stringsAsFactors = F, comment.char = "#")
  
  log <- log[-c(1:(length(log$Sample) * 0.2)),]
  
  if( ((!ESS_check(log, var = "prior")) || (!ESS_check(log, var = "likelihood")) || (!ESS_check(log, var = "posterior"))) ) {
    
    # report convergence issues
    cat(gsub(x = gsub(x = gsub(x = gsub(x = i,
                                        pattern = "\\./",
                                        replacement = "sbatch "),
                               pattern = "\\.log",
                               replacement = ".sh"),
                      pattern = "output/",
                      replacement = "resume_"), 
             pattern = "rates.sh",
             replacement = "rate.sh"), 
        "\n", 
        file = paste0("ESS_logPathsResume_", 
                      Sys.Date(),
                      ".txt"),
        append = T)
    
  } else {
    # split string
    j = strsplit(i, ".log")[[1]][1] # output files only
    
    # move output files
    system(paste0("mkdir -p ./2_scripts/new_xmls/done/", paste0( strsplit(i, "/")[[1]][3:6], collapse = "/"), "/ && mv ", j, ".* ./2_scripts/new_xmls/done/", paste0( strsplit(i, "/")[[1]][3:6], collapse = "/"), "/"))
    cat(paste0("mkdir -p ./2_scripts/new_xmls/done/", paste0( strsplit(i, "/")[[1]][3:6], collapse = "/"), "/ && mv ", j, ".* ./2_scripts/new_xmls/done/", paste0( strsplit(i, "/")[[1]][3:6], collapse = "/"), "/\n\n"),
        file = "mv_finished.txt", append = T)
  }	
}

