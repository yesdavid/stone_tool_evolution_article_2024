# !!!! in lines 246 onwards: make sure the absolute paths match the ones on your machine!



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


logf <- list.files(path = file.path("2_scripts", "new_xmls", "BEAST2_contraband"), pattern = "\\.log", recursive = T, full.names = T)

cat("total number of log files = ", length(logf), "\n", file = paste0("ESS_log_", Sys.Date(),".txt"))
cat("#!/bin/bash -l\n\n",
    file = paste0("ESS_logPathsResume_", Sys.Date(),".txt"))
cat("", file = "mv_finished.txt")

for(i in logf){
  
  cat(i, "\n")
  
  log = read.table(i, header = T, stringsAsFactors = F, comment.char = "#")
  
  log <- log[-c(1:(length(log$Sample) * 0.2)),]
  
  check_for_error <- 
    tryCatch(((!ESS_check(log, var = "prior")) || (!ESS_check(log, var = "likelihood")) || (!ESS_check(log, var = "posterior"))), 
             error=function(e) "Error")
  
  if(check_for_error == "Error"){
    cat("Error: ", i, "/\n\n", 
        file = paste0("logfile_error",
                      ".txt"),
        append = T)
  } else if( ((!ESS_check(log, var = "prior")) || (!ESS_check(log, var = "likelihood")) || (!ESS_check(log, var = "posterior"))) ) {
    
    # resume!
    cat(gsub(x = 
               gsub(x = 
                      gsub(x = 
                             gsub(x = i,
                                  # pattern = "\\./",
                                  pattern = "2_scripts/new_xmls/BEAST2_contraband/",
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


#####################################
# check HPC *.o* files for "fatal exceptions"
#####################################

df <- 
  list.files(file.path("2_scripts", "new_xmls", "BEAST2_contraband"),
             pattern = "*\\.o") %>% 
  lapply(., function(x){
    split <- 
      strsplit(x,
               split = ":")[[1]]
    data.frame(runResume = split[1],
               TaxaPc = split[2],
               version = strsplit(split[3],
                                  split = "\\.o")[[1]][1],
               O = strsplit(split[3],
                            split = "\\.o")[[1]][2],
               filename = x)
  }
  ) %>% 
  do.call(rbind.data.frame, .)

# View(df)

latest_o_files <- 
  df %>% 
  group_by(runResume,TaxaPc,version) %>%
  slice(which.max(O)) %>% 
  mutate(TaxaPcVersion = paste0(as.character(TaxaPc), "_", as.character(version)))

latest_o_files 

a_TaxaPcVersion <- 
  data.frame(TaxaPc = sapply(chains_under_ESS200$TAXA_PC,
                             function(x){
                               paste0(strsplit(x,
                                               split = "_")[[1]][c(2,4)],
                                      collapse = "")
                             }),
             version = chains_under_ESS200$current_model,
             TaxaPcVersion = paste0(sapply(chains_under_ESS200$TAXA_PC,
                                           function(x){
                                             paste0(strsplit(x,
                                                             split = "_")[[1]][c(2,4)],
                                                    collapse = "")
                                           }),
                                    "_",
                                    chains_under_ESS200$current_model)
  )

latest_o_files[which((latest_o_files$TaxaPcVersion %in% a_TaxaPcVersion$TaxaPcVersion)),]

fatal_exception_df_list <- list()
for(i in latest_o_files$filename){
  current_o_file <- read.delim(file.path(file.path("2_scripts", "new_xmls", "BEAST2_contraband", i)))
  
  if(identical(grep("Fatal exception", current_o_file), integer(0))){ # if there's no "Fatal exception" in the file, grep returns "integer(0)"
    fatal_exception = F
  } else {
    fatal_exception = T
  }
  fatal_exception_df_list[[i]] <- 
    data.frame(filename = i,
               fatal_exception = fatal_exception)
  
}
fatal_exception_df <- 
  do.call(rbind.data.frame, fatal_exception_df_list)

latest_o_files_w_fatalExc <- 
  dplyr::left_join(latest_o_files,
                   fatal_exception_df,
                   by = "filename") 

summary(latest_o_files_w_fatalExc$fatal_exception)

subset(latest_o_files_w_fatalExc, fatal_exception == T) %>% 
  View()

fatal_exception_plot <- 
  ggplot2::ggplot(data = latest_o_files_w_fatalExc) +
  ggplot2::geom_bar(ggplot2::aes(y = TaxaPc,
                                 fill = fatal_exception),
                    position = "dodge")+
  ggplot2::facet_grid(TaxaPc~version,
                      scales = "free") +
  ggplot2::scale_fill_manual(values = c("green", "red")) +
  ggplot2::theme_bw()  +
  ggplot2::scale_x_continuous(breaks = c(0,1,2))

ggplot2::ggsave(fatal_exception_plot,
                filename = file.path("3_output", "fatal_exception_plot.png"),
                width = 50, height = 30, units = "cm", device = "png")


##################################################################################
# logcombiner
## combine chains if both have converged successfully
##################################################################################
done_path <- file.path("2_scripts",
                       "new_xmls",
                       "done")
taxa_pc_folders_converged <-
  list.files(done_path)


for (current_log_or_tre in c("trees", "log")) {
  
  cat("#!bash \n",
      file = paste0("command_to_combine_", current_log_or_tre, ".sh"))
  
  for(current_folder in taxa_pc_folders_converged){
    
    already_combined_files <- 
      data.frame(filename = list.files(file.path(done_path,
                                                 current_folder),
                                       pattern = paste0("COMBINED.", current_log_or_tre))) %>% 
      dplyr::mutate(model_name = sapply(.$filename, function(x){
        strsplit(x, split = paste0("_COMBINED.",current_log_or_tre))[[1]][1]
      }))
      
    if(nrow(already_combined_files) == 10){
      
    } else {
      
      if(nrow(already_combined_files) > 0){
        all_current_files <-
          data.frame(filename = list.files(file.path(done_path,
                                                     current_folder),
                                           recursive = T,
                                           pattern = paste0(".", current_log_or_tre))) %>% 
          dplyr::mutate(model_name = sapply(.$filename, function(x){
            strsplit(x, split = "/")[[1]][3]
          }) %>% 
            lapply(., function(x){
              strsplit(x, split = paste0(".",current_log_or_tre))[[1]][1]
            })
          ) %>% 
          subset(., !is.na(model_name)) %>% 
          subset(., !(model_name %in% already_combined_files$model_name)) %>% 
          dplyr::pull(filename)
      } else {
        all_current_files <-
          list.files(file.path(done_path,
                               current_folder),
                     recursive = T,
                     pattern = paste0(".", current_log_or_tre))
      }
      
      current_file_df_list <- list()
      for(current_file in all_current_files){
        current_file_strsplit <- 
          strsplit(current_file,
                   split = "/")[[1]]
        
        current_file_df_list[[current_file]] <- 
          data.frame(chain = current_file_strsplit[1],
                     model = current_file_strsplit[3],
                     path = current_file)
      }
      
      if(length(current_file_df_list) > 0){
        # all the models where both chains have converged
        tally_converged_chains <- 
          do.call(rbind.data.frame, current_file_df_list) %>% 
          dplyr::group_by(model) %>% 
          dplyr::add_tally() %>% 
          subset(n == 2)
        
        for(current_converged in unique(tally_converged_chains$model)){
          
          command_to_combine_files <- 
            paste0("/home/Downloads/BEAST.v2.6.2.Linux/beast/bin/logcombiner", 
                   " -log /home/Documents/stone_tool_evolution_article_2022/2_scripts/new_xmls/done/", current_folder, "/independent_run_1/output/", current_converged,
                   " -log /home/Documents/stone_tool_evolution_article_2022/2_scripts/new_xmls/done/", current_folder, "/independent_run_2/output/", current_converged,
                   " -o /home/Documents/stone_tool_evolution_article_2022/2_scripts/new_xmls/done/", current_folder, "/", strsplit(current_converged,
                                                                                                                                                              split = "\\.")[[1]][1],
                   "_COMBINED.",current_log_or_tre,
                   " -b 20",
                   "\n\n")
          
          # system(command_to_combine_logs)
          cat(command_to_combine_files,
              file = paste0("command_to_combine_", current_log_or_tre, ".sh"), 
              append = T)
        }
      }
      
      
    }  
      
  }
  
  
}






