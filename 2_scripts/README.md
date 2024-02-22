# Contents

`./new_xmls` contains all BEAST2 scripts and the output MCC trees from this study.

`contraband.jar` is the Java .jar-file which has to be used to run the xml scripts using the contraband package for Brownian Motion (see `./new_xmls/run_BLANK.sh`).

`10_C14_dates_v3.R`: prepare and calibrate C14 data, plot.
	
`11_subset_outlines.R`: load outlines, combine with C14 data, harmonics calibration, PCA.
	
`11_2_plot_subset_outlines.R`: plot outlines and key sites.

(`12_00_FUNCTION_plot_different_age_scalings.R`: optional function for `12_fill_BEAST2_scripts.R` to visualize the impact of different ways to scale the age information.)

`12_fill_BEAST2_scripts.R`: script to fill the BEAST2-xml files with data.
	
`13_LOCAL_move_converged_runs.R`: script to check whether prior, likelihood, and posterior have converged/ESS>>200. If yes, moves scripts and output to `./new_xmls/done` folder, if no, creates "resume.sh" file to resume job on HPC. 

`14_logfiles.R`: script to load .log-files of converged runs and summarize evolutionary rates and variances.
	
`15_plot_Skyline.R`: plots skyline results for NCAT2 clock model.
	
`15_plot_Skyline_ULNC.R`: plots skyline results for ULNC model.
	
`15_plot_Skyline_underPrior.R`: plots skyline results under the prior.
	
`15_plot_Tree.R`: script to create MCC trees of converged runs, plot different summaries.

`delta_o_18.R`: script to prepare and plot the Delta 18O data from the Supplementary Data to _NGRIP Members, Nature 431, 147 - 151 (09 September 2004); doi:10.1038/nature02805_.
