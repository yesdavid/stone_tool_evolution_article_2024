NGRIP <- readr::read_delim(file.path("1_data",
                                     "temperature",
                                     "ngrip_50y_wo_header.txt"), 
                           delim = "\t", 
                           escape_double = FALSE, 
                           trim_ws = TRUE, 
                           skip = 1)
NGRIP$BP1950 <- NGRIP$ss09sea_age_years_BP_2000 - 50

library(ggplot2)

NGRIP_subset <- 
  dplyr::filter(NGRIP, 
                BP1950 < 16000 & BP1950 > 9500)

ngrip_delta_o_18 <-
  ggplot2::ggplot(data = NGRIP_subset,
                  ggplot2::aes(x = BP1950,
                               y = Del_18O_permille)) +
  ggplot2::annotate(xmin = c(16000, 12900), 
                    xmax = c(14600, 11700),
                    ymin = min(NGRIP_subset$Del_18O_permille)-.06, 
                    ymax = max(NGRIP_subset$Del_18O_permille)+.06,
                    "rect", 
                    alpha = 0.2, fill = c("grey60", "grey60")) +
  ggplot2::annotate("text",
                    y = -37.8,
                    x = c(15400, 13650, 12300, 10500),
                    label = c("GS-2","GI-1","GS-1","Holocene"),
                    size = 4) +
  ggplot2::geom_path() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_reverse(limits = c(16000, 
                             9500), 
                  expand = c(0, 0),
                  minor_breaks = seq(to = 16000,
                                     from = 9500,
                                     by = 100),
                  breaks = seq(to = 16000,
                               from = 10000,
                               by = 1000)) +
  theme_bw() +
  labs(y = expression(paste(δ^18, "O", " (‰)")), 
       x = "Age BP (1950)") +
  theme(legend.position = "none",
        panel.grid.minor = element_line(colour="grey95"),
        panel.grid.major = element_line(colour="grey88"))

ngrip_delta_o_18

ggsave(ngrip_delta_o_18,
       filename = file.path("3_output",
                            "ngrip_delta-o-18.png"),
       width = 15, height = 10, units = "cm", device = "png")
ggsave(ngrip_delta_o_18,
       filename = file.path("3_output",
                            "ngrip_delta-o-18.tiff"),
       width = 15, height = 10, units = "cm", device = "tiff")
