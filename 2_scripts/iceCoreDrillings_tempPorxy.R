# data from http://www.iceandclimate.nbi.ku.dk/data/
# code from https://dominicroye.github.io/en/2018/how-to-create-warming-stripes-in-r/


# Data

# GRIP, Oxygen Isotopes. 20 year averages back to 122 kyrs BP. Chronology is ss09sea or GRIP2001 (Johnsen et al 2001). 
# two measurements for each year?
GRIP <- readr::read_delim(file.path("1_data",
                                    "temperature",
                                    "grip-ss09sea-cl-20yr_wo_header.txt"), 
                          delim = "\t", 
                          escape_double = FALSE, 
                          trim_ws = TRUE)
GRIP$BP1950 <- GRIP$BP1990 - 40

# GISP2, Oxygen Isotopes. 20 years averages on GISP2 time scale, 375 - 103000 yrs BP.
GISP2 <- readr::read_table2(file.path("1_data",
                                      "temperature",
                                      "gisp2_oxygen_20years_wo_header.txt"))
GISP2$BP1950 <- GISP2$`GISP2_timescale[1950]`


# NGRIP, 50 year means of δ18O values from the NGRIP ice core (NGRIP members, 2004) 
# two measurements for each year?
NGRIP <- readr::read_delim(file.path("1_data",
                                     "temperature",
                                     "ngrip_50y_wo_header.txt"), 
                           delim = "\t", 
                           escape_double = FALSE, 
                           trim_ws = TRUE, 
                           skip = 1)
NGRIP$BP1950 <- NGRIP$ss09sea_age_years_BP_2000 - 50



## subset data to Last Glacial-Interglacial Transition (15000-11000 BP)
GRIP_subset <- subset(GRIP, BP1950 >= 11000 & BP1990 <= 15000)
GISP2_subset <- subset(GISP2, BP1950 >= 11000 & BP1950 <= 15000)
NGRIP_subset <- subset(NGRIP, BP1950 >= 9000 & BP1950 <= 16000)


## extrapolation of temperature
# Extrapolation of temperature after @epstein_revised_1953 using the formula $T = 16.5 -4.3* \delta + 0.14*\delta^2$, where $T$ equals the temperature in °C and $\delta$ equals the permil differences between the ration of masses 46 and 44 (i.e. $\delta^{18}O$).

degrees_celcius_extrapol <- function(x){16.5 - (4.3*x + 0.14*x^2)} 

GRIP_subset$degree_celsius <- degrees_celcius_extrapol(GRIP_subset$`del18O permille`)
NGRIP_subset$degree_celsius <- degrees_celcius_extrapol(NGRIP_subset$Del_18O_permille)
GISP2_subset$degree_celsius <- degrees_celcius_extrapol(GISP2_subset$GRIP_delta_O18)

ggplot() + 
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
  scale_x_reverse() +
  theme_bw() +
  ggtitle("NGRIP")


ggplot() + 
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
  scale_x_reverse() +
  theme_bw() +
  ggtitle("NGRIP")
  





