#
# Liam Lachs, Newcastle University, 12/2023 
#
# Analysis of outputs from ReefMod-PALAU_v1.1 to investigate whether thermally sensitive
# corymbose Acropora populations in Palau will be able to adapt fast enough to keep pace 
# with ocean warming. This study explores different parameterisations of inheritance of
# coral heat tolerance (from h2=0 to h2=1) and different future scenarios of selective 
# pressure from 16 global climate models (with a range of climate sensitivities) and 
# across 3 different future emissions scenarios (shared socioeconomic pathways): 
# SSP1-2.6: Paris Agreement scenario with global warming limited to 2C
# SSP2-4.5: Miodle-of-the-road scenario
# SSP5-8.5: Worst-case scenario

rm(list = ls())

library(R.matlab)
library(lubridate)
library(cowplot)
library(caTools)
library(ggtext)
library(tidyr)
library(emmeans)
library(dplyr)
library(grid)
library(ggpubr)
library(ggplot2)
library(assertthat)
library(purrr)
library(igraph)
library(sf)
library(ggraph)
library(ggmap)
library(swfscMisc)
library(ggpmisc)
library(rnaturalearth)
library(readxl)
library(rnaturalearthdata)
library(ncdf4)

doingTTacclim = F
if(doingTTacclim){TTacclim = "_TTacclim"} else {TTacclim=""}


setwd("C:/Users/dlcyli/OneDrive/PhD_Code/SLiM/Lachs et al 2024/Figshare data and code/REEFMOD-PALAU_v1.1_analysis")

sens_mcs=1 # CHANGE HERE ONLY

# OISST data from UCAR
OISST_nc = nc_open("C:/Users/dlcyli/OneDrive/PhD_Code/SLiM/Lachs et al 2024/oisst_of_coral_data_lachs.nc")
oisst_lat = ncvar_get(OISST_nc, "matched_lat")
oisst_lon = ncvar_get(OISST_nc, "matched_lon")
sst = ncvar_get(OISST_nc, "sst")
# time is wrongly saved. Let's get it back
start_date <- as.Date("1982-01-01")
end_date <- as.Date("2023-12-31")
oisst_dates <- seq(from = start_date, to = end_date, by = "days")
all_coords = data.frame(oisst_lon, oisst_lat)
colnames(all_coords)[1] = "lon"
colnames(all_coords)[2] = "lat"
unique_coordinates = unique(all_coords[c("lon","lat")])

### get max month and MMM of each unique OISST coordinate
unique_coordinates$max_month = NA
unique_coordinates$MMM = NA

for (i in 1:length(unique_coordinates$lon)) {
  # what is the index in all_coords?
  chosen_index = which(all_coords$lon == unique_coordinates$lon[i] &  all_coords$lat == unique_coordinates$lat[i])[1]
  
  # double sst[time,reefcheck_index]
  sst_of_interest = sst[,chosen_index]
  sst_of_interest = data.frame(sst = sst_of_interest, 
                               time = oisst_dates)
  
  # filter ts between 1985 and 2012 (this is from the original MMM definition)
  start <- as.Date("1985/1/1")
  end <- as.Date("2012/12/31")
  sst_of_interest <- sst_of_interest[!(sst_of_interest$time > end | sst_of_interest$time < start), ]
  
  
  # Calculate mean temperature for each month in each year
  monthly_data <- aggregate(sst ~ format(time, format = "%Y-%m"), data = sst_of_interest, mean)
  colnames(monthly_data)[1] = "time"
  monthly_data$Year <- as.integer(substr(monthly_data$time, 1, 4))
  monthly_data$Month <- as.integer(substr(monthly_data$time, 6, 7))
  
  # Perform linear regression
  final_clim <- numeric(0)
  for (j in 1:12) {
    monthly_means_filtered <- subset(monthly_data, Month == j)
    lm_result <- lm(sst ~ Year, data = monthly_means_filtered)
    predicted_mean_temp <- coef(lm_result)[1] + coef(lm_result)[2] * 1988.286
    final_clim <- c(final_clim, predicted_mean_temp)
  }
  
  unique_coordinates$MMM[i] = max(final_clim)
  unique_coordinates$max_month[i] = which.max(final_clim)
}

# got 29.04901, which is close to the value of 29.2deg obtained in the Humanes  et al. paper

# ---------------------------------------------------- #  -----------------------
#                   Figure 10 - Heat tolerance        (FINAL)    #  -----------------------
# ---------------------------------------------------- #
BMIfun = function(x,Status_levs){
  props = c()
  for(i in 1:length(Status_levs)){
    props[i] = length(which(x == Status_levs[i]))/length(x)
  }
  BMI = sum(props * Status_levs)/(length(Status_levs)-1)
  return(BMI)
}
PredictorProb50 = function(B0, B1){
  (qlogis(0.5) - B0)/ B1 
}

ST = read_xlsx("long heat stress/Long_heat_stress_05.31.22.xlsx", sheet=1)


colnames(ST)
unique(ST$Status)

Status_levs = c(0,1,2,3,4)
# status means health... healthy, bleached, dead
# healthy (c1), partially bleached (c2), fully bleached (c3), partially dead (c4) or dead (c5)
ST$Status = as.numeric(dplyr::recode(ST$Status,
                                     H=0,P=0,p=0,PB=1,B=2,PM=3,D=4))
ST$Colony = as.factor(ST$Colony)
ST$Nubbin_ID = as.factor(ST$Nubbin_ID)
ST$Tank = as.factor(ST$Tank)
ST$Cohort = as.factor(ifelse(substr(ST$Site,1,1)=="N","F1","Wild"))
ST$Treat = as.factor(ifelse(is.na(match(ST$Tank,c(9,10,13))),"Stress","Control"))

# Remove colonies with too few alive/healthy nubbins in stress tanks. 
# Note if we use a minimum of 2 healthy nubbins on the first day of the experiment, then we remove colony 35 which is a parental colony
# So we then use a minimum of 2 alive  nubbins on the first day of the experiment instead
# Also need a minimum of 1 nub/colony alive in control on day 
# and at least one nub/colony must stay alive for the whole experiment (last day)
ST = ST %>%
  group_by(Colony) %>%
  mutate(NumHealthyNubsDay1Stress = sum(Status == 0 & Treat == "Stress" & Date == min(Date)),
         NumAliveNubsDay1Stress = sum(Status < 4 & Treat == "Stress" & Date == min(Date)),
         NumAliveNubsDay1Control = sum(Status < 4 & Treat == "Control" & Date == min(Date)),
         NumAliveNubsEndControl = sum(Status < 4 & Treat == "Control" & Date == max(Date[Treat=="Control"]))) 
sort(unique(ST$Colony))

# Tally number of colonies lost from experiment due to not enough stress nubbins
Dumped = ST[ST$NumAliveNubsDay1Stress==1,]# & 
Dumped$Colony = as.character(Dumped$Colony)
length(unique(Dumped$Colony))

# Tally number of colonies lost from experiment due to no control nubbin at beginning or end
Dumped = ST[ST$NumAliveNubsDay1Control == 0 | 
              ST$NumAliveNubsEndControl == 0,]
Dumped$Colony = as.character(Dumped$Colony)
length(unique(Dumped$Colony))

ST = ST[ST$NumAliveNubsDay1Stress>=2 &
          ST$NumAliveNubsDay1Control >= 1 &
          ST$NumAliveNubsEndControl >= 1,
        -((ncol(ST)-1):ncol(ST))]

# Check the number of nubbins per coral
checks = ST %>% 
  group_by(Cohort, Colony, Nubbin_ID) %>%
  summarise(N = length(Status)) 
unique(checks$N) # should be four entries
# check N colonies alive per day - to subset end of long stress experiment
checks = ST[ST$Treat=="Stress",] %>% 
  group_by(Date) %>%
  summarise(NcolAlive = length(unique(Colony[Status<4])))
unique(checks$NcolAlive)
ST = ST[!(ST$Date > ymd("2022-05-24")),]
colnames(ST)


# Adding DHW data
DHW = read_excel("long heat stress/Heat Stress Experiment Temperatures.xlsx", sheet = "2022 Long")

# if we are doing 0.1 C per decade (run small commented code below once to generate dataset)
# DHW = read_excel("long heat stress/Heat Stress Experiment Temperatures 0.1Cpd.xlsx", sheet = "2022 Long")
# DHW$DHW=DHW$DHW_0.1pd

# DHW$Date.time = ymd_hms(DHW$Date.time)
DHW$Date = date(DHW$Date.time)
# DHW$Hour = as.character(hour(DHW$Date.time))
# DHW$Minute = as.character(minute(DHW$Date.time))
DHW$Year = as.factor(DHW$Year)

#              ASIDE                  #
# get HOBOdat from Dataset Formation Long script
# H2 = HOBOdat[!is.na(match(colnames(HOBOdat),c("Date.time","Tank","DHW_Nu0.1pd")))]
# H2$Tank = as.numeric(as.character(H2$Tank))
# colnames(H2)[3] = "DHW_0.1pd"
# class(H2$DHW_0.1pd)
# DHW2=left_join(DHW,H2)
# write.csv(DHW2, "long heat stress/Heat Stress Experiment Temperatures 0.1Cpd.csv",row.names = F)
# -------------- END ---------------- #

DHW2 = DHW[(hour(DHW$Date.time)==9 & minute(DHW$Date.time)==0),]
summary(DHW2$Temp.cal)
# Number of tanks per specdified treatment
DHW2 %>% filter(Treat=="Stress") %>% 
  summarise(paste(unique(Tank),collapse="_"))
# It seems that the MMM data used is only from Humanes et al. (2024), whereby: The local climatological
# baseline was adjusted to 29.4 °C for Maschechur reef and 29.6 °C for the
# nursery based on the relationship between satellite sea surface
# temperatures and in situ temperatures21. 

# 13 tanks... (including control) need to find out their locations
# 10 stress tanks
# The final heat stress exposure on the 42nd day was 16 °C-weeks 
# (Supplementary Text) leading to 99% mortality of colony fragments.
# This means that it was a 6-week exposure. However, the papers 
# (Lachs et al 2024) state that it was "a long-term 5-week simulated marine heatwave experiment"
unique(DHW2$Tank)

# given that the MMM adj used are all the same, let us compute MHWs at the same location: 7.375, 134.625
library(heatwaveR)
OISST_of_interest = data.frame(t = oisst_dates, temp = sst[,7])
ts = ts2clm(OISST_of_interest, climatologyPeriod = c("1983-01-01", "2012-12-31"))
for (i in 1:length(DHW2$MMM.adj)) {
  DHW2$MMM.adj[i] = ts$thresh[which(ts$t == DHW2$Date[i])]
}
DHW2$HotSpot = DHW2$Temp.cal - DHW2$MMM.adj
for (i in unique(DHW2$Tank)) {
  temp_hs = DHW2$HotSpot[which(DHW2$Tank == i)]
  dhw_vect = vector()
  for (j in temp_hs) {
    if (length(dhw_vect) == 0) {
      if (j < 0) {
        dhw_vect = append(dhw_vect, 0)
      } else {
        dhw_vect = append(dhw_vect, j/7)
      }
    } else {
      if (j < 0) {
        dhw_vect = append(dhw_vect, dhw_vect[length(dhw_vect)])
      } else {
        dhw_vect = append(dhw_vect, dhw_vect[length(dhw_vect)] + j/7)
      }
    }
  }
  DHW2$DHW[which(DHW2$Tank == i)] = dhw_vect
}

ST = ST %>%
  group_by(Nubbin_ID, Date) %>%
  mutate(DHW = ifelse(length(DHW2$DHW[DHW2$Date == Date & DHW2$Tank == Tank]) > 0,
                      DHW2$DHW[DHW2$Date == Date & DHW2$Tank == Tank],
                      mean(DHW2$DHW[DHW2$Date == Date & DHW2$Treat == Treat])))


# N genets / treatment
ST %>% filter(Treat=="Stress") %>%
  summarise(N = length(unique(Colony)))

# N genets / tank / experiment
ST %>% filter(Treat=="Stress") %>% 
  group_by(Tank) %>%
  summarise(N = length(unique(Colony)),
            range = paste(min(N),"-",max(N)))

# NOTE: there seems to be some package conflict after this step... 
# If BMI has only 1 obs, close R and re-load all packages

# Computing BMIs
BMI = ST %>%
  group_by(Cohort, Region=substr(Site,1,3),Site, Colony, Date) %>%
  summarise(BMI = BMIfun(Status[Treat=="Stress"],Status_levs),
            DHWav = mean(DHW[Treat=="Stress"]))

# raw data plot
ggplot() + 
  geom_line(data = BMI,
            aes(x=DHWav, y=BMI, 
                group=Colony),
            size=1.8, alpha = 0.1)

# we can do the below to compute a "bleaching survival probability", which will help to compute the phenoFitness later on
BMI$BMI = abs(BMI$BMI-1)

# make MCSs more resistant here
BMI$DHWav = BMI$DHWav*sens_mcs

# # we will also add a theoretical Xth_perc colony
# Xth_perc_bmi = BMI %>%
#   group_by(Date) %>%
#   summarise(Xth_perc_bmi = quantile(BMI, 0.9)) # but first, let's compute Xth_perc BMI and Xth_perc DHW for each day
# 
# Xth_perc_dhw = BMI %>%
#   group_by(Date) %>%
#   summarise(Xth_perc_dhw = quantile(DHWav, 0.9))

# BMI = rbind(BMI, BMI[1:30,])
# BMI$Colony = as.character(BMI$Colony)
# BMI$Colony[8341:8370] = "Xth_perc"
# BMI$BMI[which(BMI$Colony == "Xth_perc")] = Xth_perc_bmi$Xth_perc_bmi
# BMI$DHWav[which(BMI$Colony == "Xth_perc")] = Xth_perc_dhw$Xth_perc_dhw

mod = glm(BMI ~ DHWav + Colony, data = BMI, family="binomial") # IMPORTANT NOTE: coefficients were derived using the original eq... don't change the signs cos B1 will automatically be -ve
summary(mod)

# HT contains the intercept and slope for each 278 colony
HT = data.frame(Colony = c(unique(BMI$Colony)[is.na(match(unique(BMI$Colony), sub("Colony","",names(mod$coefficients[-1:-2]))))],
                           sub("Colony","",names(mod$coefficients[-1:-2]))),
                Intercept = c(mod$coefficients[1], mod$coefficients[1] + mod$coefficients[-1:-2]),
                Slope = mod$coefficients[2])
mean(HT$Intercept)               
sd(HT$Intercept)

# compute the DHW at which we have 50% mortality for each colony
HT$DHW50 = PredictorProb50(HT$Intercept,HT$Slope)
hist(HT$DHW50, xlab = "DHW50 (degC-weeks)")
mean(HT$DHW50)
var(HT$DHW50)^0.5
summary(HT$DHW50)
# ---------------- output the DHW50 quantiles here
quantiles = seq(0,1,by=0.01)
dhw50_given_quantile = vector()
for (i in quantiles) {
  dhw50_given_quantile = append(dhw50_given_quantile, quantile(HT$DHW50, i))
}
hist(dhw50_given_quantile)

# ---------------- output the DHW50 quantiles here, BUT using rank() which should approximate behaviour of quantile()
quantiles = seq(0,1,by=0.01)
dhw50_given_quantile = vector()
quantile_rank = function(x,qt){
  # x is the vector of numbers
  ranks <- rank(x)
  normalized_ranks <- (ranks - 1) / (length(x) - 1)
  index <- which.min(abs(normalized_ranks - i))
  return(x[index])
}

for (i in quantiles) {
  dhw50_given_quantile = append(dhw50_given_quantile, quantile_rank(HT$DHW50, i))
}
hist(dhw50_given_quantile)

dhw50_given_quantile = as.data.frame(dhw50_given_quantile)
# write.table(t(dhw50_given_quantile), file = "C:/Users/dlcyli/OneDrive/PhD_Code/SLiM/Main codes/dhw50_given_quantile.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
# test = read.table("C:/Users/dlcyli/OneDrive/PhD_Code/SLiM/Main codes/dhw50_given_quantile.txt")

# HT contains the intercept and slope for all the full distribution of DHW50
HT = data.frame(DHW50 = dhw50_given_quantile,
                Slope = mod$coefficients[2]) # used the model to get the slope and intercept
colnames(HT)[1] = "DHW50"
HT$Intercept = HT$DHW50*HT$Slope*-1

# fake colonies that span the full distribution
HT$Colony = 1:length(HT$DHW50)
HT$Colony = as.factor(HT$Colony)

# let us generate the curve at higher res
# ND = expand.grid(Colony=HT$Colony[order(HT$DHW50)],
#                  DHWav=seq(0,16,by=0.05)) 

# --------------- you want to regenerate the full curves for your fake colonies, using regular intervals of the DHW and computing the probability using the given slope and intercept for each colony
ND = data.frame(Colony = vector(), DHWav = vector(), pred = vector())
PBSpred = function(DHW, B0){
  B1 =mod$coefficients[2]
  pred = exp(B0 + B1*DHW)/(1+exp(B0 + B1*DHW))
}
# takes a few mins to run
dhw_range = 16*sens_mcs
chosen_res = 0.05
# chosen_res = 0.1 # update: 23 Sep 2025 (trying to make the model run faster)
for (i in 1:length(HT$DHW50)) {
  DHWav=seq(0,dhw_range,by=chosen_res) # change to 32 if using 2* more resistant, and 24 if using 1.5x more resistant
  for (j in DHWav) {
    ND = rbind(ND, data.frame(Colony = HT$Colony[i], DHWav = j, pred = PBSpred(j, HT$Intercept[i])))
  }
}
ND$Colony = as.factor(ND$Colony)
ND <- merge(ND, HT[, c("Colony", "DHW50")], by = "Colony")

# bsi_vs_cum_int = ND_Xth_perc[,2:3]
# bsi_only = bsi_vs_cum_int[,2]
# Save the relationship between BSI and MCS cumulative intensity as a text file, for use in SLiM model
# write.table(t(bsi_only), file = "C:/Users/dlcyli/OneDrive/PhD_Code/SLiM/Main codes/bsi_only_mcs_90th.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

P =
  ggplot() + 
  geom_hline(yintercept = 0.5,linetype="dashed",colour="grey50", size = 0.8) +
  geom_line(data = ND,
            aes(x=DHWav, y=pred, 
                group=Colony, colour=DHW50),
            size=1.8, alpha = 0.1) +
  labs(x=bquote(Cumulative~intensity~(degree*C-weeks)),
       y="Bleaching survival probability") +
  scale_colour_gradientn(colours=c("firebrick3","firebrick1","dodgerblue1","dodgerblue3"),
                         breaks = c(0,4,8,12,16,20,24), limits=c(2,24),
                         guide = guide_colourbar(direction = "horizontal", title.position = "top",
                                                 barwidth = 7, barheight = 0.5,
                                                 title = bquote(DHW[50]~(degree*C-weeks)))) +
  scale_x_continuous(expand=c(0,0),breaks = c(0,4,8,12,16,20,24)) + 
  scale_y_continuous(expand=c(0.01,0.01),breaks = c(0,0.25,0.5,0.75,1),labels = c("0","","","","1")) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.title = element_text(hjust=0.5, size=10),
        legend.background = element_blank(),
        legend.position = c(0.15,0.17))
P

rm(ND2)
for (i in 1:length(unique(ND$Colony))) {
  pred_temp = ND$pred[which(ND$Colony == unique(ND$Colony)[i])]
  if (exists("ND2")) {
    ND2 = cbind(ND2, as.data.frame(pred_temp))  
  } else {
    ND2 = as.data.frame(pred_temp)
  }
}

# 2 May 2026 - we are going to test if we can just use a single file
# and then do the niche transformation directly in the model 
# for ND2, rows represent DHW and col represent phenotypes
# say, we want the fitness when DHW = 8 and phen is middle
ND2[161,51] # 0.6894305

# now, shift the curve by sens_mcs = 1.5
# again, we want the fitness when DHW = 8 and phen is middle
ND2[161,51] # 0.9458795

# now, let's go back to the original curve, multiply DHW by 1.5 and see if we get the same answer
8*0.5
which(DHWav == 4)
ND2[81,51] # 0.9800158

# write.table(t(ND2), file = paste0("C:/Users/dlcyli/OneDrive/PhD_Code/SLiM/Main codes/bsi_given_dhw50_MCS", sens_mcs,".txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(t(ND2), file = paste0("C:/Users/dlcyli/OneDrive/PhD_Code/SLiM/Main codes/bsi_given_dhw50_MCS", sens_mcs,".txt"), sep = "\t", row.names = FALSE, col.names = FALSE) # update: 23 Sep 2025 (trying to make the model run faster by changing resolution of parameters)
# test = read.table("C:/Users/dlcyli/OneDrive/PhD_Code/SLiM/Main codes/bsi_given_dhw50_MCS11_res0.05.txt")
# write.table(t(test), file = "C:/Users/dlcyli/OneDrive/PhD_Code/SLiM/Main codes/bsi_given_dhw50_MCS11.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


