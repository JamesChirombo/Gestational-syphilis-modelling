rm(list = ls())

library(sf)
library(tmap)
library(tidyverse)
library(RColorBrewer)
library(spdep)
library(CARBayesST)
library(CARBayes)
library(MASS)
library(coda)
library(lubridate)


source("scripts/analysis_functions.R")

# load required spatial data

#maldistricts <- st_read("data/spatial/mwi_admbnda_adm2.shp")
districts <- st_read("data/spatial/malawi shapes/District.shp") 
 
lakes <- st_read("data/spatial/malawi shapes/lake.shp") %>%
  st_set_crs(32736) 

## syphilis data and population data

dat <- read_csv("data/epidata/syphilis_data_2023.csv")
denom <- read_csv("data/epidata/denominator_data_final.csv")


syphilis_data <- dat %>%
  dplyr::select(-periodid,-perioddescription,-organisationunitid,-organisationunitcode,-organisationunitdescription,-`DHS estimate of ANC 1 coverage`,
         -`RHD ANC (Pre-) Eclampsia Yes`,-`RHD ANC Albendazole 1`,-`RHD ANC CPT Women on CPT`,-`RHD ANC Fe-Fo Tablets 120+`,-`RHD ANC ITN (Bed nets) given Received ITN`,
         -`RHD ANC TTV Doses 2+`,-`RHD ANC NVP Syrup Given Received NVP`,-`RHD ANC SP Doses`,
         -periodcode,-`RHD ANC  HIV Status 1st Visit Previous Negative`,-`RHD ANC  HIV Status 1st Visit Previous Positive`,
         -`RHD ANC 1st Visit ART Status of mother`,-`RHD ANC Final ART Status of mother`) %>%
  rename(mnth_yr = periodname,syphpos = `RHD ANC Syphilis Status Positive`,syphneg = `RHD ANC Syphilis Status Negative`,syphukwn = `RHD ANC Syphilis Status Unknown`) %>%
  rename(anc_vst1 = `RHD ANC Visits per Woman Total with 1 Visit`) %>%
  rename(anc_vst2 = `RHD ANC Visits per Woman Total with 2 Visits`) %>%
  rename(anc_vst3 = `RHD ANC Visits per Woman Total with 3 Visits`) %>%
  rename(anc_vst4 = `RHD ANC Visits per Woman Total with 4 Visits`) %>%
  rename(anc_vst5 = `RHD ANC Visits per Woman Total with 5+ Visits`) %>%
  rename(anc_trm1 = `RHD ANC Week of first ANC Visit Week 0-12`) %>%
  rename(anc_reg = `RHD ANC New Women Registered`) %>%
  rename(anc_hivpos_new = `RHD ANC HIV Test Result New Positive`) %>%
  rename(anc_hivneg_new = `RHD ANC HIV Test Result New Negative`) %>%
  rename(anc_hivnotest = `RHD ANC HIV Test Result Not Done`) %>%
  rename(anc_hivpos_prev = `RHD ANC HIV Test Result Previous Positive`) %>%
  rename(anc_hivneg_prev = `RHD ANC HIV Test Result Previous Negative`) %>%
  rename(anc_hivpos_vst1 = `RHD ANC  HIV Status 1st Visit New Positive`) %>%
  rename(anc_hivneg_vst1 = `RHD ANC  HIV Status 1st Visit New Negative`) %>%
  rename(anc_hivnotest_vst1 = `RHD ANC  HIV Status 1st Visit Not Done`) %>%
  rename(district = organisationunitname) %>%
  filter(mnth_yr != 2023) %>%
  mutate(mnth_yr = lubridate::my(mnth_yr)) %>%
  separate(mnth_yr,into = c("year","month","day"),sep = "-") %>%
  filter(!grepl("MOH",district)) %>%
  filter(district != "Mzimba-North") %>%
  filter(district != "Mzimba-South") %>%
  mutate(district = gsub("-DHO","",district)) %>%
  arrange(district)

# remove Mzimba-North and Mzimba-South
syphilis_data_final <- syphilis_data[-which(syphilis_data$district == "Mzimba-North" | syphilis_data$district == "Mzimba-South"),]
syphilis_data_final$district[syphilis_data_final$district == "Nkhata-Bay"] <- "Nkhatabay"

# process population data

popdata <- denom %>%
  pivot_longer(cols = y2014:y2022, names_to = "year",values_to = "pop") %>%
  mutate(year = gsub("y","", year))

syphdata <- left_join(syphilis_data_final, popdata,by = c("year","district")) %>%
  mutate(distrcode = rep(1:28, each = 108)) %>%
  mutate(monthcode = rep(1:108, 28))


## add additional data from DHS
variableData <- read_csv("data/epidata/Syphilis_model_var1.csv")
variableData2 <- read_csv("data/epidata/Syphilis_model_var2.csv")

syphdata$improved_water <- rep(variableData$`% of population using improved water source`,each = 108)
syphdata$improved_sanit <- rep(variableData$`% of population with access to Improved sanitation1`,each = 108)
syphdata$electricity <- rep(variableData$`% of households with access to electricity`,each = 108)
syphdata$sec_educ <- rep(variableData$`% distribution of women aged 15-49 who had completed more than secondary education`,each = 108)
syphdata$median_educ <- rep(variableData$`Median years of education completed by female household membe aged 15-49r.`,each = 108)
syphdata$literacy <- rep(variableData$`% of women who can read a whole sentence`,each = 108)
syphdata$median_age_sex <- rep(variableData$`Median age of first sexual intercourse among women 20-49`,each = 108)
syphdata$median_age_birth <- rep(variableData$`Median age at first birth in women aged 20-49`,each = 108)
syphdata$perc_live_birth <- rep(variableData$`% of women aged 15-19 who have had a live birth`,each = 108)
syphdata$polygamy <- rep(variableData$`% of currently married men with one or more co-wives`,each = 108)
syphdata$employed <- rep(variableData$`% of women employed in the 12 months
preceding the survey`,each = 108)
syphdata$demand_condom <- rep(variableData$`% who can ask their husband to use a condom`,each = 108)
syphdata$neg_sex_rel <- rep(variableData$`% of women who can negotiate sexual relations`,each = 108)
syphdata$STI <- rep(variableData$`% of women who reported having an STI in last 12 months`,each = 108)
syphdata$men_condom <- rep(variableData$`% of men using condoms`,each = 108)
syphdata$sex_one_partner <- rep(variableData$`% of women aged 15-24 who had 1+ sexual partners in last 12 months`,each = 108)
syphdata$paid_sex <- rep(variableData$`% of men who have ever paid for sexual intercourse`,each = 108)
syphdata$women_more_sexPpartners <- rep(variableData$`% of women aged 15-24 who had 1+ sexual partners in last 12 months`,each = 108)
syphdata$women_HIV_pos <- rep(variableData$`Percentage women  HIV positive`,each = 108)
syphdata$hosp_deliv <- rep(variableData$`% delivered in a health facility`,each = 108)
syphdata$blood_sample_ANC <- rep(variableData$`% who had a blood sample taken during ANC`,each = 108)
syphdata$nopostnatal_check <- rep(variableData$`% with no post-natal check within 2 days of birth`,each = 108)
syphdata$facility_deliv <- rep(variableData$`% delivered in a health facility`,each = 108)
syphdata$skilled_ANC <- rep(variableData$`% receiving Anc from a skilled provider`,each = 108)
syphdata$no_ANC <- rep(variableData$`% who received no ANC`,each = 108)
syphdata$problem_health_care <- rep(variableData$`% of women who reported at least  one problem accessing health care`,each = 108)
syphdata$TFR <- rep(variableData$`Total fertility rate`,each = 108)
syphdata$HTC_ANC <- rep(variableData2$`% of women who received HIV testing and counselling during ANC and delivery`,each = 108)
syphdata$med_months_lastBirth <- rep(variableData2$`Median number of months since last birth`,each = 108)
syphdata$low_BMI <- rep(variableData2$`% of women with BMI<18.5`,each = 108)
syphdata$severe_anaemia <- rep(variableData2$`% of women with severe anaemia`,each = 108)
syphdata$past_sec_educ <- rep(variableData$`% distribution of women aged 15-49 who had completed more than secondary education`,each = 108)


#dat$allwomen <- dat$anc_vst1
syphdata$allwomen <- syphdata$anc_reg # those that registered will be our denominator
#dat$allbirths <- dat$numbabies
syphdata$prev <- syphdata$syphpos/syphdata$allwomen*100
syphdata$syphlisTestingCoverage <- (syphdata$syphpos + syphdata$syphneg)/syphdata$allwomen * 100


casesDistr <- syphdata %>%
  group_by(year,month) %>%
  summarise(totalCases = sum(syphpos,na.rm = TRUE),
            prevdistr = mean(prev,na.rm = TRUE),
            totalwomen = sum(allwomen,na.rm = TRUE)) %>%
  mutate(prev = totalCases/totalwomen*100)

## yearly syphilis prevalence

# summarize differently
case_summary <- syphdata %>%
  group_by(monthcode,year) %>%
  summarize(n = sum(syphpos, na.rm = TRUE)/sum(allwomen, na.rm = TRUE) * 100)

# summarize by time and district
case_summ_district <- syphdata %>%
  group_by(monthcode,district) %>%
  summarise(n = sum(syphpos)/sum(allwomen) * 100)

ggplot(case_summary, aes(x = monthcode, y = n)) + 
  geom_line(size = 0.8) +
  theme_bw() +
  scale_x_continuous(breaks = seq(1,108,12), labels = 2014:2022) +
  ylim(0,10) +
  labs(x = "Year",y = "Prevalence (%)") +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))
ggsave("images/syphilis_temporal_trend.png", width = 250, height = 200, units = "mm")

ggplot(case_summ_district, aes(x = monthcode, y = n, group = district)) +
  geom_line() +
  facet_wrap(~ district) +
  theme_bw() +
  labs(x = "Month", y = "Prevalence") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 17))
ggsave("images/district_syphilis_trend.png", width = 450, height = 450, units = "mm")


# calculation of expected syphilis cases
syphdata$overallProb <- sum(syphdata$syphpos,na.rm = T)/sum(syphdata$pop,na.rm = T) # use population of child bearing age
syphdata$expectedCases <- ceiling(syphdata$overallProb * syphdata$pop)
syphdata$smr <- syphdata$syphpos/syphdata$expectedCases

# summarize and plot smr
smrData <- syphdata %>%
  group_by(year) %>%
  summarise(smrYear = mean(smr,na.r = T))
#plot(smrData$year,smrData$smrYear,type="l",xlab="",ylab="SMR",lwd=2)

# summarize smr by month
smrDataMonth <- syphdata %>%
  group_by(month) %>%
  summarise(smrMonth = mean(smr,na.r = T))

# summarize smr by district
smrDistrict <- syphdata %>%
  group_by(district) %>%
  summarise(smrAvg = sum(syphpos,na.rm = T)/sum(expectedCases,na.rm = T))

#tiff("images/malawi_syphilis_SMR_month.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
#par(mfrow = c(1,1),cex.lab = 1.5,cex.axis = 1.5)
#plot(1:ntime,smrDataMonth$smrMonth,type = "l",axes = F,xlab = "",ylab = "SMR")
#axis(1,at = seq(1,ntime,7),labels = month.abb)
#axis(2,lwd = 1)
#abline(h = 1,lty = 2)
#box(lwd = 1,bty = "o")
#dev.off()

#SMR=rep(NA,28)
#for(i in 1:28){
#  SMR[i]=sum(dat$syphpos[dat$district==i],na.rm=T)/sum(dat$expectedCases[dat$district==i],na.rm=T)
#}

## plot distribution of smr values
# change to modern methods - tmap or ggplot
districts$DISTRICT[districts$DISTRICT == "Nkhata Bay"] <- "Nkhatabay"

# join spatial data with smr data

mw_all_smr <- left_join(districts,smrDistrict,by = c("DISTRICT" = "district"))

overall_smr <- tm_shape(mw_all_smr) +
  tm_polygons(col = "smrAvg", border.col = "black", title = "SMR",textNA = NA, palette = brewer.pal(9,"YlGnBu")) +
  tm_shape(lakes) +
  tm_polygons(col = "lightblue",border.col = "black") +
  tm_layout(legend.text.size = 1.1,
            legend.title.size = 1.2,
            panel.label.size = 1.1,
            legend.outside = FALSE,
            legend.position = c("left","bottom"),
            asp = 0,
            frame.lwd = 1)
tmap_save(overall_smr, "images/SMR_overall.png", width = 200, height = 250, units = "mm")
### map syphilis prevalence among the women cohort ####

prev_summ <- syphdata %>%
  group_by(district,year) %>%
  summarise(syph_prev = sum(syphpos,na.rm = TRUE)/sum(allwomen,na.rm = TRUE)*100)

mw_districts_prev <- left_join(districts,prev_summ,by = c("DISTRICT" = "district")) %>%
  filter(!is.na(year))

# plot prevalence
malawi_prev <- plot_district_prevalence(summary.district.poly = mw_districts_prev, lake.poly = lakes)
tmap_save(malawi_prev,"images/syphilis_district_prevalence.png",width = 400,height = 250,units = "mm")



# summarize smr by year and district
smr_year_all <- syphdata %>%
  group_by(district,year) %>%
  summarize(smr = sum(syphpos,na.rm = TRUE)/sum(expectedCases,na.rm = TRUE))

# merge with spatial data

mw_districs_smr <- left_join(districts,smr_year_all,by = c("DISTRICT" = "district")) %>%
  filter(!is.na(year))

# plot smr
malawi_smr <- plot_district_smr(summary.district.poly = mw_districs_smr,lake.poly = lakes)
tmap_save(malawi_smr, "images/syphilis_district_SMR.png", width = 400, height = 250, units = "mm")

## exploratory model fitting


#############################################################################################
################################ question one ###############################################
## Question 1: Are there any population characteristics associated with geographical differences in syphilis prevalence ##

# non-spatial model
m1 <- glm(syphpos ~ offset(log(expectedCases)) +
            sec_educ + literacy + employed + perc_live_birth + median_age_sex +
            median_age_birth + women_more_sexPpartners + polygamy + men_condom + paid_sex + women_HIV_pos+
            STI + improved_sanit + improved_water + electricity + skilled_ANC +no_ANC + problem_health_care + nopostnatal_check + facility_deliv,
          data = syphdata,family = poisson())
summary(m1)
exp(cbind(IRR = coef(m1),CI = confint(m1)))
AIC(m1)
sjPlot::tab_model(m1)

m2 <- glm(syphpos ~ offset(log(expectedCases)) + employed +
            sec_educ + syphlisTestingCoverage  + electricity +
            median_age_birth + women_more_sexPpartners + women_HIV_pos,
          data = syphdata,family = poisson())
summary(m2)
exp(cbind(IRR = coef(m2),CI = confint(m2)))
AIC(m2)
sjPlot::tab_model(m2)

m3 <- glm(syphpos ~ offset(log(expectedCases)) + sec_educ + no_ANC + median_age_sex + women_more_sexPpartners +
            women_HIV_pos + paid_sex + employed + facility_deliv,
          data = syphdata,family = poisson())
summary(m3)
exp(cbind(IRR = coef(m3),CI = confint(m3)))
AIC(m3)
sjPlot::tab_model(m3)



## regression modelling for lattice data - spatio-temporal modelling ### 
### CAR models 

## first remove cities - they have no data
districts_CAR <- districts %>%
  filter(DISTRICT != "Likoma") %>%
  filter(!(DISTRICT %in% c("Lilongwe City","Blantyre City","Mzuzu City","Zomba City")))

# create neighbourhood matrix
MWnb <- poly2nb(districts_CAR)
mwMat <- nb2mat(MWnb,style = "B", zero.policy = TRUE)

## also remove Likoma Island as it has no neighbours
finaldata <- syphdata %>%
  filter(district != "Likoma")
#finaldata <- filter(dat,district != "likoma")

## remove missing values

finaldata$expectedCases[is.na(finaldata$expectedCases)] <- mean(finaldata$expectedCases,na.rm = TRUE) 
finaldata$TFR[is.na(finaldata$TFR)] <- mean(finaldata$TFR,na.rm = TRUE)
finaldata$syphlisTestingCoverage[is.na(finaldata$syphlisTestingCoverage)] <- mean(finaldata$syphlisTestingCoverage,na.rm = TRUE)

fitdata <- finaldata
fitdata$syphpos[is.na(fitdata$syphpos)] <- 0

formula <- syphpos ~ offset(log(expectedCases)) + employed + sec_educ + syphlisTestingCoverage  + electricity + median_age_birth + women_more_sexPpartners + women_HIV_pos

# fit 3 chains
# chain 1
set.seed(0012)
fit1 <- ST.CARanova(formula = formula,
                  family = "poisson",
                  data = fitdata, 
                  W = mwMat,
                  burnin = 20000,
                  n.sample = 620000,
                  thin = 100,
                  verbose = TRUE)
colnames(fit1$samples$beta) <- c("Intercept","employed","Education","Testing coverage","electricity","Age birth","Sex partners","HIV positive")


# chain 2
set.seed(124)
fit2 <- ST.CARanova(formula = formula,
                    family = "poisson",
                    data = fitdata, 
                    W = mwMat,
                    burnin = 20000,
                    n.sample = 620000,
                    thin = 100,
                    verbose = TRUE)
colnames(fit2$samples$beta) <- c("Intercept","employed","Education","Testing coverage","electricity","Age birth","Sex partners","HIV positive")

# model 3
set.seed(024)
fit3 <- ST.CARanova(formula = formula,
                    family = "poisson",
                    data = fitdata, 
                    W = mwMat,
                    burnin = 20000,
                    n.sample = 620000,
                    thin = 100,
                    verbose = TRUE)
colnames(fit3$samples$beta) <- c("Intercept","employed","Education","Testing coverage","electricity","Age birth","Sex partners","HIV positive")



# incident rate ratios for model 1
#params <- summarise.samples(exp(fit1$samples$beta[,-1]),quantiles = c(0.5,0.025,0.975))
#round(params$quantiles,2)
#summary(fit1)
print(fit1)
#colnames(fit1$samples$beta) <- c("Intercept","employed","Education","Testing coverage","electricity","Age birth","Sex partners","HIV positive")
plot(exp(fit1$samples$beta[,-1]))

# summarize rate ratios
#params <- summarise.samples(exp(fit1$samples$beta[,-1]),quantiles = c(0.5,0.025,0.975))
#round(params$quantiles,2)

# combine chains for further inference

model_samples <- coda::mcmc.list(fit1$samples$beta,fit2$samples$beta,fit3$samples$beta)
plot(model_samples)
coda::gelman.diag(model_samples)

# further processing
model_samples_all <- rbind(fit1$samples$fitted,fit2$samples$fitted,fit3$samples$fitted)
n_samples <- nrow(model_samples_all)
n_all <- ncol(model_samples_all)
risk_samples_combined <- model_samples_all /
  matrix(rep(finaldata$expectedCases, n_samples), nrow = n_samples, ncol = n_all, byrow = TRUE)

# check spatial trends
N <- length(table(finaldata$year))
risk_trends <- array(NA,c(n_samples,N))
for (i in 1:n_samples) {
  risk_trends[i,] <- tapply(risk_samples_combined[i,],finaldata$year,mean)
}

time_trends <- as.data.frame(t(apply(risk_trends,2,quantile,c(0.5, 0.025, 0.975))))
time_trends <- time_trends %>% mutate(Year = names(table(finaldata$year)))
colnames(time_trends)[1:3] <- c("Median","LCI", "UCI")

# plot temporal trends
ggplot(time_trends, aes(x = factor(Year), y = Median, group = 1)) +
  geom_line(col = "black",size = 1) +
  geom_line(aes(x = factor(Year), y = LCI),col = "black",lty = 2) +
  geom_line(aes(x = factor(Year), y = UCI),col = "black",lty = 2) +
  scale_x_discrete(name = "Year", breaks = c(2014,2015,2016,2017,2018,2019,2020,2021,2022),
                   labels = c(2014,2015,2016,2017,2018,2019,2020,2021,2022)) +
  scale_y_continuous(name = "Risk", limits = c(0,2)) +
  theme_bw()  +
  theme(text = element_text(size = 20), 
        plot.title = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 19),
        axis.text.y = element_text(size = 19))
ggsave("images/temporal_risk_trend.png", width = 210, height = 200, units = "mm")


# function to plot exceedance probabilities
yr <- (2014:2022)
exceed_prob_list <- list()
for (i in seq_along(yr)) {
  exceed_prob_list[[i]] <- yearly_posterior_exceedance_prob(year = yr[i])
}

exceed_prob_df <- bind_rows(exceed_prob_list) %>%
  mutate(year = rep(2014:2022,each = 28))

# add exceedance probabilities to the spatial data frame
district_exceed_df <- left_join(districts,exceed_prob_df, by = c("DISTRICT" = "district")) %>%
  filter(!is.na(year))

# plot exceedance probabilities
exprob <- plot_exceedance_probabilities(summary.district.poly = district_exceed_df, lake.poly = lakes)
tmap_save(exprob,"images/exceedance_prob_malawi.png", width = 400, height = 250, units = "mm")


# plot model predicted values
# use fitdata and one of the chains

fitdata$fitted_values <- fit2$fitted.values 
fitdata$fit_smr <- fitdata$fitted_values/fitdata$expectedCases
plot_fit_data <- fitdata %>%
  group_by(district) %>%
  summarise(avgIRR = mean(fit_smr,na.rm = TRUE)) %>%
  add_row(district = "likoma",avgIRR = NA,.after = 9)

# merge plot fit data with spatial data
# this will provide the overall fitted values


# process yearly estimates of fitted values

yearly_fit <- fitdata %>%
  group_by(district,year) %>%
  summarise(irr = mean(fit_smr,na.rm = TRUE)) 

# join with spatial data
yearly_fit_sf <- left_join(districts,yearly_fit, by = c("DISTRICT" = "district")) %>%
  filter(!is.na(year))

predict_irr <- plot_predicted_irr(summary.district.poly = yearly_fit_sf, lake.poly = lakes, break.style = "equal")
tmap_save(predict_irr, "images/yearly_predicted_irr.png", width = 400, height = 250, units = "mm")




# temporal trends 
temporalTrendData <- finaldata %>%
  group_by(year) %>%
  summarize(m1 = mean(r_st[month == 1],na.rm = T),
            m2 = mean(r_st[month == 2],na.rm = T),
            m3 = mean(r_st[month == 3],na.rm = T),
            m4 = mean(r_st[month == 4],na.rm = T))
nmonth <- 12
nyear <- 7
y1_risk <- rep(NA,nmonth)
for (i in 1:nmonth) {
  y1_risk[i] <- mean(finaldata$r_st[finaldata$month == i & finaldata$year == 2019],na.rm = T)
}

# posterior estimates of SIR
y.fit <- fit1$samples$fitted
SIR <- t(t(y.fit)/finaldata$expectedCases)
finaldata$SIR.50 <- apply(SIR, 2, median)
finaldata$SIR.025 <- apply(SIR, 2, quantile, 0.025)
finaldata$SIR.975 <- apply(SIR, 2, quantile, 0.975)
#map$PP <- apply(SIR, 2, function(x) length(which(x > 1))) / M
RR <- finaldata %>%
  group_by(district) %>%
  summarise(irr_2014 = mean(SIR.50[year == 2014],na.rm = TRUE),
            irr_2015 = mean(SIR.50[year == 2015],na.rm = TRUE),
            irr_2016 = mean(SIR.50[year == 2016],na.rm = TRUE),
            irr_2017 = mean(SIR.50[year == 2017],na.rm = TRUE),
            irr_2018 = mean(SIR.50[year == 2018],na.rm = TRUE),
            irr_2019 = mean(SIR.50[year == 2019],na.rm = TRUE),
            irr_2020 = mean(SIR.50[year == 2020],na.rm = TRUE)) %>%
  add_row(district = "likoma",irr_2014 = NA,irr_2015 = NA,irr_2016 = NA,irr_2017 = NA,irr_2018 = NA,irr_2019 = NA,irr_2020 = NA,.after = 9)
## add yearly irr to spatial dataframe
mwdistr$irr_2014 <- RR$irr_2014
mwdistr$irr_2015 <- RR$irr_2015
mwdistr$irr_2016 <- RR$irr_2016
mwdistr$irr_2017 <- RR$irr_2017
mwdistr$irr_2018 <- RR$irr_2018
mwdistr$irr_2019 <- RR$irr_2019
mwdistr$irr_2020 <- RR$irr_2020


# temporal estimates

temporalTrendData <- finaldata %>%
  group_by(year,month) %>%
  summarise(irrMean = median(SIR.50),
            L_IRR = median(SIR.025),
            U_IRR = median(SIR.975))

tiff("images/temporal_RR_estim.tif",width = (30*0.39),height = (20*0.39),units = "in",res = 450,compression = "lzw")
par(mfrow = c(1,1),mar = c(2,4.5,2,2),cex.axis = 1.4,cex.lab = 1.4)
plot(temporalTrendData$month,temporalTrendData$irrMean,type = "l",ylim = c(0,4),lwd = 2,axes = F,ylab = "Risk")
axis(1,at = seq(1,84,12),labels = c(2014:2020))
axis(2,lwd = 1)
polygon(c(temporalTrendData$month,rev(temporalTrendData$month)),c(temporalTrendData$U_IRR,rev(temporalTrendData$L_IRR)),col = "#e8daef",border = NA)
lines(temporalTrendData$month,temporalTrendData$irrMean,lwd = 2,col = "#6f03a5")
#lines(temporalTrendData$month,temporalTrendData$L_IRR,type = "l",lty=2)
#lines(temporalTrendData$month,temporalTrendData$U_IRR,type = "l",lty=2)
abline(h = 1,lty = 2)
box(lwd = 1)
dev.off()

# temporal data for selected districts
temporalTrendData_district <- finaldata %>%
  group_by(year,district) %>%
  summarise(irrMean = median(SIR.50),
            L_IRR = median(SIR.025),
            U_IRR = median(SIR.975))

plot(temporalTrendData_district$year[temporalTrendData_district$district == "balaka"],temporalTrendData_district$irrMean[temporalTrendData_district$district == "balaka"],type = "l",ylim = c(0,2.5),ylab = "Relative Risk")
lines(2014:2020,temporalTrendData_district$L_IRR[temporalTrendData_district$district == "balaka"],type = "l",lty = 2)
lines(2014:2020,temporalTrendData_district$U_IRR[temporalTrendData_district$district == "balaka"],type = "l",lty = 2)


# district trends

border_districts = c("karonga","mchinji","dedza","mwanza","machinga","mulanje")
plot_district_level_risk("karonga")
plot_district_level_risk("balaka")

tiff("images/border_district_risk.tif",width = 35*0.39,height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow = c(3,2),mar = c(4,4.5,2,2))
for(i in 1:length(border_districts)){
  plot_district_level_risk(border_districts[i])
  title(main = simpleCap(border_districts[i]),cex.main = 1.7)
}
# predicted estimates of relative risk
n_samples <- nrow(finaldata)
n_vars <- ncol(fit1$samples$fitted)

RR_samples <- fit1$samples$fitted/matrix(finaldata$expectedCases,nrow = n_samples,ncol = n_vars,byrow = T)



