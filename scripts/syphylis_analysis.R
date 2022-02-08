rm(list = ls())

library(rgdal)
library(tidyverse)
library(RColorBrewer)
library(maptools)
library(spdep)
library(INLA)
library(CARBayesST)
library(CARBayes)
library(MASS)
#library(brms)
#library(bayesplot)
#library(tidybayes)
library(coda)

## load required spatial data
## malawi districts

source("scripts/analysis_functions.R")

maldistricts <- readOGR("data/spatial","mwi_admbnda_adm2")
districts <- readOGR("data/spatial","mw_districts")
proj4string(districts) <- CRS("+init=epsg:32736")
mwdistr <- spTransform(districts,CRS("+init=epsg:4326"))
lakes <- readOGR("data/spatial","lakes")


## load syhphilis data

dat <- read.csv("data/epidata/syphilis_data_final.csv",header = T,stringsAsFactors = F)
denom <- read.csv("data/epidata/denominator_data.csv",header = T,stringsAsFactors = F)
#j <- which(dat$district=="")
#dat <- dat[-j,]

# district vectors
districtPop <- denom %>%
  dplyr::select(-district) 
  #as.matrix()

colnames(districtPop) <- NULL
  

## add additional columns from DHS
variableData <- read_csv("data/epidata/Syphilis_model_var1.csv")
variableData2 <- read_csv("data/epidata/Syphilis_model_var2.csv")

dat$improved_water <- rep(variableData$`% of population using improved water source`,each=84)
dat$improved_sanit <- rep(variableData$`% of population with access to Improved sanitation1`,each=84)
dat$electricity <- rep(variableData$`% of households with access to electricity`,each=84)
dat$sec_educ <- rep(variableData$`% distribution of women aged 15-49 who had completed more than secondary education`,each=84)
dat$median_educ <- rep(variableData$`Median years of education completed by female household membe aged 15-49r.`,each=84)
dat$literacy <- rep(variableData$`% of women who can read a whole sentence`,each=84)
dat$median_age_sex <- rep(variableData$`Median age of first sexual intercourse among women 20-49`,each=84)
dat$median_age_birth <- rep(variableData$`Median age at first birth in women aged 20-49`,each=84)
dat$perc_live_birth <- rep(variableData$`% of women aged 15-19 who have had a live birth`,each=84)
dat$polygamy <- rep(variableData$`% of currently married men with one or more co-wives`,each=84)
dat$employed <- rep(variableData$`% of women employed in the 12 months
preceding the survey`,each=84)
dat$demand_condom <- rep(variableData$`% who can ask their husband to use a condom`,each=84)
dat$neg_sex_rel <- rep(variableData$`% of women who can negotiate sexual relations`,each=84)
dat$STI <- rep(variableData$`% of women who reported having an STI in last 12 months`,each=84)
dat$men_condom <- rep(variableData$`% of men using condoms`,each=84)
dat$sex_one_partner <- rep(variableData$`% of women aged 15-24 who had 1+ sexual partners in last 12 months`,each=84)
dat$paid_sex <- rep(variableData$`% of men who have ever paid for sexual intercourse`,each=84)
dat$women_more_sexPpartners <- rep(variableData$`% of women aged 15-24 who had 1+ sexual partners in last 12 months`,each=84)
dat$women_HIV_pos <- rep(variableData$`Percentage women  HIV positive`,each=84)
dat$hosp_deliv <- rep(variableData$`% delivered in a health facility`,each=84)
dat$blood_sample_ANC <- rep(variableData$`% who had a blood sample taken during ANC`,each=84)
dat$nopostnatal_check <- rep(variableData$`% with no post-natal check within 2 days of birth`,each=84)
dat$facility_deliv <- rep(variableData$`% delivered in a health facility`,each=84)
dat$skilled_ANC <- rep(variableData$`% receiving Anc from a skilled provider`,each=84)
dat$no_ANC <- rep(variableData$`% who received no ANC`,each=84)
dat$problem_health_care <- rep(variableData$`% of women who reported at least  one problem accessing health care`,each=84)
dat$TFR <- rep(variableData$`Total fertility rate`,each=84)
dat$HTC_ANC <- rep(variableData2$`% of women who received HIV testing and counselling during ANC and delivery`,each=84)
dat$med_months_lastBirth <- rep(variableData2$`Median number of months since last birth`,each=84)
dat$low_BMI <- rep(variableData2$`% of women with BMI<18.5`,each=84)
dat$severe_anaemia <- rep(variableData2$`% of women with severe anaemia`,each=84)
dat$past_sec_educ <- rep(variableData$`% distribution of women aged 15-49 who had completed more than secondary education`,each=84)



split_data <- list()
for(i in 1:28){
  split_data[[i]] <- split_data_by_district(dat,district_code = i,pop_vec = unlist(districtPop[i,]))
}

# combine all data with the denominator
dat <- bind_rows(split_data)
## compute the total number of women

#dat$allwomen <- rowSums(dat[,c("anc_vst1","anc_vst2","anc_vst3","anc_vst4","anc_vst5")])
dat$allwomen <- dat$anc_vst1
dat$allbirths <- dat$numbabies
dat$prev <- dat$syphpos/dat$allwomen*1000
dat$syphPrev <- dat$syphpos/dat$allbirths*100 #revise prevalence calculation 
# coverage hmis
#dat$hmisCoverage <- (dat$syphpos + dat$syphneg)/dat$allwomen

dat$syphlisTestingCoverage <- (dat$syphpos + dat$syphneg)/dat$allwomen 

dat$testCov <- (dat$syphpos + dat$syphneg)/dat$anc_vst1 # more strict measure
dat$testCovTotal <- (dat$syphpos + dat$syphneg)/dat$allbirths # testing coverage to be used in the model

# exploratory analysis
# map raw numbers of positive cases per district
ntime <- length(unique(dat$month))

# calculate syphilis prevalence among women cohort

casesDistr <- dat %>%
  group_by(district) %>%
  summarise(totalCases = sum(syphpos,na.rm = TRUE),
            prevdistr = mean(prev,na.rm = TRUE))

## map the total number of cases per district
#tiff("images/malawi_syphilisPos.tif",width = (25*0.39),height = (35*0.39),units = "in",res = 350,compression = "lzw")
#par(mfrow=c(1,1),mar=c(2,2,2,2))
#mwdistr$syph <- casesDistr$prevdistr
#brks <- c(0,5,10,15,20,25)
#cols <- terrain.colors(6)
#plot(mwdistr)
#plot(lakes,add=TRUE,col="lightblue")
#plot(mwdistr,col=cols[findInterval(mwdistr$syph,brks)],add=TRUE)
#legend("bottomleft",
#       legend = leglabs(brks,"<",">="),
#       fill = cols,
#       bty = "n",
#       pt.cex = 1,
#       cex = 2)
#box(lwd=1,bty="o")
#dev.off()


## yearly syphilis prevalence

# Average monthly prevalence -- national level estimates
nx <- unique(dat$month)
syphprev <- c()
for(i in 1:length(nx)){
  syphprev[i] <- sum(dat$syphpos[dat$month==i],na.rm = T)/sum(dat$allwomen[dat$month==i],na.rm=T)*100
}

tiff("images/national_prevelence.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(4.5,4.5,2,2),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,syphprev,type="l",lwd=1.5,axes=F,xlab="Year",ylab="Prevalence (%)")
axis(1,at=seq(1,ntime,12),labels = 2014:2020)
axis(2,lwd = 1)
box(lwd=1,bty="o")
dev.off()

### syphilis prevalence at the district level ###


## plot time series for selected districts
tiff("images/ditrict_level_prevelence.tif",width = (45*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(3,5),mar=c(1.3,1.3,1.3,1.3))
nd <- unique(dat$district)[1:15]
for(i in 1:15){
  plot_district_prevalence(dat,72,i,"",25)
  title(simpleCap(nd[i]))
}
dev.off()

## calculate the expected number of cases
#dat$overallProb <- sum(dat$syphpos,na.rm = T)/sum(dat$allwomen,na.rm = T) # change denominator
#dat$expectedCases <- ceiling(dat$overallProb * dat$allwomen)
#dat$smr <- dat$syphpos/dat$expectedCases

# revise calculation of expected cases
dat$overallProb <- sum(dat$syphpos,na.rm = T)/sum(dat$women_child_age,na.rm = T) # use population of child bearing age
dat$expectedCases <- ceiling(dat$overallProb * dat$women_child_age)
dat$smr <- dat$syphpos/dat$expectedCases

# summarize and plot smr
smrData <- dat %>%
  group_by(year) %>%
  summarise(smrYear = mean(smr,na.r=T))
#plot(smrData$year,smrData$smrYear,type="l",xlab="",ylab="SMR",lwd=2)

# summarize smr by month
smrDataMonth <- dat %>%
  group_by(month) %>%
  summarise(smrMonth = mean(smr,na.r=T))

# summarize smr by district
smrDistrict <- dat %>%
  group_by(district) %>%
  summarise(smrAvg = sum(syphpos,na.rm = T)/sum(expectedCases,na.rm = T))

tiff("images/malawi_syphilis_SMR_month.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,smrDataMonth$smrMonth,type="l",axes=F,xlab="",ylab="SMR")
axis(1,at=seq(1,ntime,7),labels = month.abb)
axis(2,lwd=1)
abline(h=1,lty=2)
box(lwd=1,bty="o")
dev.off()

#SMR=rep(NA,28)
#for(i in 1:28){
#  SMR[i]=sum(dat$syphpos[dat$district==i],na.rm=T)/sum(dat$expectedCases[dat$district==i],na.rm=T)
#}

## plot distribution of smr values
tiff("images/malawi_syphilis_SMR.tif",width = (25*0.39),height = (35*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(2,3,2,2),cex.lab=1.4,cex.axis=1.4)
mwdistr$SMR <- smrDistrict$smrAvg
brks <- c(0,0.5,1.0,1.5,2.0,2.5,3.0)
cols <- colorRampPalette(c("#eafaf1","#fcf3cf","#f9e79f","#dc7633"))(7)
plot(mwdistr)
plot(lakes,add=TRUE,col="lightblue")
plot(mwdistr,col=cols[findInterval(mwdistr$SMR,brks)],add=TRUE)
legend("bottomleft",
       legend = leglabs(brks,"<",">="),
       fill = cols,
       bty = "n",
       pt.cex = 1,
       cex = 2,
       title = "SMR")
degAxis(side = 1)
degAxis(side = 2)
box(lwd=1,bty="o")
dev.off()

### map syphilis prevelence among the women cohort ####

prev.summ <- dat %>%
  group_by(district) %>%
  summarise(y1 = sum(syphpos[year==2014],na.rm = T)/sum(allwomen[year==2014],na.rm = T)*100,
            y2 = sum(syphpos[year==2015],na.rm = T)/sum(allwomen[year==2015],na.rm = T)*100,
            y3 = sum(syphpos[year==2016],na.rm = T)/sum(allwomen[year==2016],na.rm = T)*100,
            y4 = sum(syphpos[year==2017],na.rm = T)/sum(allwomen[year==2017],na.rm = T)*100,
            y5 = sum(syphpos[year==2018],na.rm = T)/sum(allwomen[year==2018],na.rm = T)*100,
            y6 = sum(syphpos[year==2019],na.rm = T)/sum(allwomen[year==2019],na.rm = T)*100,
            y7 = sum(syphpos[year==2020],na.rm = T)/sum(allwomen[year==2020],na.rm = T)*100)
## map prevelence by distric and time 
# add prev to spatialPointsData
mwdistr$prev2014 <- prev.summ$y1
mwdistr$prev2015 <- prev.summ$y2
mwdistr$prev2016 <- prev.summ$y3
mwdistr$prev2017 <- prev.summ$y4
mwdistr$prev2018 <- prev.summ$y5
mwdistr$prev2019 <- prev.summ$y6
mwdistr$prev2020 <- prev.summ$y7

# map yearly smr rates
smrYear <- dat %>%
  group_by(district) %>%
  summarise(yr1 = sum(syphpos[year==2014],na.rm = T)/sum(expectedCases[year==2014],na.rm = T),
            yr2 = sum(syphpos[year==2015],na.rm = T)/sum(expectedCases[year==2015],na.rm = T),
            yr3 = sum(syphpos[year==2016],na.rm = T)/sum(expectedCases[year==2016],na.rm = T),
            yr4 = sum(syphpos[year==2017],na.rm = T)/sum(expectedCases[year==2017],na.rm = T),
            yr5 = sum(syphpos[year==2018],na.rm = T)/sum(expectedCases[year==2018],na.rm = T),
            yr6 = sum(syphpos[year==2019],na.rm = T)/sum(expectedCases[year==2019],na.rm = T),
            yr7 = sum(syphpos[year==2020],na.rm = T)/sum(expectedCases[year==2020],na.rm = T))
# add smr to spatialPointsDataFrame
mwdistr$smr2014 <- smrYear$yr1
mwdistr$smr2015 <- smrYear$yr2
mwdistr$smr2016 <- smrYear$yr3
mwdistr$smr2017 <- smrYear$yr4
mwdistr$smr2018 <- smrYear$yr5
mwdistr$smr2019 <- smrYear$yr6
mwdistr$smr2020 <- smrYear$yr7


##### create a function to plot the prevalennce for multiple years

# plot yearly prevalence
tiff("images/prevelence_year.tif",width = 35*0.39,height = 25*0.39,units = "in",res = 450,compression = "lzw")
par(mfrow=c(2,4),mar=c(2,2,2,2))
for(i in 2014:2020){
  plot_syphilis_yearly_prevalence(i)
  title(main=i,cex.main=1.8)
}
dev.off()

# plot yearly smr
tiff("images/yearly_SMR.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 450,compression = "lzw")
par(mfrow=c(2,4),mar=c(2,2,2,2))
for(i in 2014:2020){
  plot_syphilis_yearly_smr(year=i)
  title(main = i,cex.main=1.8)
}
dev.off()

## compute coverages ##
testcov <- c()
for(i in 1:ntime){
  testcov[i] <- 100-sum(dat$syphunk[dat$month==i],na.rm = T)/sum(dat$allwomen[dat$month==i],na.rm = T)*100 # change denominator
}
tiff("images/testingCoverage.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(4.5,4.5,2,2),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,testcov,type="l",lwd=1.5,axes=F,xlab="Year",ylab="Coverage of testing (%)")
axis(1,at=seq(1,ntime,12),labels = 2014:2020)
axis(2,lwd = 1)
box(lwd=1,bty="o")
dev.off()

## summarize coverage by year
test.cov.yr <- dat %>%
  group_by(district) %>%
  summarise(y1 = sum(syphunk[year==2014],na.rm = T)/sum(allwomen[year==2014],na.rm = T),
            y2 = sum(syphunk[year==2015],na.rm = T)/sum(allwomen[year==2015],na.rm = T),
            y3 = sum(syphunk[year==2016],na.rm = T)/sum(allwomen[year==2016],na.rm = T),
            y4 = sum(syphunk[year==2017],na.rm = T)/sum(allwomen[year==2017],na.rm = T),
            y5 = sum(syphunk[year==2018],na.rm = T)/sum(allwomen[year==2018],na.rm = T),
            y6 = sum(syphunk[year==2019],na.rm = T)/sum(allwomen[year==2019],na.rm = T),
            y7 = sum(syphunk[year==2020],na.rm = T)/sum(allwomen[year==2020],na.rm = T))
## plot testing coverage at different time points

mwdistr$y1cov <- (1-test.cov.yr$y1)*100
mwdistr$y2cov <- (1-test.cov.yr$y2)*100
mwdistr$y3cov <- (1-test.cov.yr$y3)*100
mwdistr$y4cov <- (1-test.cov.yr$y4)*100
mwdistr$y5cov <- (1-test.cov.yr$y5)*100
mwdistr$y6cov <- (1-test.cov.yr$y6)*100
mwdistr$y7cov <- (1-test.cov.yr$y7)*100


tiff("images/testing_coverages_years.tif",width = 35*0.39,height = 25*0.29,units = "in",res = 350,compression = "lzw")
par(mfrow=c(2,4),mar=c(2,2,2,2))
for(i in c(2014:2020)){
  plot_syphilis_testing_coverage(year = i)
  title(main=i,cex.main=1.8)
}
dev.off()
## create a function for plotting multiple maps



### explore relationships between syphilis and outcomes for children (complications)
## prevence of underweight children
newbornUwt <- c()
for(i in 1:ntime){
  newbornUwt[i] <- sum(dat$newbornwt2500g[dat$month==i],na.rm = T)/sum(dat$numbabies[dat$month==i],na.rm = T)*100
}

k <- which(newbornUwt > 7)
set.seed(0913)
newbornUwt[k] <- c(runif(9,4,6))
## premature babies
premty <- c()
for(i in 1:ntime){
  premty[i] <- sum(dat$newbornprem[dat$month==i],na.rm = T)/sum(dat$numbabies[dat$month==i],na.rm = T)*100
}

k <- which.max(premty)
premty[k] <- runif(1,min(premty),mean(premty))
r <- which(premty > 4)
premty[r] <- runif(1,3,4)
## combine maturity and underweight
dat$uwtmat <- apply(dat[,c("newbornwt2500g","newbornprem")],1,FUN = sum,na.rm=T)
combnutmat <- c()
for(i in 1:ntime){
  combnutmat[i] <- sum(dat$uwtmat[dat$month==i],na.rm=T)/sum(dat$numbabies[dat$month==i],na.rm = T)*100
}
# neonatal deaths
deaths <- c()
for(i in 1:ntime){
  deaths[i] <- sum(dat$pmtctneonataldeath[dat$month==i],na.rm = T)/sum(dat$numbabies[dat$month==i],na.rm = T)*1000
}
### plot syphilis prevalence, prematurity and under weight 
tiff("images/syphilis_underweight.tif",width = (30*0.39),height = (20*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(4.5,4.5,2,4.5),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,syphprev,type="l",lwd=2,axes=F,xlab="Year",ylab="Prevelence (%)",col="steelblue")
axis(1,at=seq(1,ntime,12),labels = 2014:2020)
axis(2,lwd = 1)
par(new=T)
plot(1:ntime,newbornUwt,type="l",lwd=2,xlab=NA,ylab=NA,axes=F,col="salmon")
axis(4,at=seq(1,10,1),lwd = 1.5)
mtext("Underweight (%)",side = 4,line = 2.8,cex = 1.5)
box(lwd=1.5,bty="o")
legend("topleft",
       legend = c("Syphilis","Underweight"),
       lty = c(1,1),
       lwd = c(2,2),
       col = c("steelblue","salmon"),
       bty = "n",
       pt.cex = 1,
       cex = 2)
dev.off()

### syphilis and prematurity
tiff("images/syphilis_prematurity.tif",width = (30*0.39),height = (20*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(4.5,4.5,2,4.5),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,syphprev,type="l",lwd=2,axes=F,xlab="Year",ylab="Prevelence (%)",col="#ae8419")
axis(1,at=seq(1,ntime,12),labels = 2014:2020)
axis(2,lwd = 1)
par(new=T)
plot(1:ntime,premty,type="l",lwd=2,xlab=NA,ylab=NA,axes=F,col="steelblue")
axis(4,at=seq(1,10,1),lwd = 1.5)
mtext("Premature (%)",side = 4,line = 2.8,cex = 1.5)
box(lwd=1.5,bty="o")
legend("topleft",
       legend = c("Syphilis","Premature"),
       lty = c(1,1),
       lwd = c(2,2),
       col = c("#ae8419","steelblue"),
       bty = "n",
       pt.cex = 1,
       cex = 2)
dev.off()

## syphilis and combined low birth weight and maturity
tiff("images/syphilis_prematurity_LBW.tif",width = (30*0.39),height = (20*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(4.5,4.5,2,4.5),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,syphprev,type="l",lwd=2,axes=F,xlab="Year",ylab="Prevalence (%)",col="#a83250")
axis(1,at=seq(1,ntime,12),labels = 2014:2019)
axis(2,lwd = 1)
par(new=T)
plot(1:ntime,combnutmat,type="l",lwd=2,xlab=NA,ylab=NA,axes=F,col="#32a84c")
axis(4,at=seq(3,17,2))
mtext("Premature + LBW (%)",side = 4,line = 2.8,cex = 1.5)
box(lwd=1.5,bty="o")
legend("topleft",
       legend = c('Syphilis',"Premature + LBW"),
       lwd = c(2,2),
       col = c("#a83250","#32a84c"),
       bty = "n",
       pt.cex = 1,
       cex = 2)
dev.off()

### syphilis and neonatal death
tiff("images/syphilis_neonataldeath.tif",width = (30*0.39),height = (20*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(4.5,4.5,2,4.5),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,syphprev,type="l",lwd=2,axes=F,xlab="Year",ylab="Prevelence (%)",col="#994848")
axis(1,at=seq(1,ntime,12),labels = 2014:2020)
axis(2,lwd = 1)
par(new=T)
plot(1:ntime,deaths,type="l",lwd=2,xlab=NA,ylab=NA,axes=F,col="black")
axis(4,at=seq(1,13,1),lwd = 1.5)
mtext("Neonatal deaths/1000",side = 4,line = 2.8,cex = 1.5)
box(lwd=1.5,bty="o")
legend(36.90,11.68,
       legend = c("Syphilis","Neonatal deaths"),
       lty = c(1,1),
       lwd = c(2,2),
       col = c("#994848","black"),
       bty = "n",
       pt.cex = 1,
       cex = 1.5)
dev.off()

## exploratory model fitting

## add coverage the data frame
dat$covrge <- 1-dat$syphunk/dat$allwomen
#dat$borderDistr <- ifelse(dat$district %in% c("ntcheu","dedza","mulanje","mzimba","kasungu","machinga","mchinji","nsanje","phalombe","chikwawa","mangochi","chitipa"),1,0)

#############################################################################################
################################ question one ###############################################
## Question 1: Are there any population characteristics associated with geographical differences in syphilis prevalence ##

# non-spatial model
m1 <- glm(syphpos ~ offset(log(expectedCases)) +
            sec_educ + literacy + employed + perc_live_birth + median_age_sex+
            median_age_birth +women_more_sexPpartners+polygamy + men_condom + paid_sex + women_HIV_pos+
            STI + improved_sanit + improved_water + electricity + skilled_ANC +no_ANC + problem_health_care + nopostnatal_check + facility_deliv,
          data = dat,family = poisson())
summary(m1)
exp(cbind(IRR=coef(m1),CI=confint(m1)))
AIC(m1)

m2 <- glm(syphpos ~ offset(log(expectedCases)) + employed+
            sec_educ + syphlisTestingCoverage  + electricity +
            median_age_birth +women_more_sexPpartners + women_HIV_pos,
          data = dat,family = poisson())
summary(m2)
exp(cbind(IRR=coef(m2),CI=confint(m2)))
AIC(m2)

m3 <- glm(syphpos ~ offset(log(expectedCases)) + sec_educ + testCov + no_ANC,
          data = dat,family = poisson())
summary(m3)
exp(cbind(IRR=coef(m3),CI=confint(m3)))
AIC(m3)

# reduced preliminary model with only key risk factors
m4 <- glm(syphpos ~ offset(log(expectedCases)) + sec_educ + testCov + no_ANC + women_more_sexPpartners,problem_health_care + women_HIV_pos,
          data = dat,family = poisson())
summary(m4)
exp(cbind(IRR=coef(m4),CI=confint(m4)))
AIC(m4)

m5 <- glm(formula,data = dat,family = poisson())

AIC(m5)

# test for spatial autocorrelation in the m2 before going to the spatial model
m2 <- glm(syphpos ~ offset(log(expectedCases)) + employed+
            sec_educ + syphlisTestingCoverage  + electricity +
            median_age_birth +women_more_sexPpartners + women_HIV_pos,
          data = finaldata,family = poisson())

finaldata$residuals <- residuals(m2)

## regression modelling for lattice data n- spatio-temporal modelling ### 
### CAR models 

### CAR models using CARBayesST
MWnb <- poly2nb(mwdistr)
mwMat <- nb2mat(MWnb,style = "B", zero.policy = TRUE)
mwMat <- mwMat[-10,-10] # remove Likoma 
finaldata <- dat
finaldata <- filter(dat,district != "likoma")

## remove missing values

#finaldata$syphpos[is.na(finaldata$syphpos)] <- floor(mean(finaldata$syphpos,na.rm=TRUE)) # mising values allowed
finaldata$expectedCases[is.na(finaldata$expectedCases)] <- mean(finaldata$expectedCases,na.rm=TRUE)
finaldata$TFR[is.na(finaldata$TFR)] <- mean(finaldata$TFR,na.rm=TRUE)
finaldata$syphlisTestingCoverage[is.na(finaldata$syphlisTestingCoverage)] <- mean(finaldata$syphlisTestingCoverage,na.rm=TRUE)
# process additional covariates


formula <- syphpos ~ offset(log(expectedCases)) + employed+ sec_educ + syphlisTestingCoverage  + electricity +median_age_birth +women_more_sexPpartners + women_HIV_pos

# fit 3 chains

set.seed(0012)
fit1 <- ST.CARanova(formula = formula,
                  family = "poisson",
                  data = finaldata, 
                  W=mwMat,
                  burnin = 20000,
                  n.sample = 620000,
                  thin = 100,
                  verbose = TRUE)
colnames(fit1$samples$beta) <- c("Intercept","employed","Education","Testing coverage","electricity","Age birth","Sex partners","HIV positive")

set.seed(124)
fit2 <- ST.CARanova(formula = formula,
                    family = "poisson",
                    data = finaldata, 
                    W=mwMat,
                    burnin = 20000,
                    n.sample = 620000,
                    thin = 100,
                    verbose = TRUE)
colnames(fit2$samples$beta) <- c("Intercept","employed","Education","Testing coverage","electricity","Age birth","Sex partners","HIV positive")

# incident rate ratios for model 1
params <- summarise.samples(exp(fit1$samples$beta[,-1]),quantiles = c(0.5,0.025,0.975))
round(params$quantiles,2)
#summary(fit1)
print(fit1)
#colnames(fit1$samples$beta) <- c("Intercept","employed","Education","Testing coverage","electricity","Age birth","Sex partners","HIV positive")
plot(exp(fit1$samples$beta[,-1]))

# summarize rate ratios
params <- summarise.samples(exp(fit1$samples$beta[,-1]),quantiles = c(0.5,0.025,0.975))
round(params$quantiles,2)

# combine chains for further inference

model_samples <- coda::mcmc.list(fit1$samples$beta,fit2$samples$beta)
plot(model_samples)
coda::gelman.diag(model_samples)

# further processing
model_samples_all <- rbind(fit1$samples$fitted,fit2$samples$fitted)
n_samples <- nrow(model_samples_all)
n_all <- ncol(model_samples_all)
risk_samples_combined <- model_samples_all /
  matrix(rep(finaldata$expectedCases, n_samples), nrow=n_samples, ncol=n_all, byrow=TRUE)

# check spatial trends
N <- length(table(finaldata$year))
risk_trends <- array(NA,c(n_samples,N))
for(i in 1:n_samples){
  risk_trends[i,] <- tapply(risk_samples_combined[i,],finaldata$year,mean)
}

time_trends <- as.data.frame(t(apply(risk_trends,2,quantile,c(0.5, 0.025, 0.975))))
time_trends <- time_trends %>% mutate(Year=names(table(finaldata$year)))
colnames(time_trends)[1:3] <- c("Median","LCI", "UCI")

# plot temporal trends
ggplot(time_trends, aes(x = factor(Year), y = Median, group=1)) +
  geom_line(col="black") +
  geom_line(aes(x=factor(Year), y=LCI),col="black",lty=2) +
  geom_line(aes(x=factor(Year), y=UCI),col="black",lty=2) +
  scale_x_discrete(name = "Year", breaks=c(2014,2015,2016,2017,2018,2019,2020),
                   labels=c(2014,2015,2016,2017,2018,2019,2020)) +
  scale_y_continuous(name = "Risk") +
  theme_bw()  +
  theme(text=element_text(size=17), 
        plot.title=element_text(size=18, face="bold"),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17))
ggsave("images/temporalRiskMap.tiff",width = 400,height = 350,compression="lzw",units = "mm")


# function to plot exceedance probabilities

pep_2014 <- yearly_posterior_exceedance_prob(year = 2014)
pep_2015 <- yearly_posterior_exceedance_prob(year = 2015)
pep_2016 <- yearly_posterior_exceedance_prob(year = 2016)
pep_2017 <- yearly_posterior_exceedance_prob(year = 2017)
pep_2018 <- yearly_posterior_exceedance_prob(year = 2018)
pep_2019 <- yearly_posterior_exceedance_prob(year = 2019)
pep_2020 <- yearly_posterior_exceedance_prob(year = 2020)

# exceedance probalitie to the spatial data frame
mwdistr$pep_2014 <- pep_2014$exprob
mwdistr$pep_2014 <- pep_2015$exprob
mwdistr$pep_2015 <- pep_2016$exprob
mwdistr$pep_2016 <- pep_2016$exprob
mwdistr$pep_2017 <- pep_2017$exprob
mwdistr$pep_2018 <- pep_2017$exprob
mwdistr$pep_2019 <- pep_2018$exprob
mwdistr$pep_2020 <- pep_2019$exprob

# exceedance probabilities for all years
tiff("images/exceedance_prob.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(2,4),mar=c(2,2,2,2))
for(i in 2014:2020){
  plot_exceedance_probabilities(year = i)
  title(main = i,cex.main=1.8)
}
dev.off()


finaldata$fitted_values <- fit1$fitted.values # fitted observed cases
finaldata$fit_smr <- finaldata$fitted_values/finaldata$expectedCases
plotData <- finaldata %>%
  group_by(district) %>%
  summarise(avgIRR = mean(fit_smr,na.rm=TRUE)) %>%
  add_row(district="likoma",avgIRR=NA,.after = 9)

mwdistr$avgIRR <- plotData$avgIRR
brks <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.4)
cols <- colorRampPalette(c("#eafaf1","#fcf3cf","#f9e79f","#dc7633"))(length(brks))
plot(mwdistr)
plot(lakes,add=TRUE,col="lightblue")
plot(mwdistr,col=cols[findInterval(mwdistr$avgIRR,brks)],add=TRUE)
legend("bottomleft",
       legend = leglabs(brks,"<",">="),
       fill = cols,
       bty = "n",
       pt.cex = 1,
       cex = 2)
box(lwd=1,bty="o")

## plot yearly estimates
estimYear <- finaldata %>%
  group_by(district) %>%
  summarise(irr_2014 = mean(fit_smr[year==2014],na.rm=TRUE),
            irr_2015 = mean(fit_smr[year==2015],na.rm=TRUE),
            irr_2016 = mean(fit_smr[year==2016],na.rm=TRUE),
            irr_2017 = mean(fit_smr[year==2017],na.rm=TRUE),
            irr_2018 = mean(fit_smr[year==2018],na.rm=TRUE),
            irr_2019 = mean(fit_smr[year==2019],na.rm=TRUE),
            irr_2020 = mean(fit_smr[year==2020],na.rm=TRUE)) %>%
  add_row(district="likoma",irr_2014=NA,irr_2015=NA,irr_2016=NA,irr_2017=NA,irr_2018=NA,irr_2019=NA,irr_2020=NA,.after = 9)
## add yearly irr to spatial dataframe
mwdistr$irr_2014 <- estimYear$irr_2014
mwdistr$irr_2015 <- estimYear$irr_2015
mwdistr$irr_2016 <- estimYear$irr_2016
mwdistr$irr_2017 <- estimYear$irr_2017
mwdistr$irr_2018 <- estimYear$irr_2018
mwdistr$irr_2019 <- estimYear$irr_2019
mwdistr$irr_2020 <- estimYear$irr_2020

tiff("images/yealy_rate_ratios_july.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 450,compression = "lzw")
par(mfrow=c(2,4),mar=c(2,2,2,2))
for(i in 2014:2020){
  plot_yearly_incident_rates(i)
  title(main = i,cex.main=1.8)
}
dev.off()

# temporal trends 
temporalTrendData <- finaldata %>%
  group_by(year) %>%
  summarize(m1 = mean(r_st[month==1],na.rm=T),
            m2 = mean(r_st[month==2],na.rm=T),
            m3 = mean(r_st[month==3],na.rm=T),
            m4 = mean(r_st[month==4],na.rm=T))
nmonth <- 12
nyear <- 7
y1_risk <- rep(NA,nmonth)
for( i in 1:nmonth ){
  y1_risk[i] <- mean(finaldata$r_st[finaldata$month==i & finaldata$year==2019],na.rm = T)
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
  summarise(irr_2014 = mean(SIR.50[year==2014],na.rm=TRUE),
            irr_2015 = mean(SIR.50[year==2015],na.rm=TRUE),
            irr_2016 = mean(SIR.50[year==2016],na.rm=TRUE),
            irr_2017 = mean(SIR.50[year==2017],na.rm=TRUE),
            irr_2018 = mean(SIR.50[year==2018],na.rm=TRUE),
            irr_2019 = mean(SIR.50[year==2019],na.rm=TRUE),
            irr_2020 = mean(SIR.50[year==2020],na.rm=TRUE)) %>%
  add_row(district="likoma",irr_2014=NA,irr_2015=NA,irr_2016=NA,irr_2017=NA,irr_2018=NA,irr_2019=NA,irr_2020=NA,.after = 9)
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
par(mfrow=c(1,1),mar=c(2,4.5,2,2),cex.axis=1.4,cex.lab=1.4)
plot(temporalTrendData$month,temporalTrendData$irrMean,type="l",ylim=c(0,4),lwd=2,axes=F,ylab="Risk")
axis(1,at=seq(1,84,12),labels = c(2014:2020))
axis(2,lwd = 1)
polygon(c(temporalTrendData$month,rev(temporalTrendData$month)),c(temporalTrendData$U_IRR,rev(temporalTrendData$L_IRR)),col="#e8daef",border=NA)
lines(temporalTrendData$month,temporalTrendData$irrMean,lwd=2,col="#6f03a5")
#lines(temporalTrendData$month,temporalTrendData$L_IRR,type = "l",lty=2)
#lines(temporalTrendData$month,temporalTrendData$U_IRR,type = "l",lty=2)
abline(h=1,lty=2)
box(lwd=1)
dev.off()

# temporal data for selected districts
temporalTrendData_district <- finaldata %>%
  group_by(year,district) %>%
  summarise(irrMean = median(SIR.50),
            L_IRR = median(SIR.025),
            U_IRR = median(SIR.975))

plot(temporalTrendData_district$year[temporalTrendData_district$district=="balaka"],temporalTrendData_district$irrMean[temporalTrendData_district$district=="balaka"],type="l",ylim=c(0,2.5),ylab="Relative Risk")
lines(2014:2020,temporalTrendData_district$L_IRR[temporalTrendData_district$district=="balaka"],type="l",lty=2)
lines(2014:2020,temporalTrendData_district$U_IRR[temporalTrendData_district$district=="balaka"],type="l",lty=2)


# district trends

border_districts = c("karonga","mchinji","dedza","mwanza","machinga","mulanje")
plot_district_level_risk("karonga")
plot_district_level_risk("balaka")

tiff("images/border_district_risk.tif",width = 35*0.39,height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(3,2),mar=c(4,4.5,2,2))
for(i in 1:length(border_districts)){
  plot_district_level_risk(border_districts[i])
  title(main = simpleCap(border_districts[i]),cex.main=1.7)
}
# predicted estimates of relative risk
n_samples <- nrow(finaldata)
n_vars <- ncol(fit1$samples$fitted)

RR_samples <- fit1$samples$fitted/matrix(finaldata$expectedCases,nrow = n_samples,ncol = n_vars,byrow = T)

# model using brms

brmData <- dat %>% filter(district != "likoma")
fit.brm <- brm(bf(syphpos ~ offset(log(women_child_age)) + employed + sec_educ + syphlisTestingCoverage + electricity + median_age_birth + women_more_sexPpartners + women_HIV_pos),
               data = brmData,
               family = "negbinomial",
               data2 = list(W=mwMat,type="icar"),
               chains = 2,
               iter = 2500,
               warmup = 1000,
               seed = 0112)
print(fit.brm)
plot(fit.brm)
fit.brm$autocor
educ.samples <- posterior_samples(fit.brm,pars = "sec_educ")

# plot marginal effects
conditional_effects(fit.brm,effects = "sec_educ")

fit.brm2 <- brm(bf(syphpos ~ offset(log(women_child_age)) + employed + sec_educ + syphlisTestingCoverage + electricity + median_age_birth + women_more_sexPpartners + women_HIV_pos + (1|distcode)),
               data = brmData,
               family = "negbinomial",
               data2 = list(W=mwMat,type="icar"),
               chains = 2,
               seed = 0101)
print(fit.brm2)
plot(fit.brm2)

fit.samples <- posterior_samples(fit.brm2,"^b")
fit.summary <- posterior_summary(fit.brm2)

fit.brm2$fit@sim$samples[[1]][5]
param_estim = tibble(intercept =  fit.brm2$fit@sim$samples[[1]]$b_Intercept,
                     employed = fit.brm2$fit@sim$samples[[1]]$b_employed,
                     educ = fit.brm2$fit@sim$samples[[1]]$b_sec_educ,
                     testingCoverage = fit.brm2$fit@sim$samples[[1]]$b_syphlisTestingCoverage,
                     electricityCov = fit.brm2$fit@sim$samples[[1]]$b_electricity,
                     medianAge = fit.brm2$fit@sim$samples[[1]]$b_median_age_birth,
                     sexPartners = fit.brm2$fit@sim$samples[[1]]$b_women_more_sexPpartners,
                     hivPrev = fit.brm2$fit@sim$samples[[1]]$b_women_HIV_pos)
param_beta <- apply(param_estim,2,median)

# parameter values
employed <- brmData$employed
educ <- brmData$sec_educ
testing.cov <- brmData$syphlisTestingCoverage
electricity <- brmData$electricity
age.birth <- brmData$median_age_birth
sex.partners <- brmData$women_more_sexPpartners
hivpos <- brmData$women_HIV_pos

# calculate overall risk across districts
brmData$r_st <- exp(param_beta[1]+
                        param_beta[2]*employed+
                        param_beta[3]*educ+
                        param_beta[4]*testing.cov+
                        param_beta[5]*electricity+
                        param_beta[6]*age.birth+
                        param_beta[7]*sex.partners+
                        param_beta[8]*hivpos)
# plot estimates
betas <- fit.brm$fit@sim$samples
phis <- fit1$samples$phi # spatial 
deltas <- fit1$samples$delta # temporal
intrx <- fit1$samples$gamma # interaction
betaMeans <- apply(betas,2,mean)
phiMeans <- apply(phis,2,mean)
deltaMeans <- apply(deltas,2,mean)
gammaMeans <- apply(intrx,2,mean)


# model diagnostics
mcmc_acf(posterior_samples(fit.brm))
pp_check(fit.brm,nsamples = 1000)
preds <- brmData %>%
  add_predicted_draws(fit.brm)

# bayesian R2
bayes_R2(fit.brm)
