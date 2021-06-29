rm(list = ls())

library(rgdal)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(maptools)
library(spdep)
library(INLA)
library(CARBayesST)
library(CARBayes)
library(MASS)
library(brms)
library(bayesplot)
library(tidybayes)
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
denom <- read.csv("data/epidata/denominator_data.csv",header = T)
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

# function to add denoninator - women of child bearing age
split_data_by_district <- function(syphilis_data,district_code,pop_vec){
  df <- filter(syphilis_data,distcode == district_code)
  df$women_child_age <- rep(pop_vec,each=12)
  return(df)
}

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
tiff("images/malawi_syphilisPos.tif",width = (25*0.39),height = (35*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(2,2,2,2))
mwdistr$syph <- casesDistr$prevdistr
brks <- c(0,5,10,15,20,25)
cols <- terrain.colors(6)
plot(mwdistr)
plot(lakes,add=TRUE,col="lightblue")
plot(mwdistr,col=cols[findInterval(mwdistr$syph,brks)],add=TRUE)
legend("bottomleft",
       legend = leglabs(brks,"<",">="),
       fill = cols,
       bty = "n",
       pt.cex = 1,
       cex = 2)
box(lwd=1,bty="o")
dev.off()


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
Plot.district.prevelence <- function(data,timeLen,districtCode,plotLabel,ymark){
  newdata <- dat[dat$distcode==districtCode,]
  newdataSeries <- rep(NA,timeLen)
  for(i in 1:timeLen){
    newdataSeries[i] <- sum(newdata$syphpos[newdata$month==i],na.rm = T)/(sum(newdata$allwomen[newdata$month==i],na.rm=T))*100 # edit demoninator
  }
  par(cex.lab=1.5,cex.axis=1.5,mar=c(4,4.5,2,2))
  plot(1:timeLen,newdataSeries,type='l',col="black",lwd=1.5,xlab="Year",ylab="Prevalence (%)", ylim=c(0,8),axes=F)
  box(lwd=1,bty="o")
  axis(1,at=seq(1,ntime,12),labels=2014:2020,lwd=1)
  axis(2,lwd = 1)
  mtext(plotLabel,side = 2,line = 1,cex = 1.7,at=ymark,0,las=2)
}

# function to capitalize first letter of a character
simpleCap <- function(x){
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s,2),
        sep="", collapse=" ")
}


## plot time series for selected districts
tiff("images/ditrict_level_prevelence.tif",width = (45*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(3,5),mar=c(1.3,1.3,1.3,1.3))
nd <- unique(dat$district)[1:15]
for(i in 1:15){
  Plot.district.prevelence(dat,72,i,"",25)
  title(simpleCap(nd[i]))
}
dev.off()

## calculate the expected number of cases
dat$overallProb <- sum(dat$syphpos,na.rm = T)/sum(dat$allwomen,na.rm = T) # change denominator
dat$expectedCases <- ceiling(dat$overallProb * dat$allwomen)
dat$smr <- dat$syphpos/dat$expectedCases

# summarize and plot smr
smrData <- dat %>%
  group_by(year) %>%
  summarise(smrYear = mean(smr,na.r=T))
plot(smrData$year,smrData$smrYear,type="l",xlab="",ylab="SMR",lwd=2)

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
par(mfrow=c(1,1),mar=c(2,2,2,2))
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
       cex = 2)
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
##### function uses the data from the data frame G summarized above 
plot.syphilis.yearly.prevelence <- function(year){
  if(year==2014){
    x <- "prev2014"
  } else if (year==2015){
    x <- "prev2015"
  } else if (year==2016){
    x <- "prev2016"
  } else if (year==2017){
    x <- "prev2017"
  } else if (year==2018){
    x <- "prev2018"
  } else if (year==2019){
    x <- "prev2019"
  } else {
    x <- "prev2020"
  }
  brks <- c(0.5,1,1.5,2,2.5,3,4,4.5,5,5.5)
  #mycol <- colorRampPalette(c("lightgreen","yellow","gold"))(length(brks))
  mycol <- colorRampPalette(c("#f5eef8","#6c3483"))(length(brks))
  plot(mwdistr)
  plot(lakes,add=T,col="lightblue")
  plot(mwdistr,col=mycol[findInterval(mwdistr@data[,x],brks)],add=T)
  legend("bottomleft",
         legend = leglabs(brks,"<",">="),
         fill = mycol,
         bty = "n",
         pt.cex = 1,
         cex = 1)
  box(lwd=1,bty="o")
}

# plot yearly prevalence
tiff("images/prevelence_year.tif",width = 35*0.39,height = 25*0.39,units = "in",res = 450,compression = "lzw")
par(mfrow=c(2,4),mar=c(2,2,2,2))
for(i in 2014:2020){
  plot.syphilis.yearly.prevelence(i)
  title(main=i)
}
dev.off()

# plot yearly smr
tiff("images/yearly_SMR.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 450,compression = "lzw")
par(mfrow=c(2,4),mar=c(2,2,2,2))
for(i in 2014:2020){
  plot.syphilis.yearly.smr(year=i)
  title(main = i)
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


plot.syphilis.testing.coverage <- function(year){
  if(year==2014){
    x <- "y1cov"
  } else if (year==2015){
    x <- "y2cov"
  } else if (year==2016){
    x <- "y3cov"
  } else if (year==2017){
    x <- "y4cov"
  } else if (year==2018){
    x <- "y5cov"
  } else if (year==2019){
    x <- "y6cov"
  } else {
    x <- "y7cov"
  }
  brks <- c(40,50,60,70,80,90)
  cols <- colorRampPalette(rev(c("gold","yellow","green")))(9)
  plot(mwdistr)
  plot(lakes,add=T,col="lightblue")
  plot(mwdistr,col=cols[findInterval(mwdistr@data[,x],brks)],add=T)
  legend("bottomleft",
         legend = leglabs(brks,"<",">="),
         fill = cols,
         bty = "n",
         pt.cex = 1,
         cex = 1)
  box(lwd=1,bty="o")
}
tiff("images/testing_coverages_years.tif",width = 35*0.39,height = 25*0.29,units = "in",res = 350,compression = "lzw")
par(mfrow=c(2,4),mar=c(2,2,2,2))
for(i in c(2014:2020)){
  plot.syphilis.testing.coverage(year = i)
  title(main=i)
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
m1 <- glm(syphpos ~ offset(log(women_child_age)) +
            sec_educ + literacy + employed + perc_live_birth + median_age_sex+
            median_age_birth +women_more_sexPpartners+polygamy + men_condom + paid_sex + women_HIV_pos+
            STI + improved_sanit + improved_water + electricity + skilled_ANC +no_ANC + problem_health_care + nopostnatal_check + facility_deliv,
          data = dat,family = quasipoisson())
summary(m1)
exp(cbind(IRR=coef(m1),CI=confint(m1)))

m2 <- glm(syphpos ~ offset(log(women_child_age)) + employed+
            sec_educ + syphlisTestingCoverage  + electricity +
            median_age_birth +women_more_sexPpartners + women_HIV_pos,
          data = dat,family = poisson())
summary(m2)
exp(cbind(IRR=coef(m2),CI=confint(m2)))


m3 <- glm(syphpos ~ offset(log(women_child_age)) + sec_educ + testCov + no_ANC,
          data = dat,family = quasipoisson())
summary(m3)
exp(cbind(IRR=coef(m3),CI=confint(m3)))

# reduced preliminary model with only key risk factors
m4 <- glm(syphpos ~ offset(log(women_child_age)) + sec_educ + testCov + no_ANC + women_more_sexPpartners,problem_health_care + women_HIV_pos,
          data = dat,family = quasipoisson())
summary(m4)
exp(cbind(IRR=coef(m4),CI=confint(m4)))


### regression modelling for lattice data n- spatio-temporal modelling ### 
### CAR models 

### CAR models using CARBayesST
MWnb <- poly2nb(mwdistr)
mwMat <- nb2mat(MWnb,style = "B", zero.policy = TRUE)
mwMat <- mwMat[-10,-10] # remove Likoma 
finaldata <- dat
finaldata <- filter(dat,district != "likoma")

## remove missing values

finaldata$syphpos[is.na(finaldata$syphpos)] <- floor(mean(finaldata$syphpos,na.rm=TRUE))
finaldata$expectedCases[is.na(finaldata$expectedCases)] <- mean(finaldata$expectedCases,na.rm=TRUE)
finaldata$TFR[is.na(finaldata$TFR)] <- mean(finaldata$TFR,na.rm=TRUE)
finaldata$syphlisTestingCoverage[is.na(finaldata$syphlisTestingCoverage)] <- mean(finaldata$syphlisTestingCoverage,na.rm=TRUE)
# process additional covariates

formula <- syphpos ~ offset(log(women_child_age)) + employed+ sec_educ + syphlisTestingCoverage  + electricity +
  median_age_birth +women_more_sexPpartners + women_HIV_pos


set.seed(0012)
fit1 <- ST.CARanova(formula = formula,
                  family = "poisson",
                  data = finaldata, 
                  W=mwMat,
                  burnin = 3000,
                  n.sample = 500000,
                  thin = 50,
                  verbose = TRUE)
#summary(fit1)
print(fit1)
colnames(fit1$samples$beta) <- c("Intercept","employed","Education","Testing coverage","electricity","Age birth","Sex partners","HIV positive")
plot(exp(fit1$samples$beta[,-1]))

# model diagnostics
plot.ecdf(ecdf(fit1$samples$beta[,2][1:2950])) # border
lines(ecdf(fit1$samples$beta[,2][2951:5900]))

# function to plot ecdf
plot_ecdf_diagnostics <- function(modelFit,varpos,midpt,nSample){
  ecdfLower <- ecdf(modelFit$samples$beta[,varpos][1:midpt])
  ecdfUpper <- ecdf(modelFit$samples$beta[,varpos][midpt+1:nSample])
  plot(ecdfLower,col="black",lwd=2)
  lines(ecdfUpper,col="red",lwd=2)
  legend("topleft",legend = c("Lower","Upper"),
         lwd = c(2,2),
         pt.cex = 1,
         cex = 2,
         col = c("black","red"),
         bty = "n")
}
plot_ecdf_diagnostics(modelFit = fit1,varpos = 2,midpt = 250000,nSample = 500000)

# incident rate ratios
params <- summarise.samples(exp(fit1$samples$beta[,-1]),quantiles = c(0.5,0.025,0.975))
round(params$quantiles,2)

# map fitted values
betas <- fit1$samples$beta
phis <- fit1$samples$phi # spatial 
deltas <- fit1$samples$delta # temporal
intrx <- fit1$samples$gamma # interaction
betaMeans <- apply(betas,2,mean)
phiMeans <- apply(phis,2,mean)
deltaMeans <- apply(deltas,2,mean)
gammaMeans <- apply(intrx,2,mean)

# lower quantiles
betaLower <- apply(betas,2,function(x)quantile(x,0.025))
phiLower <- apply(phis,2,function(x)quantile(x,0.025))
deltaLower <- apply(deltas,2,function(x)quantile(x,0.025))
gammaLower <- apply(intrx,2,function(x)quantile(x,0.025))

# upper quantiles
betaUpper <- apply(betas,2,function(x)quantile(x,0.975))
phiUpper <- apply(phis,2,function(x)quantile(x,0.975))
deltaUpper <- apply(deltas,2,function(x)quantile(x,0.975))
gammaUpper <- apply(intrx,2,function(x)quantile(x,0.975))

# covariate values
employed <- finaldata$employed
educ <- finaldata$sec_educ
testing.cov <- finaldata$syphlisTestingCoverage
electricity <- finaldata$electricity
age.birth <- finaldata$median_age_birth
sex.partners <- finaldata$women_more_sexPpartners
hivpos <- finaldata$women_HIV_pos

# calculate overall risk across districts
finaldata$r_st <- exp(betaMeans[1]+
              betaMeans[2]*employed+
              betaMeans[3]*educ+
              betaMeans[4]*testing.cov+
              betaMeans[5]*electricity+
              betaMeans[6]*age.birth+
              betaMeans[7]*sex.partners+
              betaMeans[8]*hivpos+
              phiMeans+deltaMeans+gammaMeans)
# lower quantile risk
finaldata$r_st_lower <- exp(betaLower[1]+
                              betaLower[2]*employed+
                              betaLower[3]*educ+
                              betaLower[4]*testing.cov+
                              betaLower[5]*electricity+
                              betaLower[6]*age.birth+
                              betaLower[7]*sex.partners+
                              betaLower[8]*hivpos+
                              phiMeans+deltaMeans+gammaMeans)
# upper quantile risk
finaldata$r_st_upper <- exp(betaUpper[1]+
                              betaUpper[2]*employed+
                              betaUpper[3]*educ+
                              betaUpper[4]*testing.cov+
                              betaUpper[5]*electricity+
                              betaUpper[6]*age.birth+
                              betaUpper[7]*sex.partners+
                              betaUpper[8]*hivpos+
                              phiMeans+deltaMeans+gammaMeans)

plotData <- finaldata %>%
  group_by(district) %>%
  summarise(avgIRR = mean(r_st,na.rm=TRUE)) %>%
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
  summarise(irr_2014 = mean(r_st[year==2014],na.rm=TRUE),
            irr_2015 = mean(r_st[year==2015],na.rm=TRUE),
            irr_2016 = mean(r_st[year==2016],na.rm=TRUE),
            irr_2017 = mean(r_st[year==2017],na.rm=TRUE),
            irr_2018 = mean(r_st[year==2018],na.rm=TRUE),
            irr_2019 = mean(r_st[year==2019],na.rm=TRUE),
            irr_2020 = mean(r_st[year==2020],na.rm=TRUE)) %>%
  add_row(district="likoma",irr_2014=NA,irr_2015=NA,irr_2016=NA,irr_2017=NA,irr_2018=NA,irr_2019=NA,irr_2020=NA,.after = 9)
## add yearly irr to spatial dataframe
mwdistr$irr_2014 <- estimYear$irr_2014
mwdistr$irr_2015 <- estimYear$irr_2015
mwdistr$irr_2016 <- estimYear$irr_2016
mwdistr$irr_2017 <- estimYear$irr_2017
mwdistr$irr_2018 <- estimYear$irr_2018
mwdistr$irr_2019 <- estimYear$irr_2019
mwdistr$irr_2020 <- estimYear$irr_2020

tiff("images/yealy_rate_ratios_april.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 450,compression = "lzw")
par(mfrow=c(3,3),mar=c(2,2,2,2))
for(i in 2014:2020){
  plot.yearly.incident.rates(i)
  title(main = i)
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

# create a matrix
irrMat <- matrix(nrow = nyear,ncol = nmonth)
for(i in 1:nmonth){
  for(j in 1:nyear){
    irrMat[j,i] <- mean(finaldata$r_st[finaldata$year==j | finaldata$month==i])
  }
}

filled.contour(x = temporalTrendData$month,y = temporalTrendData$year ,z = temporalTrendData$mean_irr)

x <- 1:dim(M)[1]
y <- 1:ncol(M)
filled.contour(x,y,M)

## summary for fit2
print(fit2)
colnames(fit2$samples$beta) <- c("Intercept","border","Education","HIMS coverage","TFR","ANC","Sex partners","HC problem","HIV")
plot(exp(fit2$samples$beta[,-1]))

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

plot(temporalTrendData$month,temporalTrendData$irrMean,type="l",ylim=c(0,4),lwd=2,axes=F)
axis(1,at=seq(1,84,12),labels = c(2014:2020))
axis(2,lwd = 1)
lines(temporalTrendData$month,temporalTrendData$L_IRR,type = "l",lty=2)
lines(temporalTrendData$month,temporalTrendData$U_IRR,type = "l",lty=2)
abline(h=1,lty=2)
box(lwd=1)

# predicted estimates of relative risk

# model using brms
not_1_car <- brm(bf(total_notif ~ 1 + 
                      offset(log(population))),
                 data=dat_scale_ln, 
                 family=poisson,
                 control = list(adapt_delta = 0.99, max_treedepth=10),
                 autocor=cor_car(w4, ~ 1 | scale_cluster_area, type = "icar"),
                 inits = 0,
                 #prior = prior_not_1,
                 cores=3,
                 iter=15000, warmup=1000,
                 seed = 1237,
                 chains=3)

brmData <- dat %>% filter(district != "likoma")
fit.brm <- brm(bf(syphpos ~ offset(log(women_child_age)) + employed + sec_educ + syphlisTestingCoverage + electricity + median_age_birth + women_more_sexPpartners + women_HIV_pos),
               data = brmData,
               family = "negbinomial",
               data2 = list(W=mwMat,type="icar"),
               chains = 2,
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
