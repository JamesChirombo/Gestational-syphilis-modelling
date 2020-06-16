rm(list = ls())

library(rgdal)
library(tidyverse)
library(RColorBrewer)
library(maptools)

## load required spatial data
## malawi districts

maldistricts <- readOGR("data/spatial","mwi_admbnda_adm2")
districts <- readOGR("data/spatial","mw_districts")
proj4string(districts) <- CRS("+init=epsg:32736")
mwdistr <- spTransform(districts,CRS("+init=epsg:4326"))
lakes <- readOGR("data/spatial","lakes")


## load syhphilis data

dat <- read.csv("data/epidata/syphilis_data.csv",header = T,stringsAsFactors = F)

## compute the total number of women

dat$allwomen <- rowSums(dat[,c("anc_vst1","anc_vst2","anc_vst3","anc_vst4","anc_vst5")])
dat$prev <- dat$syphpos/dat$allwomen*1000
# exploratory analysis
# map raw numbers of positive cases per district
ntime <- length(unique(dat$month))

# calculate syphilis prevalence among women cohort

casesDistr <- dat %>%
  group_by(district) %>%
  summarise(totalCases = sum(syphpos,na.rm = TRUE),
            prevdistr = mean(prev,na.rm = TRUE))

## map the total nuumber of cases per district
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

## plot overal changes in cases over time
casesTime <- dat %>%
  group_by(month) %>%
  summarise(syphprev = mean(dat,na.rm = T))

### calculate average incidence over time

## yearly syphilis prevalence

# Average monthly prevalence -- national level estimates
nx <- unique(dat$month)
syphprev <- c()
for( i in 1:length(nx)){
  syphprev[i] <- sum(dat$syphpos[dat$month==i],na.rm = T)/sum(dat$allwomen[dat$month==i],na.rm=T)*100
}

tiff("images/national_prevelence.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(4.5,4.5,2,2),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,syphprev,type="l",lwd=1.5,axes=F,xlab="Year",ylab="Prevalence (%)")
axis(1,at=seq(1,ntime,12),labels = 2014:2019)
axis(2,lwd = 1)
box(lwd=1,bty="o")
dev.off()

### syphilis prevalence at the district level ###
Plot.district.prevelence <- function(data,timeLen,districtCode,plotLabel,ymark){
  newdata <- dat[dat$distcode==districtCode,]
  newdataSeries <- rep(NA,timeLen)
  for( i in 1:timeLen ){
    newdataSeries[i] <- sum(newdata$syphpos[newdata$month==i],na.rm = T)/(sum(newdata$allwomen[newdata$month==i],na.rm=T))*100
  }
  par(cex.lab=1.5,cex.axis=1.5,mar=c(4,4.5,2,2))
  plot(1:timeLen,newdataSeries,type='l',col="black",lwd=1.5,xlab="Year",ylab="Prevalence (%)", ylim=c(0,8),axes=F)
  box(lwd=1,bty="o")
  axis(1,at=seq(1,ntime,12),labels=2014:2019,lwd=1)
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


### map syphilis prevelence among the women cohort ####


G <- dat %>%
  group_by(district) %>%
  summarise(y1 = sum(syphpos[year==2014],na.rm = T)/sum(allwomen[year==2014],na.rm = T)*100,
            y2 = sum(syphpos[year==2015],na.rm = T)/sum(allwomen[year==2015],na.rm = T)*100,
            y3 = sum(syphpos[year==2016],na.rm = T)/sum(allwomen[year==2016],na.rm = T)*100,
            y4 = sum(syphpos[year==2017],na.rm = T)/sum(allwomen[year==2017],na.rm = T)*100,
            y5 = sum(syphpos[year==2018],na.rm = T)/sum(allwomen[year==2018],na.rm = T)*100,
            y6 = sum(syphpos[year==2019],na.rm = T)/sum(allwomen[year==2019],na.rm = T)*100)
## map prevelence by distric and time 

mwdistr$prev2014 <- G$y1
mwdistr$prev2015 <- G$y2
mwdistr$prev2016 <- G$y3
mwdistr$prev2017 <- G$y4
mwdistr$prev2018 <- G$y5
mwdistr$prev2019 <- G$y6



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
  } else {
    x <- "prev2019"
  }
  brks <- c(0.5,1,1.5,2,2.5,3,4)
  mycol <- colorRampPalette(c("lightgreen","yellow","gold"))(length(brks))
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

tiff("images/prevelence_year.tif",width = 35*0.39,height = 25*0.39,units = "in",res = 350,compression = "lzw")
par(mfrow=c(2,3),mar=c(2,2,2,2))
for(i in 2014:2019){
  plot.syphilis.yearly.prevelence(i)
  title(main=i)
}
dev.off()

## compute coverages ##
testcov <- c()
for(i in 1:ntime){
  testcov[i] <- 100-sum(dat$syphunk[dat$month==i],na.rm = T)/sum(dat$allwomen[dat$month==i],na.rm = T)*100
}
tiff("images/testingCoverage.tif",width = (35*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(4.5,4.5,2,2),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,testcov,type="l",lwd=1.5,axes=F,xlab="Year",ylab="Coverage of testing (%)")
axis(1,at=seq(1,ntime,12),labels = 2014:2019)
axis(2,lwd = 1)
box(lwd=1,bty="o")
dev.off()

## summarize coverage by year
U <- dat %>%
  group_by(district) %>%
  summarise(y1 = sum(syphunk[year==2014],na.rm = T)/sum(allwomen[year==2014],na.rm = T),
            y2 = sum(syphunk[year==2015],na.rm = T)/sum(allwomen[year==2015],na.rm = T),
            y3 = sum(syphunk[year==2016],na.rm = T)/sum(allwomen[year==2016],na.rm = T),
            y4 = sum(syphunk[year==2017],na.rm = T)/sum(allwomen[year==2017],na.rm = T),
            y5 = sum(syphunk[year==2018],na.rm = T)/sum(allwomen[year==2018],na.rm = T),
            y6 = sum(syphunk[year==2019],na.rm = T)/sum(allwomen[year==2019],na.rm = T))
## plot testing coverage at different time points

mwdistr$y1cov <- (1-U$y1)*100
mwdistr$y2cov <- (1-U$y2)*100
mwdistr$y3cov <- (1-U$y3)*100
mwdistr$y4cov <- (1-U$y4)*100
mwdistr$y5cov <- (1-U$y5)*100
mwdistr$y6cov <- (1-U$y6)*100


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
  } else {
    x <- "y6cov"
  }
  brks <- c(10,20,30,40,50,60,70,80,90)
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
par(mfrow=c(2,3),mar=c(2,2,2,2))
for(i in c(2014:2019)){
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
plot(1:ntime,syphprev,type="l",lwd=2,axes=F,xlab="Year",ylab="Prevelence (%)",col="red")
axis(1,at=seq(1,ntime,12),labels = 2014:2019)
axis(2,lwd = 1)
par(new=T)
plot(1:ntime,newbornUwt,type="l",lwd=2,xlab=NA,ylab=NA,axes=F,col="blue")
axis(4,at=seq(1,10,1),lwd = 1.5)
mtext("Underweight (%)",side = 4,line = 2.8,cex = 1.5)
box(lwd=1.5,bty="o")
legend("topleft",
       legend = c("Syphilis","Underweight"),
       lty = c(1,1),
       lwd = c(2,2),
       col = c("red","blue"),
       bty = "n",
       pt.cex = 1,
       cex = 2)
dev.off()

### syphilis and prematurity
tiff("images/syphilis_prematurity.tif",width = (30*0.39),height = (20*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),mar=c(4.5,4.5,2,4.5),cex.lab=1.5,cex.axis=1.5)
plot(1:ntime,syphprev,type="l",lwd=2,axes=F,xlab="Year",ylab="Prevelence (%)",col="#ae8419")
axis(1,at=seq(1,ntime,12),labels = 2014:2019)
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
axis(1,at=seq(1,ntime,12),labels = 2014:2019)
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
dat$borderDistr <- ifelse(dat$district %in% c("ntcheu","dedza","mulanje","mzimba","kasungu","machinga","mchinji","nsanje","phalombe","chikwawa","mangochi","chitipa"),1,0)
# add population density

m1 <- glm(syphpos ~ offset(log(allwomen))+covrge + borderDistr,data = dat,family = poisson())
summary(m1)
exp(cbind(IRR=coef(m1),CI=confint(m1)))




