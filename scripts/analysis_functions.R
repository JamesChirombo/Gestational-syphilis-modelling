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

# syphilis testing coverage
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

## function to plot yearly smr rates
plot.syphilis.yearly.smr <- function(year){
  if(year==2014){
    x <- "smr2014"
  } else if (year==2015){
    x <- "smr2015"
  } else if (year==2016){
    x <- "smr2016"
  } else if (year==2017){
    x <- "smr2017"
  } else if (year==2018){
    x <- "smr2018"
  } else if (year==2019){
    x <- "smr2019"
  } else {
    x <- "smr2020"
  }
  brks <- c(0,0.5,1,1.5,2,2.5,3)
  mycol <- colorRampPalette(c("#d6eaf8","#d2b4de","#7d3c98"))(length(brks))
  plot(mwdistr)
  plot(lakes,add=T,col="lightblue",density=50)
  plot(mwdistr,col=mycol[findInterval(mwdistr@data[,x],brks)],add=T)
  legend("bottomleft",
         legend = leglabs(brks,"<",">="),
         fill = mycol,
         bty = "n",
         pt.cex = 1,
         cex = 1)
  box(lwd=1,bty="o")
}

plot.yearly.incident.rates <- function(year){
  if(year==2014){
    x <- "irr_2014"
  } else if (year==2015){
    x <- "irr_2015"
  } else if (year==2016){
    x <- "irr_2016"
  } else if (year==2017){
    x <- "irr_2017"
  } else if (year==2018){
    x <- "irr_2018"
  } else if (year==2019){
    x <- "irr_2019"
  } else {
    x <- "irr_2020"
  }
  #brks <- c(0,0.3,0.6,0.9,1.2,1.5,1.8,2.1)
  #brks <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.4)
  brks <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)
  #mycol <- colorRampPalette(c("lightgreen","yellow","gold"))(length(brks))
  mycol <- colorRampPalette(c("#d6eaf8","#d2b4de","#7d3c98"))(length(brks))
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
