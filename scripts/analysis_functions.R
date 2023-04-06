# functions for data analysis
# create a function for the figures
# district level prevalence
plot_district_prevalence <- function(summary.district.poly,lake.poly){
  malawi_prev <- tm_shape(summary.district.poly) +
    tm_polygons(col = "syph_prev",border.col = "black",title = "Prevalence", palette = "YlGnBu") +
    tm_facets(by = "year",ncol = 5, nrow = 2) +
    tm_shape(lake.poly) +
    tm_polygons(col = "lightblue", border.col = "black") +
    tm_layout(legend.position = c("left","bottom"),
              legend.outside = FALSE,
              legend.text.size = 1.1,
              legend.title.size = 1.2,
              panel.label.size = 1.3,
              panel.label.fontface = "bold",
              asp = 0)
  return(malawi_prev)
}

# district level smr
plot_district_smr <- function(summary.district.poly,lake.poly){
  malawi_smr <- tm_shape(summary.district.poly) +
    tm_polygons(col = "smr",border.col = "black", title = "SMR", palette = "YlGnBu") +
    tm_facets(by = "year", ncol = 5, nrow = 2) +
    tm_shape(lake.poly) +
    tm_polygons(col = "lightblue", border.col = "black") +
    tm_layout(legend.outside = FALSE,
              legend.position = c("left","bottom"),
              legend.text.size = 1.1,
              legend.title.size = 1.2,
              panel.label.size = 1.3,
              panel.label.fontface = "bold")
  return(malawi_smr)
}

# exceedance probabilities
plot_exceedance_probabilities <- function(summary.district.poly,lake.poly){
  exprob <- tm_shape(summary.district.poly) +
    tm_polygons(border.col = "black",col = "exprob", title = "Exceedance prob.", palette = "YlGnBu") +
    tm_facets(by = "year", nrow = 2, ncol = 5) +
    tm_shape(lake.poly) +
    tm_polygons(col = "lightblue", border.col = "black") +
    tm_layout(legend.outside = FALSE,
              legend.position = c("left","bottom"),
              legend.text.size = 1.1,
              legend.title.size = 1.2,
              panel.label.size = 1.3,
              panel.label.fontface = "bold")
  return(exprob)
}

plot_predicted_irr <- function(summary.district.poly, lake.poly, break.style){
  irr_yr <- tm_shape(summary.district.poly) +
    tm_polygons(border.col = "black",col = "irr", title = "IRR", n = 7, style = break.style, palette = "YlGnBu") +
    tm_facets(by = c("year"), nrow = 2, ncol = 5) +
    tm_shape(lake.poly) +
    tm_polygons(col = "lightblue", border.col = "black") +
    tm_layout(legend.outside = FALSE,
              legend.position = c("left","bottom"),
              legend.text.size = 1.1,
              legend.title.size = 1.2,
              legend.format = list(fun = function(x) formatC(x, digits = 2, format = "f")),
              panel.label.size = 1.3,
              panel.label.fontface = "bold")
}
# function to add denoninator - women of child bearing age
split_data_by_district <- function(syphilis_data,district_code,pop_vec){
  df <- filter(syphilis_data,distcode == district_code)
  df$women_child_age <- rep(pop_vec,each = 12)
  return(df)
}

# function to capitalize first letter of a character
simpleCap <- function(x){
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s,1,1)), substring(s,2),
        sep = "", collapse = " ")
}

# function to plot ecdf
plot_ecdf_diagnostics <- function(modelFit,varpos,midpt,nSample){
  ecdfLower <- ecdf(modelFit$samples$beta[,varpos][1:midpt])
  ecdfUpper <- ecdf(modelFit$samples$beta[,varpos][midpt + 1:nSample])
  plot(ecdfLower,col = "black",lwd = 2)
  lines(ecdfUpper,col = "red",lwd = 2)
  legend("topleft",legend = c("Lower","Upper"),
         lwd = c(2,2),
         pt.cex = 1,
         cex = 2,
         col = c("black","red"),
         bty = "n")
}

# function to create and summarize yearly exceedance probabilities
yearly_posterior_exceedance_prob <- function(year,prob_threshold = 1){
  estim_risk_samples <- risk_samples_combined[,finaldata$year == year]
  risk_estim_year <- apply(estim_risk_samples,2,median)
  pep_estim_year <- apply(estim_risk_samples > prob_threshold,2,mean)
  fit_data <- finaldata[finaldata$year == year,] %>%
    mutate(pep = pep_estim_year,
           r_st = risk_estim_year,
           year = year)
  fit_data_summ <- fit_data %>%
    group_by(district) %>%
    summarise(exprob = mean(pep,na.rm = T),
              risk = mean(r_st,na.rm = T)) %>%
    add_row(district = "likoma",exprob = NA,.after = 9)
  return(fit_data_summ)
}
