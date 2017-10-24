############################################
# Script to calculate CO2 flux
# from Vaisala IRGA, GMP343.
# Gives flux as 'mg C m-2 h-1'
# Read all comments as some things need to
# be adapted for specific system/chamber/etc
#
# Fit linear and quadratic models
#
# Gustaf.Granath@gmail.com
#
# This function DO NOT account 
# for changes in chamber humidity
############################################

# Usage: save Vaisala output as csv, one file for each sample.
#        put this file in the same catalog and run it source("soil_flux_calc.R")

require(ggplot2) # for plotting

files <- list.files(pattern = "*.csv") # list all files in the directory

# Read the files and fix column names
varNames<-c("time","co2","temp") # Fix Vaisala column names
myfiles = lapply(files, read.csv, header=F, sep=",", dec=".",col.names=varNames, skip = 1) 
# NOTE: CHANGE sep and dec depending on csv file type 
# sep = ',' or ';' as separator
# dec = "." or "," as decimal indicator

#myfiles2<-lapply(myfiles,FUN=t)
lis.name <- substr(unlist(files), 1, nchar(unlist(files)) - 4) # store file names
names(myfiles)<- lis.name # name each data file after file name

# time for each sample
dates.tim <- strptime(lis.name, "%Y-%m-%d %H_%M") # extract time from file name

# time for each sample
meas.sec <- lapply(myfiles, function (x) as.numeric(strptime(x[,1], "%Y-%m-%d %H:%M:%S")) - 
                     as.numeric(strptime(x[1,1], "%Y-%m-%d %H:%M:%S"))) # new column with seconds since first measurement

myfiles <- Map(cbind, myfiles, seconds = meas.sec) # add 'seconds' column to each sample

# OPTIONAL
# add treatment, site name, id... or whatever you want
# each vector should have the same length as the number of files
#treat = xx
#site = xx
#id = xx
#myfiles <- Map(cbind, myfiles, treat = treat, site = site, id = id) 

# functions to estiamte slope
# first set the data point range to include in the calculations
int <- 2:12 # here point number 2 to 12
slope.fun.quad <- function (x) {
  if(nrow(x) >1) {ll <-lm(x[int, "co2"]~x[int, "seconds"] + I(x[int, "seconds"]^2) )}
  else {NA} # give NA if there is only one data point
  return(ll)
}

slope.fun.lin <- function (x) {
  if(nrow(x) >1) {ll <-lm(x[int, "co2"]~x[int, "seconds"] )}
  else {NA} # give NA if there is only one data point
  return(ll)
}

# Estimate slopes
  #linear
mods.lin <- lapply(myfiles,  FUN = slope.fun.lin) # fit linear slopes
err.lin <- lapply(mods.lin, function (x) if(!(is.na(x[1]))) {summary(x)$sigma}) #rmse
r2.lin <- lapply(mods.lin, function (x) if(!(is.na(x[1]))) {summary(x)$r.squared})#R2  
coef.mod.lin <- lapply(mods.lin, function (x)  if(!(is.na(x[1]))) {coef(x)[2]}) # extract linear slope
  #quadratic
mods.quad <- lapply(myfiles,  FUN = slope.fun.quad) # fir quadratic model
err.quad <- lapply(mods.quad, function (x) if(!(is.na(x[1]))) {summary(x)$sigma}) #rmse
r2.quad <- lapply(mods.quad, function (x) if(!(is.na(x[1]))) {summary(x)$r.squared})#R2  
coef.mod.quad <- lapply(mods.quad, function (x)  if(!(is.na(x[1]))) {coef(x)[2]}) # extract slope at start of time

airTemp <-lapply(myfiles, function (x) mean(x[int, "temp"], na.rm = TRUE)) # cala mean temp over measuring period
  #h2o <-lapply(myfiles, function (x) mean(x[2:13, "humidity"]))
co2 <-  lapply(myfiles, function (x) x[int[2], "co2"]) # CO2 value at first data point
  #site <-  unlist(lapply(myfiles, function (x) x[1, c("site")]))
  #treat <-  unlist(lapply(myfiles, function (x) x[1, c("treat")]))
  #samp.id <-  unlist(lapply(myfiles, function (x) x[1, c("id")]))
dates.tim  <-  unlist(lapply(myfiles, function (x) x[1, c("time")])) # store time stamp of measurement
  
  # get flux (partly copied from R package 'flux')
  M = 12.01  #for C mass, use 44.01 if you want CO2 
  t = 1/3600 # transform seconds into hours 
  p = 101325 # pressure. Assumed to be the same for all measurements. Not OK for altitudinal gradients!
  R <- 8.314 # gas constant for the used units
    # chamber size
  area <- pi*0.125^2 # area of chamber, m2. Here the radius is 12.5 cm 
  v <- 0.12*area # chamber volume, m3. Here the height is 12 cm
  v <- v - (0.05 * pi*0.02^2) # IMPORTANT: the volume of the CO2 unit must be accounted for!!
                        # you need to measure this, numbers here are just examples
  
  flux.lin <- ((unlist(coef.mod.lin) * v * M * p) / (t * R * (unlist(airTemp) + 273.15) * area)) #micro gram C m-2 h-1
  flux.lin <- flux.lin/1e+3 # mg C m-2 h-1
  flux.quad <- ((unlist(coef.mod.quad) * v * M * p) / (t * R * (unlist(airTemp) + 273.15) * area)) #micro gram C m-2 h-1
  flux.quad <- flux.quad/1e+3 # mg C m-2 h-1

# save data
  names(flux.lin)<-NULL; names(flux.quad)<-NULL
  names(err.lin)<-NULL; names(err.quad)<-NULL
  names(airTemp)<-NULL; names(co2)<-NULL 
  
  flux.dat <- data.frame("time"=dates.tim,  
                      "flux.linear"=flux.lin, "airTemp"=unlist(airTemp), "co2_at_start" = unlist(co2),
                      "rmse.linear"=unlist(err.lin),"R2.lin"=unlist(r2.lin), #site = site, "treat" = treat, samp.id = samp.id, 
                      "flux.quad"=flux.quad, unlist(err.quad),"R2.lin"=unlist(r2.quad))

# Plot all samples and save as pdf
# Optional to run
myfiles.plot <- Map(cbind, myfiles, sample = names(myfiles)) # add seconds column to each sample
myfiles.plot <- do.call(rbind.data.frame, myfiles.plot)
plots <- list()
for (i in 1:length(unique(myfiles.plot$sample))) {
  id <- unique(myfiles.plot$sample)[i]
  samp.dat <- myfiles.plot[myfiles.plot$sample == id,]
  lin.n <- round(unlist(coef.mod.lin[names(coef.mod.lin[]) == id]),3)
  lin <- paste("linear=",lin.n,sep="")
  quad.n <- round(unlist(coef.mod.quad[names(coef.mod.quad[]) == id]),3)
  quad <- paste("quad=",quad.n,sep="")
  
  plots[[i]] <- ggplot(samp.dat, aes(x=seconds, y=co2)) +
   geom_point() +
  geom_smooth(data = samp.dat[int,], method = "lm", se = FALSE, color = "blue") +
  geom_smooth(data = samp.dat[int,], method = "lm",formula = y ~ x + I(x^2), se = FALSE, color ="red") +
  geom_text(data = samp.dat[1,], aes(x = seconds+100, y = co2+40), label = lin,size=8, color = "blue" ) +
  geom_text(data = samp.dat[1,], aes(x = seconds+100, y = co2+20), label = quad,size=8, color = "red" )
  
  # Run if you have eg site id and treatments and want label the plot with such info
  #   + geom_text(data = myfiles.plot[!(duplicated(myfiles.plot$sample)), c("site", "id", "time", "sample")] ,
  #       aes(x=3,y=500,label=paste(site,id, sep="-")), color = "red")
}

# save all plots. One per page
pdf("co2VStime.pdf",width = 8, height=5)
plots
dev.off()
