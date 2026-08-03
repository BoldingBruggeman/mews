# Compare model output with observations from monitoring stations.

#Initialization and File Paths
project <- "kinneret"

# get time zone
ltz <- Sys.timezone()

vars <- c("temp") # variables names in GETM output. This does not include depth (zct) which is always needed.
vars0 <- c("zct",vars) # add depth

#stn <- c('A','G','K')
depth_layers <- 20
mo <- c("model","obs")

# input and output pathes
folder0 <- ifelse (Sys.info()[6] == 'shaja', "C:/Users/shaja/OneDrive - IOLR/MEWS", "C:/Users/mestr/OneDrive - IOLR/MEWS") 
path_Bath <- paste0(folder0,"/Git_MEWS/mews/",project,"/Bathymetry/bathymetry.nc")
path_MO <- paste0(folder0,"/Output_example/",project)
path_Out <- paste0(folder0,"/Output_analysis/")
path_Obs <- paste0(folder0,"/Observations/")

library(pygetmtools)
#library(ncdf4)
#library(lubridate)
#library(abind)


#Read Model Data  - adjust to your output structure - many files in one folder or many folders
MOfiles <- list.files(path=path_MO, pattern=paste0(project,"_3d.nc"), full.names=TRUE, recursive = TRUE)
stdt_files <- list.files(path=path_MO)

Obs <- paste0(path_Obs,"Obs1.csv")


cm2o <- pair_sim_and_obs(MOfiles,Obs, list_vars = list(vars), col_date = "Date", col_depth = "Depth", col_x = "lon",col_y = "lat")










# read last file to get start date, end date, time interval. This assumes that the units of the time dimension contains the start date.
nc <- nc_open(MOfiles[length(MOfiles)])
stdt <- as.POSIXct(substr(nc$dim$time$units,15,33),tz=ltz)
time0 <- ncvar_get(nc, 'time') 
out_time_step <- time0[2]-time0[1]
endt <- stdt + time0[length(time0)]
nc_close(nc)
### optional to set start and end date manually: ###
#stdt <- as.POSIXct("2022-01-01 00:00:00", tz=ltz)
#endt <- as.POSIXct("2023-01-01 00:00:00", tz=ltz)
#out_time_step <- 60*60*6 # 6 hours in seconds
  
# get lat lon of center grids from bathymetry
nc <- nc_open(path_Bath)
yLat <- data.frame(lat= ncvar_get(nc, 'latitude')[1,],Y=ncvar_get(nc, 'Y'))
xLon <- data.frame(lon= ncvar_get(nc, 'longitude')[,1],X=ncvar_get(nc, 'X'))
nc_close(nc)

# get monitoring stations coordinates
station <- read.csv(paste0(path_Obs,"station_lat_lon.csv"))

Obs1 <- merge(Obs,station, by.x="Station", by.y="stn")
Obs1$Depth <- (-1)*Obs1$Depth
write.csv(Obs1, paste0(path_Obs,"Obs1.csv"), row.names = FALSE)

# find nearest grids to monitoring stations 
station$Xi <- sapply(station$lon, clv, xLon$lon)
station$Yi <- sapply(station$lat, clv, yLat$lat)

# load observations - should contain: station, depth, time, var1, var2 etc.. and use var names same as in GETM
Obs <- read.csv(paste0(path_Obs,"interpTemp_allSt_2015_2023.csv"))
colnames(Obs1)[5] <- "temp"
# create time vector of unique date-time observations
Obs$DateTime <- ymd_hms(paste0(Obs$Date," ",Obs$time,":00"),tz=ltz)
# reduce to start to end date
Obs <- Obs[Obs$DateTime <= endt & Obs$DateTime >= stdt,]
# Get unique values to build the array
stn <- unique(Obs$Station)
obs_time_vector <- unique(Obs$DateTime)
# create time vector of output. 
out_time_vector <- seq(stdt,endt,by=out_time_step)

# find nearest date-time to observations, write time index.
Obs$ti <- sapply(Obs$DateTime,clv,out_time_vector)

# indexes of needed dates from the model output that correspond to observation dates.
ti1 <- as.integer(unique(Obs$ti))

# prepare array for model results for the whole periods.
ar1 <- array(data=NA,dim=c(length(mo),length(vars0), length(stn),depth_layers,length(out_time_vector[ti1])),
             dimnames = list(mo,vars0,stn,(0:(depth_layers-1)),ti1))  #time dim-name to be date-time: out_time_vector[ti1]

time_running_index <- 1
# get all data from output files and immediately reduce to observed dates and coordinates to reduce memory usage
for (i in MOfiles) {
  #i <- MOfiles[2] # for testing limit to 1 file
  nc <- nc_open(i)
  time0 <- ncvar_get(nc, 'time') 
  
  # prepare array for model output
  ar0 <- array(data=NA,dim=c(length(vars0), nrow(xLon),nrow(yLat),depth_layers,length(time0)),
               dimnames = list(vars0,xLon$X,yLat$Y,(0:(depth_layers-1)),time0))
  
  # get all variables into an array
  for (v in vars0) {
    ar0[v,,,,] <- ncvar_get(nc, v) 
  }
  nc_close(nc)
  
  # turn time to date-time format
  time0 <- as.POSIXct(time0, origin=stdt, format="%Y-%m-%d %H:%M:%S")
  # reduce dates to this file (tv=time vector, 10800 seconds to adjust time difference - should change to using time zone?)
  obs_tv <- data.frame(obs_time = unique(Obs$DateTime[Obs$DateTime <= time0[nrow(time0)] & Obs$DateTime >= time0[1]]), ti=NA)
  # find nearest date-time to observations, write time index.
  obs_tv$ti <- sapply(obs_tv$obs_time,clv,time0)
  ti0 <- unique(obs_tv$ti) # indexes of time
  #time indexes of this file to put in the array (that collects results from all files)
  t_fr <- time_running_index
  t_to <- t_fr + length(ti0) -1
  # now we have indexes of time , X and Y so we can minimize the array to a minimum
  for (v in vars0) {
    for (k in seq_along(stn)) {
      ar1["model",v,stn[k],,t_fr:t_to] <- ar0[v,station$Xi[k],station$Yi[k],,ti0] 
    }
  }
  
  time_running_index <- t_to +1
}

# turn depth to positive
ar1["model","zct",,,] <- (-1)*ar1["model","zct",,,]

# now we have array of model results with closest dates and grid points for all model depth layers.
# I choose here to keep the model layer structure and find for the observation depth, the closest model layer.
# Usually there will be more than one observation per model layer so these are averaged. (the nearest depth is also possible but a bit more complex)
# the observation values are written into the model results array in the 'obs' dimension.

#  add layer no. to observation file
# Order obs so clv function can be applied correctly.
Obs <- Obs[order(Obs$ti,Obs$Station,Obs$Depth),]
# find layer number in the model which is closest to observation depth - for each station and date-time. 
# **** Is there a more efficient way than loop? ***
for (i in stn){
  ti2 <- unique(Obs$ti[Obs$Station==i])
  for (k in ti2) {
    Obs$Di[Obs$Station == i & Obs$ti ==k] <- sapply(Obs$Depth[Obs$Station == i & Obs$ti ==k], clv, ar1["model","zct",i,,as.character(k)])
  }
}

# remove NA's
Obs <- na.omit(Obs)
# aggregate mean per depth layer
Obs1 <- aggregate(list(Obs[,vars]), by=list(Station=Obs$Station,Depth_layer=Obs$Di,time_index=Obs$ti), FUN=mean)
colnames(Obs1)[4:ncol(Obs1)] <- vars
# write into the model array. Check that vars are in the right order if more than 1 var
# *** Another option to avoid loop is to create all combinations of station-depth-time and create a full array ***
for (i in 1:nrow(Obs1)){
  ar1["obs",vars,Obs1$Station[i],Obs1$Depth_layer[i],as.character(Obs1$time_index[i])] <- Obs1[i,vars]
}


# set names for time dimension instead of index
dimnames(ar1)[[5]] <- as.character(out_time_vector[as.integer(dimnames(ar1)[[5]])])

# fit measures of the observation-model result array
# save data for each variable separately (1) for all the lake (2) per monitoring station (3) per depth (4) per station and depth

nvars0 <- dim(ar1)[2]

fitm_total <- apply(ar1,2,fitm)[,2:nvars0]
fitm_station <- apply(ar1,c(2,3),fitm)[,2:nvars0,]
fitm_depth <- apply(ar1,c(2,4),fitm)[,2:nvars0,]
fitm_station_depth <- apply(ar1,c(2,3,4),fitm)[,2:nvars0,,]
layer_depth_station <- round(apply(ar1['model','zct',,,],c(1,2),mean,na.rm=TRUE),1)


# save results as a list of arrays
compare_model2obs <- list(metadata= paste(project, "model-to-observation comparison, between", stdt, "and",endt),
                          variables=vars,
                          observation_model_result_array=ar1,
                          mean_layer_depth_station = layer_depth_station,
                          fit_measures_total = fitm_total,
                          fit_measures_station = fitm_station,
                          fit_measures_depth = fitm_depth,
                          fit_measures_station_depth = fitm_station_depth)

saveRDS(compare_model2obs, paste0(path_Out,"compare_model2obs.RDS"))




