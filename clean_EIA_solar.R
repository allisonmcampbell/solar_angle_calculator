setwd("/Users/camp426/OneDrive - PNNL/Documents/PNNL/Projects/FY21/HydroWiresC/data/2020")
library(lubridate)
require(ggplot2)
library(tidyquant)
library(dplyr)
library(forecast)
df<-read.csv("solar_full.csv",header=T)
df$date<-as.POSIXct(df$date, format="%Y-%m-%d %H:%M:%S", tz="UTC")

#################
# functions
#################
ssin <- function(x){
  return(sin(pi/180*x))
}
ccos <- function(x){
  return(cos(pi/180*x))
}

# a = altitude
# t = tilt of PV panel
# g = gamma -- orientation east/west of PV panel
# d = date -- need to pull out declination
# w = omega -- hour/time in degrees from solar noon
cos_theta_i <- function(a,t,g,x) {
  aasin <- round(ccos(d(x))*ssin(w(x))/ccos(a),digits=7)
  gs <- 180/pi*aasin
  out <- ssin(a)*ccos(t) + ccos(a)*ssin(t)*ccos(g-gs)
  return(out)
}
d <- function(x){
  doy <- as.numeric(strftime(x,format="%j"))
  out <- 23.45*ssin(360*(284+doy)/365)
  return(out)
}
# omega == hour_angle
w <- function(x){
  hod <- hour(x)
  h<-ifelse(hod==0,23,hod-1)
  out<- h*15 - 180
  return(out)
}

# d = declination
# l = latitude
# w = omega
solar_altitude <- function(l,x){
  sin_a <- ssin(d(x))*ssin(l) + ccos(d(x))*ccos(w(x))*ccos(l)
  return(180/pi*asin(sin_a))
}
ideal_gen <- function(x,i){
  capacity <- ifelse(solar_altitude(lat[i],x)<0,0,cap)
  cos_incidence <- ifelse(cos_theta_i(solar_altitude(lat[i],x),tilt[i],gamma[i],x)<0,
                          0,cos_theta_i(solar_altitude(lat[i],x),tilt[i],gamma[i],x))
  ideal <- capacity*cos_incidence
  return(ideal)
}

##############
# loop for all BAs
##############
set.seed(1900)
N<- 100

BAnames <- colnames(df)[-1]
j = 1
dat <-vector("list",length(BAnames))
cleandf <-vector("list",length(BAnames))
zero_idx <-vector("list",length(BAnames))
toohigh_idx <-vector("list",length(BAnames))
toolow_idx <-vector("list",length(BAnames))
cleandf_after <-vector("list",length(BAnames))
for (BA in BAnames){
  cap <- df[1,BA]
  sprintf("BA = %s",BA)
  sprintf("capacity = %f",cap)
  lat_min <- df[2,BA]
  lat_max <- df[3,BA]
  solar <- df[-1:-3,c("date",BA)]
  lat<-runif(N,lat_min,lat_max)
  gamma<-runif(N,-90,90)
  tilt<-runif(N,0,lat)
  
  ideal_dfs <-vector("list",N)
  for (i in 1:N) {
    ideal_dfs[[i]] <- data.frame(datetime = solar$date,ideal = ideal_gen(solar$date,i))
  }

  dat[[j]] <- do.call(rbind,ideal_dfs)
  cleandf[[j]] <- data.frame(datetime = solar$date, PV = solar[,BA])
  dat2<- do.call(cbind,ideal_dfs)
  dat2 <- dat2[!duplicated(as.list(dat2))]
  
  dat2$PV <- solar[,BA]
  dat2$max <- do.call(pmax,dat2[2:101])
  dat2$fixedPV <- dat2$PV
  
  # identify data points that are positive at night
  zero_idx[[j]] <- which((dat2$PV > dat2$max) & (dat2$max == 0.0))
  # replace positive PV generation at night with zero
  dat2$fixedPV[zero_idx[[j]]] <- dat2$max[zero_idx[[j]]]
  # identify data points that are higher than "reasonable" 
  # meaning higher than any simulated PV values
  toohigh_idx[[j]] <- which((dat2$PV > dat2$max) & (dat2$max > 0.0))
  dat2$fixedPV[toohigh_idx[[j]]] <- NA
  #which(is.na(dat2$fixedPV))
  # replace identified values with interpolated values
  tt <- zoo(dat2$fixedPV,dat2$datetime)
  tt<- as.ts(tt)
  tt<- na.interp(tt)
  dat2$fixedPV[toohigh_idx[[j]]] <- tt[toohigh_idx[[j]]]
  
  # identify data points which are negative 
  # -- more than 5% of installed capacity
  toolow_idx[[j]] <- which((dat2$PV < -0.05*cap) )
  dat2$fixedPV[toolow_idx[[j]]] <- 0.
  
  cleandf_after[[j]] <- data.frame(datetime = solar$date, PV = dat2$fixedPV)
  #cleandf_after[-c(zero_idx,toohigh_idx,toolow_idx),]
  
  j = j + 1
}

j = 8
ggplot()+
  geom_point(data=dat[[j]],aes(x=datetime,y=ideal,color="Ideal"),alpha=0.25) +
  geom_point(data=cleandf[[j]][-c(zero_idx[[j]],toohigh_idx[[j]],toolow_idx[[j]]),],aes(x=datetime,y=PV,color="Actual"),alpha=0.7) +
  geom_point(data=cleandf_after[[j]][c(zero_idx[[j]],toohigh_idx[[j]],toolow_idx[[j]]),],aes(x=datetime,y=PV,color="Fixed"),alpha=0.7) +
  geom_point(data=cleandf[[j]][c(zero_idx[[j]],toohigh_idx[[j]],toolow_idx[[j]]),],aes(x=datetime,y=PV,color="Bad"),alpha=0.7) +
  labs(x="Date",y="PV Generation") +

  ggtitle(BAnames[[j]]) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual("",breaks=c("Ideal","Actual","Fixed","Bad"),values = c("black", "red","blue","green")) 

+
  coord_x_datetime(xlim=c("2020-12-24 04:00:00","2020-12-25 01:00:00"))

+
  coord_x_datetime(xlim=c("2020-11-06 01:00:00","2020-11-10 01:00:00"))
#ylim(-10,230) +
# CFE, j=7, is bad --> do not use

# manual fix for IID, j = 10
cleandf_after[[10]][2929:(2929+23),"PV"] <- cleandf_after[[10]][(2929-24*2):(2928-24),"PV"]
toohigh_idx[[10]]<-append(toohigh_idx[[10]],c(2929:(2929+23)))
cleandf_after[[10]][5544:(5544+23),"PV"] <- cleandf_after[[10]][5496:(5496+23),"PV"]
toohigh_idx[[10]]<-append(toohigh_idx[[10]],c(5544:(5544+23)))
cleandf_after[[10]][7488:(7488+23),"PV"] <- cleandf_after[[10]][7464:(7464+23),"PV"]
toohigh_idx[[10]]<-append(toohigh_idx[[10]],c(7488:(7488+23)))


# compare CAISO solar with actual
caiso<-read.csv("CAISO_solar.csv",header=T)
caiso$date<-as.POSIXct(caiso$date, format="%Y-%m-%d %H:%M:%S", tz="UTC")
j = 8
ggplot()+
  geom_point(data=dat[[j]],aes(x=datetime,y=ideal,color="Ideal"),alpha=0.25) +
  geom_point(data=cleandf[[j]][-c(zero_idx[[j]],toohigh_idx[[j]],toolow_idx[[j]]),],aes(x=datetime,y=PV,color="Actual"),alpha=0.7) +
  geom_point(data=cleandf_after[[j]][c(zero_idx[[j]],toohigh_idx[[j]],toolow_idx[[j]]),],aes(x=datetime,y=PV,color="Fixed"),alpha=0.7) +
  geom_point(data=cleandf[[j]][c(zero_idx[[j]],toohigh_idx[[j]],toolow_idx[[j]]),],aes(x=datetime,y=PV,color="Bad"),alpha=0.7) +
  geom_point(data=caiso,aes(x=date,y=Solar,color="CAISO"),alpha=0.7) +
  labs(x="Date",y="PV Generation") +
  ggtitle(BAnames[[j]]) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual("",breaks=c("Ideal","Actual","Fixed","Bad","CAISO"),values = c("black", "red","blue","green","purple")) +
  coord_x_datetime(xlim=c("2020-11-30 01:00:00","2020-12-03 01:00:00"))

#########
# write to file
#########
dat3<- do.call(cbind,cleandf_after)
dat3 <- dat3[!duplicated(as.list(dat3))]

library(data.table)
setnames(dat3,old=colnames(dat3)[-1],new=BAnames)
write.csv(dat3,"solar_clean.csv",row.names=F)
