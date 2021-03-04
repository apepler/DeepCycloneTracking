# THis is additional code from responding to reviews which assesses statistical significance
# Of latitudinal mean changes

# It relies on annual versions of the rain grids - these have not been uploaded to figshare 
# due to the larger data size, but can be provided on request - acacia.pepler@bom.gov.au

# Initial set up of libraries and functions

rm(list=ls())
library(raster)
library(sp)
library(ncdf4)
library(abind)
library(viridisLite)
library(RColorBrewer)
library(fields)
library(maps)
library(oz)

# A function that makes nice colour bars
ColorBar <- function(brks,cols,vert=T,subsampleg=1,skip1=T)
{
  if(vert) {
    par(mar = c(2, 1, 2, 4), mgp = c(1, 1, 0), las = 1, cex = 1)
    image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    if(skip1)
      axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE,
           labels = brks[seq(2, length(brks)-1, subsampleg)]) else
             axis(4, at = seq(0.5, length(brks) - 1.5, subsampleg), tick = TRUE,
                  labels = brks[seq(1, length(brks)-1, subsampleg)])
  } else {
    par(mar = c(1.5, 1, 1, 1), mgp = c(1.5, 0.3, 0), las = 1, cex = 1)
    image(1:length(cols), 1, t(t(1:length(cols))), axes = FALSE, col = cols,
          xlab = '', ylab = '')
    box()
    if(skip1)
      axis(1, at = seq(1.5, length(brks) - 1.5, subsampleg),
           labels = brks[seq(2, length(brks)-1, subsampleg)]) else
             axis(1, at = seq(0.5, length(brks) - 1.5, subsampleg),
                  labels = brks[seq(1, length(brks)-1, subsampleg)])
  }
}

# Set up data directories, thresholds etc

basedir="/g/data/eg3/asp561/CycloneTracking/CMIP5/final_netCDF/"

cmiplist=c("ACCESS1-0","ACCESS1-3","CCSM4","GFDL-ESM2M","MRI-CGCM3","CNRM-CM5","MIROC5","NorESM1-M","CanESM2","HadGEM2-CC")
memlist=c("r1i1p1","r1i1p1","r6i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1")

colseq=c("orange","magenta3","red4","blue4")
catnames=c("Other","Deep cyclone","Shallow surface","Shallow upper")

cmip=rep(T,10)
thresh="rad2cv0.5"
threshD="rad2cv1"
proj=100
lev=500
years1=1979:2005
years2=2070:2099
dist=500
pthresh=99
pthresh2=99.7

rcp=c("historical","rcp85","rcp45")
i=1
type1="cv1var_dist500event"

breaks1=c(-100,seq(-60,60,10),1000)
cols1=colorRampPalette(brewer.pal(11,"RdBu"))(length(breaks1)-1)
dnames=c("No","Deep","Shallow surface","Shallow upper")

snames=c("MJJASO","NDJFMA","Annual")
mlist=list(5:10,c(11:12,1:4),1:12)
s=3

#### Using the reggridded versions

lat=seq(-89.5,89.5)

histlat_freq<-array(NaN,c(length(cmiplist),1,length(years1),4,length(lat),12))
dimnames(histlat_freq)[[1]]=cmiplist
dimnames(histlat_freq)[[2]]=rcp
dimnames(histlat_freq)[[4]]=c("No","Deep","Shallow surface","Shallow upper")
histlat_raindays<-histlat_raindayrain<-histlat_rain<-histlat_threshdays<-histlat_freq

projlat_freq<-array(NaN,c(length(cmiplist),2,length(years2),4,length(lat),12))
dimnames(projlat_freq)[[1]]=cmiplist
dimnames(projlat_freq)[[2]]=rcp
dimnames(projlat_freq)[[4]]=c("No","Deep","Shallow surface","Shallow upper")
projlat_raindayrain<-projlat_raindays<-projlat_rain<-projlat_threshdays<-projlat_freq


for(i in 1:length(cmiplist))
{
  # Load all historical data for one dataset
  r=1
  print(cmiplist[i])
  raindays1<-raindayrain1<-threshdays1<-freq1<-rain1<-array(NaN,c(length(years1),4,length(lon),length(lat),12))

  sdir=paste0(basedir,"CMIP5/",cmiplist[i],"/",rcp[r],"/",memlist[i],"/proj",proj,"_lows_",thresh,"/")

  for(y in 1:length(years1))
  {
    fname=paste0(cmiplist[i],"_",rcp[r],"_",memlist[i],"_deepcyclones_dailygrid_",type1,"_",
                 years1[y],"_rain_Q",pthresh,"_10deg_monthly_regrid.nc")
    a=nc_open(paste0(sdir,fname))
    rain1[y,,,,]=ncvar_get(a,"cyc_rain")
    freq1[y,,,,]=ncvar_get(a,"cyc_freq")
    nc_close(a)

    fname=paste0(cmiplist[i],"_",rcp[r],"_",memlist[i],"_deepcyclones_dailygrid_",type1,"_",
                 years1[y],"_rain_Q",pthresh2,"_10deg_monthly_regrid.nc")

    a=nc_open(paste0(sdir,fname))
    threshdays1[y,,,,]=ncvar_get(a,"cyc_threshdays")
    nc_close(a)

    fname=paste0(cmiplist[i],"_",rcp[r],"_",memlist[i],"_deepcyclones_dailygrid_",type1,"_",
                 years1[y],"_rain_Above1mm_10deg_monthly_regrid.nc")

    a=nc_open(paste0(sdir,fname))
    raindayrain1[y,,,,]=ncvar_get(a,"cyc_rain")
    raindays1[y,,,,]=ncvar_get(a,"cyc_threshdays")
    nc_close(a)
  }

  ## Make the necessary bits of the datasets
  histlat_freq[i,1,,,,]=apply(freq1,c(1,2,4,5),sum)
  histlat_rain[i,1,,,,]=apply(rain1,c(1,2,4,5),sum)
  histlat_threshdays[i,1,,,,]=apply(threshdays1,c(1,2,4,5),sum)
  histlat_raindays[i,1,,,,]=apply(raindays1,c(1,2,4,5),sum)
  histlat_raindayrain[i,1,,,,]=apply(raindayrain1,c(1,2,4,5),sum)

  for(r in 2:3)
  {
    raindays2<-raindayrain2<-threshdays2<-freq2<-rain2<-array(NaN,c(length(years2),4,length(lon),length(lat),12))

    sdir=paste0(basedir,"CMIP5/",cmiplist[i],"/",rcp[r],"/",memlist[i],"/proj",proj,"_lows_",thresh,"/")

    for(y in 1:length(years2))
    {
      fname=paste0(cmiplist[i],"_",rcp[r],"_",memlist[i],"_deepcyclones_dailygrid_",type1,"_",
                   years2[y],"_rain_Q",pthresh,"_10deg_monthly_regrid.nc")
      a=nc_open(paste0(sdir,fname))
      rain2[y,,,,]=ncvar_get(a,"cyc_rain")
      freq2[y,,,,]=ncvar_get(a,"cyc_freq")
      nc_close(a)

      fname=paste0(cmiplist[i],"_",rcp[r],"_",memlist[i],"_deepcyclones_dailygrid_",type1,"_",
                   years2[y],"_rain_Q",pthresh2,"_10deg_monthly_regrid.nc")

      a=nc_open(paste0(sdir,fname))
      threshdays2[y,,,,]=ncvar_get(a,"cyc_threshdays")
      nc_close(a)

      fname=paste0(cmiplist[i],"_",rcp[r],"_",memlist[i],"_deepcyclones_dailygrid_",type1,"_",
                   years2[y],"_rain_Above1mm_10deg_monthly_regrid.nc")

      a=nc_open(paste0(sdir,fname))
      raindayrain2[y,,,,]=ncvar_get(a,"cyc_rain")
      raindays2[y,,,,]=ncvar_get(a,"cyc_threshdays")
      nc_close(a)
    }

    ## Make the necessary bits of the datasets
    projlat_freq[i,r-1,,,,]=apply(freq2,c(1,2,4,5),sum)
    projlat_rain[i,r-1,,,,]=apply(rain2,c(1,2,4,5),sum)
    projlat_threshdays[i,r-1,,,,]=apply(threshdays2,c(1,2,4,5),sum)
    projlat_raindays[i,r-1,,,,]=apply(raindays2,c(1,2,4,5),sum)
    projlat_raindayrain[i,r-1,,,,]=apply(raindayrain2,c(1,2,4,5),sum)

  }
}



### Look at regional mean change and statistical significance

latreg=c("SH_polar","SH_extratropics","SH_midlat","Tropics",
         "NH_midlat","NH_extratropics","NH_polar","All_extratropics","All_midlat","Global")
latlist=list(which(lat<=-75),which(lat>-75 & lat<=-45),which(lat>-45 & lat<=-15),
             which(abs(lat)<15),which(lat>=15 & lat<45),which(lat>=45 & lat<75),which(lat>=75),
             which(abs(lat)>=45 & abs(lat)<75),which(abs(lat)>=15 & abs(lat)<45),which(abs(lat)<=90))

cycchange=array(0,c(length(latreg),10,5,5,7,3))
dimnames(cycchange)[[1]]=latreg
dimnames(cycchange)[[2]]=cmiplist
dimnames(cycchange)[[3]]=c("No cyclone","Deep cyclone","Shallow surface","Shallow upper","All")
dimnames(cycchange)[[4]]=c("Days","Total rainfall","Rain days","Rain rate","P99.7 days")
dimnames(cycchange)[[5]]=c("1979-2005","2070-2099 rcp85","2070-2099 rcp45","change rcp85","change rcp45","sig rcp85","sig rcp45")
dimnames(cycchange)[[6]]=snames

for(i in 1:10)
  for(j in 1:5)
    for(s in 1:3)
    {
      tmp1=switch(j,apply(histlat_freq[i,1,,,,mlist[[s]]],c(1,2,3),sum),
                  apply(histlat_rain[i,1,,,,mlist[[s]]],c(1,2,3),sum),
                  apply(histlat_raindays[i,1,,,,mlist[[s]]],c(1,2,3),sum),
                  apply(histlat_raindayrain[i,1,,,,mlist[[s]]],c(1,2,3),sum)/apply(histlat_raindays[i,1,,,,mlist[[s]]],c(1,2,3),sum),
                  apply(histlat_threshdays[i,1,,,,mlist[[s]]],c(1,2,3),sum))
      
      tmp2=switch(j,apply(projlat_freq[i,,,,,mlist[[s]]],c(1,2,3,4),sum),
                  apply(projlat_rain[i,,,,,mlist[[s]]],c(1,2,3,4),sum),
                  apply(projlat_raindays[i,,,,,mlist[[s]]],c(1,2,3,4),sum),
                  apply(projlat_raindayrain[i,,,,,mlist[[s]]],c(1,2,3,4),sum)/apply(projlat_raindays[i,,,,,mlist[[s]]],c(1,2,3,4),sum),
                  apply(projlat_threshdays[i,,,,,mlist[[s]]],c(1,2,3,4),sum))
      
      for(l in 1:length(latreg))
      {
        for(k in 1:4)
        {
          cycchange[l,i,k,j,1,s]=mean(apply(tmp1[,k,latlist[[l]]],1,sum,na.rm=T))
          cycchange[l,i,k,j,2:3,s]=apply(apply(tmp2[,,k,latlist[[l]]],c(1,2),sum,na.rm=T),1,mean)

          for(r in 1:2)
            {
            a=t.test(apply(tmp2[r,,k,latlist[[l]]],1,sum,na.rm=T),apply(tmp1[,k,latlist[[l]]],1,sum,na.rm=T))
            cycchange[l,i,k,j,5+r,s]=a$p.value
          }
        }
        cycchange[l,i,5,j,1,s]=mean(apply(tmp1[,,latlist[[l]]],1,sum,na.rm=T))
        cycchange[l,i,5,j,2:3,s]=apply(apply(tmp2[,,,latlist[[l]]],c(1,2),sum,na.rm=T),1,mean)
        
        for(r in 1:2)
        {
          a=t.test(apply(tmp2[r,,,latlist[[l]]],1,sum,na.rm=T),apply(tmp1[,,latlist[[l]]],1,sum,na.rm=T))
          cycchange[l,i,5,j,5+r,s]=a$p.value
        }
      }
    }

cycchange[,,,,4,]=100*((cycchange[,,,,2,]/cycchange[,,,,1,])-1)
cycchange[,,,,5,]=100*((cycchange[,,,,3,]/cycchange[,,,,1,])-1)
cycchange[,"MIROC5",,,c(3,5,7),]=NaN

