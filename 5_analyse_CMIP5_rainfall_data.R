# THis has a set of code that I run interactively to analyse cyclone rainfall

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

#Model details
cmiplist=c("ERA5","ACCESS1-0","ACCESS1-3","CCSM4","GFDL-ESM2M","MRI-CGCM3","CNRM-CM5","MIROC5","NorESM1-M","CanESM2","HadGEM2-CC")
memlist=c(NaN,"r1i1p1","r1i1p1","r6i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1")
cmip=c(F,rep(T,10))

# Details to call the correct files
type1="cv1var_dist500event" # Intensity threshold
pthresh=99 
pthresh2=99.7
rcp=c("historical","rcp85","rcp45")
ylist=list(1979:2005,2070:2099,2070:2099) # Years in each file

# Some plotting details
colseq=c("orange","magenta3","red4","blue4")
catnames=c("Other","Deep cyclone","Shallow surface","Shallow upper")
years=1950:2100

snames=c("MAM","JJA","SON","DJF","MJJASO","NDJFMA","Annual","AprSept","OctMar")
mlist=list(3:5,6:8,9:11,c(12,1:2),5:10,c(11:12,1:4),1:12,4:9,c(10:12,1:3))
s=7

#### Set up the arrays

lat=seq(-89.5,89.5)
lon=seq(0.5,359.5)

cyc_freq<-array(NaN,c(length(cmiplist),3,4,length(lon),length(lat),12))
dimnames(cyc_freq)[[1]]=cmiplist
dimnames(cyc_freq)[[2]]=rcp
dimnames(cyc_freq)[[3]]=c("No","Deep","Shallow surface","Shallow upper")
cyc_threshdays2<-cyc_raindayrain<-cyc_raindays<-cyc_rain<-cyc_threshdays<-cyc_freq

# Load all the data - available at YYYY

for(i in 1:length(cmiplist))
{
  print(cmiplist[i])
    for(r in 1:3) 
    {
      if(i==1 & r>1) break
      if(r==3 & i==8) break
      
      years=ylist[[r]]

      if(i==1) fname=paste0(cmiplist[i],"_deepcyclones_dailygrid_",type1,"_rain_Q",pthresh,"_",min(years),max(years),"_10deg_v2_regrid_monthly.nc") else
        fname=paste0(cmiplist[i],"_",rcp[r],"_",memlist[i],"_deepcyclones_dailygrid_",type1,"_",min(years),max(years),"_rain_Q",pthresh,"_10deg_monthly_regrid.nc")
      
      a=nc_open(paste0(basedir,fname))
      cyc_rain[i,r,,,,]=ncvar_get(a,"cyc_rain")/length(years)
      cyc_freq[i,r,,,,]=ncvar_get(a,"cyc_freq")/length(years)
      cyc_threshdays[i,r,,,,]=ncvar_get(a,"cyc_threshdays")/length(years)
      nc_close(a)
      
      if(i==1) fname=paste0(cmiplist[i],"_deepcyclones_dailygrid_",type1,"_rain_Q",pthresh2,"_",min(years),max(years),"_10deg_v2_regrid_monthly.nc") else
        fname=paste0(cmiplist[i],"_",rcp[r],"_",memlist[i],"_deepcyclones_dailygrid_",type1,"_",min(years),max(years),"_rain_Q",pthresh2,"_10deg_monthly_regrid.nc")
      
      a=nc_open(paste0(basedir,fname))
      cyc_threshdays2[i,r,,,,]=ncvar_get(a,"cyc_threshdays")/length(years)
      nc_close(a)
      
      if(i==1) fname=paste0(cmiplist[i],"_deepcyclones_dailygrid_",type1,"_rain_Above1mm_",min(years),max(years),"_10deg_v2_regrid_monthly.nc") else
        fname=paste0(cmiplist[i],"_",rcp[r],"_",memlist[i],"_deepcyclones_dailygrid_",type1,"_",min(years),max(years),"_rain_Above1mm_10deg_monthly_regrid.nc")
      
      a=nc_open(paste0(basedir,fname))
      cyc_raindayrain[i,r,,,,]=ncvar_get(a,"cyc_rain")/length(years)
      cyc_raindays[i,r,,,,]=ncvar_get(a,"cyc_threshdays")/length(years)
      nc_close(a)
    }
  }


## Gridded figures

#Supplementary Figure S6 - Mean frequency and bias

breaks1=c(-1000,seq(-20,20,5),1000)
cols1=colorRampPalette(brewer.pal(11,"RdBu"))(length(breaks1)-1)
breaks2=c(0,2.5,seq(5,20,5),seq(30,60,10),1000)
cols2=rev(viridis(length(breaks2)-1))
dnames=c("No","Deep","Shallow surface","Shallow upper")
rnames=c("ERA5","CMIP5 1979-2005","CMIP5 2070-2099")

pnum=1
for(s in 7)
{
  pdf(file=paste0("FigureS6_CMIP5_global_cycfreq_",type1,"_19792005_cyclonefreqprop_global_",snames[s],"_ensmean_bias_10models.pdf"),width=16,height=9,pointsize=16)
  par(mar=c(2,2,3,1))
  layout(rbind(c(1,2,10,3,11),c(4,5,10,6,11),c(7,8,10,9,11)),width=c(1,1,0.35,1,0.35))
  
  for(t in 2:4)
  {
    # Count of a given type divided by count of all days
    tmp2=apply(cyc_freq[,1,t,,,mlist[[s]]],c(1,2,3),sum)/apply(cyc_freq[,1,,,,mlist[[s]]],c(1,3,4),sum)
    
    image(lon,lat,100*tmp2[1,,],
          breaks=breaks2,col=cols2,main=paste0(letters[pnum],") ",dnames[t]," cyclones: ERA5"),
          xlab="",ylab="",cex.main=1.5)
    map("world2",add=T)
    pnum=pnum+1
    image(lon,lat,100*apply(tmp2[2:length(cmiplist),,],c(2,3),mean),
          breaks=breaks2,col=cols2,main=paste0(letters[pnum],") ",dnames[t]," cyclones: CMIP5"),
          xlab="",ylab="",cex.main=1.5)
    map("world2",add=T)
    pnum=pnum+1
    for(i in 2:length(cmiplist)) tmp2[i,,]=tmp2[i,,]-tmp2[1,,]
    image(lon,lat,100*apply(tmp2[2:length(cmiplist),,],c(2,3),mean),
          breaks=breaks1,col=cols1,main=paste0(letters[pnum],") ",dnames[t]," cyclone bias"),
          xlab="",ylab="",cex.main=1.5)
    map("world2",add=T)
    pnum=pnum+1
  }
  ColorBar(brks=paste0(breaks2,"%"),cols=cols2,vert=T,subsampleg=1)
  ColorBar(brks=paste0(breaks1,"%"),cols=cols1,vert=T,subsampleg=1)
  dev.off()
  
}

# Figure 2 - Spatial plot of mean change in frequency for rcp8.5 (and also FIgure S1 with minor edits)

breaks1=c(-100,seq(-60,60,10),100000)
cols1=colorRampPalette(brewer.pal(11,"RdBu"))(length(breaks1)-1)
pnames=c("a) Change in deep cyclones",
         "b) Change in shallow surface cyclones",
         "c) Change in shallow upper cyclones")

for(s in 7)
{

pdf(file=paste0("Figure2_CMIP5_global_cycfreq_",type1,"_19792005_cyclonechangePC_20702099_global_",snames[s],"_ensmean_points_10models.pdf"),width=16,height=3,pointsize=14)
par(mar=c(2,2,3,1))
layout(cbind(1,2,3,4),width=c(1,1,1,0.3))

for(t in 2:4)
  {
    tmp1=apply(cyc_freq[2:length(cmiplist),,t,,,mlist[[s]]],c(1,2,3,4),sum)
    
    tmp1a=tmp1[,2,,]/tmp1[,1,,] # Use tmp1[,3,,]/tmp1[,1,,] for rcp4.5
    tmp1a[tmp1[,1,,]<1]=NaN
    
    image(lon,lat,100*(apply(tmp1a,c(2,3),mean,na.rm=T)-1),
          breaks=breaks1,col=cols1,main=pnames[t-1],
          xlab="",ylab="",cex.main=1.5)
    map("world2",add=T,resolution=0)
    
    postrend<-negtrend<-array(0,c(length(lon),length(lat)))
    for(j in 1:10) # Use c(1:6,8:10)) for rcp4.5 as MIROC5 is missing
    {
      I=which(!is.na(tmp1a[j,,]) & tmp1a[j,,]>1)
      postrend[I]=postrend[I]+1
      I=which(!is.na(tmp1a[j,,]) & tmp1a[j,,]<1)
      negtrend[I]=negtrend[I]+1
    }
    I=which(!is.na(apply(tmp1a,c(2,3),mean,na.rm=T)) & postrend<8 & negtrend<8,arr.ind=T)
    points(lon[I[,1]],lat[I[,2]],cex=0.1,pch=16)
  }
ColorBar(brks=paste0(breaks1,"%"),cols=cols1,vert=T,subsampleg=1)
dev.off()
}

# Figure S2 - Two seasons version

  breaks1=c(-100,seq(-60,60,10),100000)
  cols1=colorRampPalette(brewer.pal(11,"RdBu"))(length(breaks1)-1)
  pnames=c("change in deep cyclones",
           "change in shallow surface cyclones",
           "change in shallow upper cyclones")
  
  pnum=1
  pdf(file=paste0("FigureS2_CMIP5_global_cycfreq_",type1,"_19792005_cyclonechangePC_20702099_global_vseason_ensmean_points_10models.pdf"),width=16,height=6,pointsize=14)
  par(mar=c(2,2,3,1))
  layout(rbind(c(1:3,7),4:7),width=c(1,1,1,0.3))
  
  for(s in 5:6)
  for(t in 2:4)
  {
    tmp1=apply(cyc_freq[2:length(cmiplist),,t,,,mlist[[s]]],c(1,2,3,4),sum)
    
    tmp1a=tmp1[,2,,]/tmp1[,1,,] # Use tmp1[,3,,]/tmp1[,1,,] for rcp4.5
    tmp1a[tmp1[,1,,]<1]=NaN
    
    image(lon,lat,100*(apply(tmp1a,c(2,3),mean,na.rm=T)-1),
          breaks=breaks1,col=cols1,main=paste0(letters[pnum],") ",snames[s]," ",pnames[t-1]),
          xlab="",ylab="",cex.main=1.5)
    map("world2",add=T,resolution=0)
    pnum=pnum+1
    
    postrend<-negtrend<-array(0,c(length(lon),length(lat)))
    for(j in 1:10) # Use c(1:6,8:10)) for rcp4.5 as MIROC5 is missing
    {
      I=which(!is.na(tmp1a[j,,]) & tmp1a[j,,]>1)
      postrend[I]=postrend[I]+1
      I=which(!is.na(tmp1a[j,,]) & tmp1a[j,,]<1)
      negtrend[I]=negtrend[I]+1
    }
    I=which(!is.na(apply(tmp1a,c(2,3),mean,na.rm=T)) & postrend<8 & negtrend<8,arr.ind=T)
    points(lon[I[,1]],lat[I[,2]],cex=0.1,pch=16)
  }
  ColorBar(brks=paste0(breaks1,"%"),cols=cols1,vert=T,subsampleg=1)
  dev.off()

# Figure 4 - Plot of absolute change in days with extreme rain

for(s in 7)
{
  breaks1=c(-100,seq(-1,1,0.2),100000)
  cols1=colorRampPalette(brewer.pal(11,"RdBu"))(length(breaks1)-1)
  pnames=c("d) Change in other days with extreme rain",
           "a) Change in deep cyclones with extreme rain",
           "b) Change in shallow surface cyclones with extreme rain",
           "c) Change in shallow upper cyclones with extreme rain")
  
  pdf(file=paste0("Figure4_CMIP5_global_cycrain_P99.7_",type1,"_19792005_cyclonechangeabs_20702099_global_",snames[s],"_ensmean_points_10models.pdf"),width=11,height=6,pointsize=14)
  par(mar=c(2,2,3,1))
  layout(rbind(c(1,2,5),3:5),width=c(1,1,0.3))
  
  for(t in c(2:4,1))
  {
    # Average annual number of days above percentile threshold
    tmp1=apply(cyc_threshdays2[2:length(cmiplist),,t,,,mlist[[s]]],c(1,2,3,4),sum)
    tmp1a=tmp1[,2,,]-tmp1[,1,,]
    
    image(lon,lat,apply(tmp1a,c(2,3),mean,na.rm=T),
          breaks=breaks1,col=cols1,main=pnames[t],
          xlab="",ylab="",cex.main=1.2)
    map("world2",add=T,resolution=0)
    
    postrend<-negtrend<-array(0,c(length(lon),length(lat)))
    for(j in c(1:6,8:10))
    {
      I=which(!is.na(tmp1a[j,,]) & tmp1a[j,,]>0)
      postrend[I]=postrend[I]+1
      I=which(!is.na(tmp1a[j,,]) & tmp1a[j,,]<0)
      negtrend[I]=negtrend[I]+1
    }
    I=which(!is.na(apply(tmp1a,c(2,3),mean,na.rm=T)) & postrend<8 & negtrend<8,arr.ind=T)
    points(lon[I[,1]],lat[I[,2]],cex=0.1,pch=16)
  }
  ColorBar(brks=breaks1,cols=cols1,vert=T,subsampleg=1)
  dev.off()
}

# Supplementary Figure S4 - Change in mean rain rate
  
for(s in 7)
{
  breaks1=c(-100,seq(-50,50,10),100000)
  cols1=colorRampPalette(brewer.pal(11,"RdBu"))(length(breaks1)-1)
  pnames=c("d) Change in rain rate on other days",
           "a) Change in rain rate for deep cyclones",
           "b) Change in rain rate for shallow surface cyclones",
           "c) Change in rain rate for shallow upper cyclones")
  
  pdf(file=paste0("FigureS4_CMIP5_global_cycrainrate_",type1,"_19792005_cyclonechangePC_20702099_global_",snames[s],"_ensmean_points_10models.pdf"),width=11,height=6,pointsize=14)
  par(mar=c(2,2,3,1))
  layout(rbind(c(1,2,5),3:5),width=c(1,1,0.3))
  
  for(t in c(2:4,1))
  {
    # Average rain rate, defined as total rain divided by number of days in category
    tmp1=apply(cyc_rain[2:length(cmiplist),,t,,,mlist[[s]]],1:4,sum)/apply(cyc_freq[2:length(cmiplist),,t,,,mlist[[s]]],1:4,sum)
    tmp1a=tmp1[,2,,]/tmp1[,1,,]
    
    freq=apply(cyc_freq[2:length(cmiplist),1,t,,,mlist[[s]]],1:3,sum)
    tmp1a[freq<1]=NaN
    
    image(lon,lat,100*(apply(tmp1a,c(2,3),mean)-1),
          breaks=breaks1,col=cols1,main=pnames[t],
          xlab="",ylab="",cex.main=1.2)
    map("world2",add=T,resolution=0)
    
    postrend<-negtrend<-array(0,c(length(lon),length(lat)))
    for(j in 1:10)
    {
      I=which(!is.na(tmp1a[j,,]) & tmp1a[j,,]>1)
      postrend[I]=postrend[I]+1
      I=which(!is.na(tmp1a[j,,]) & tmp1a[j,,]<1)
      negtrend[I]=negtrend[I]+1
    }
    I=which(!is.na(apply(tmp1a,c(2,3),mean)) & postrend<8 & negtrend<8,arr.ind=T)
    points(lon[I[,1]],lat[I[,2]],cex=0.1,pch=16)
  }
  ColorBar(brks=paste0(breaks1,"%"),cols=cols1,vert=T,subsampleg=2)
  dev.off()
}

  # How has the rain rate changed on average across the globe?
tmp1=apply(cyc_rain[2:11,,,,,mlist[[s]]],1:5,sum)/apply(cyc_freq[2:11,,,,,mlist[[s]]],1:5,sum)
tmp1a=tmp1[,2,,,]/tmp1[,1,,,]
freq=apply(cyc_freq[2:11,1,,,,mlist[[s]]],1:4,sum)
tmp1a[freq<1]=NaN # Only calculate where at least one cyclone day per year in the historical period
b=apply(tmp1a,c(1,2),mean,na.rm=T)
apply(b,2,mean)

I=which(abs(lat)>=15 & abs(lat)<=75) # What about excluding tropics & poles?
b=apply(tmp1a[,,,I],c(1,2),mean,na.rm=T)
apply(b,2,mean)

## Plots averaged across latitudes

## Figure 1 - Plot frequency/rainfall by latitude
s=7
pnames=c("a) Proportion of all days","b) Proportion of total rainfall","c) Proportion of extreme rain days")
pdf(file=paste0("Figure1_CMIP5_global_cycprop_P99.7_",type1,"_19792005_",snames[s],"_vlat_multipanel2_10models.pdf"),height=5,width=20,pointsize=18)
layout(cbind(1,2,3))
par(mar=c(2,2,4,1))

for(k in 1:3)
{
  tmp=switch(k,apply(apply(cyc_freq[,1,,,,mlist[[s]]],1:4,sum),c(1,2,4),mean),
             apply(apply(cyc_rain[,1,,,,mlist[[s]]],1:4,sum),c(1,2,4),mean),
             apply(apply(cyc_threshdays2[,1,,,,mlist[[s]]],1:4,sum),c(1,2,4),mean))
  
  tmp2=apply(tmp,c(1,3),sum)
  for(i in 1:4) tmp[,i,]=100*tmp[,i,]/tmp2
  
  plot(NA,xlim=c(-90,90),ylim=c(0,110),xlab="",ylab="",main=pnames[k],axes=F,cex.main=1.5)
  axis(1,at=seq(-80,80,20),labels=c(paste0(seq(80,20,-20),"S"),"0",paste0(seq(20,80,20),"N")))
  axis(2,at=seq(0,100,20),labels=paste0(seq(0,100,20),"%"))
  box()
  
  for(i in 1:4) {
    crgb=col2rgb(colseq[i])
    polygon(c(lat,rev(lat)),c(apply(tmp[2:length(cmiplist),i,],2,min),rev(apply(tmp[2:length(cmiplist),i,],2,max))),
            col=rgb(crgb[1],crgb[2],crgb[3],64,maxColorValue =255),border=NA)
    lines(lat,apply(tmp[2:length(cmiplist),i,],2,mean),lwd=3,col=colseq[i],lty=2)
    lines(lat,tmp[1,i,],lwd=3,col=colseq[i])
    
  }
  abline(h=0,lty=2)
  legend("topleft",catnames[c(2,3,4,1)],#c("No cyclone","Deep cyclone","Shallow surface","Shallow upper"),
         bty="n",lwd=2,col=colseq[c(2,3,4,1)],ncol=2)
  legend("topright",c("ERA5","CMIP5"),lwd=2,lty=1:2,col="black",bty="n")
}
dev.off()

## Figure 3 - Change vs latitude, including total rainfall/total number of days

s=7
pnames=c("a) Change in rain days","b) Change in total rainfall","c) Change in extreme rain days")
pdf(file=paste0("Figure3_CMIP5_global_cycfreq_",type1,"_19792005_changedays_Q99.7_absolute_20702099_",snames[s],"_vlat_multipanel_10models.pdf"),height=5,width=20,pointsize=18)
layout(cbind(1,2,3))
par(mar=c(2,2,4,1))
ylims=rbind(c(-30,70),c(-250,600),c(-1,6))
yaxs=list(seq(-40,80,20),seq(-400,600,200),seq(-1,6))

for(k in 1:3)
{
  tmp=switch(k,apply(cyc_raindays[2:length(cmiplist),,,,,mlist[[s]]],1:5,sum),
             apply(cyc_rain[2:length(cmiplist),,,,,mlist[[s]]],1:5,sum),
             apply(cyc_threshdays2[2:length(cmiplist),,,,,mlist[[s]]],1:5,sum))
  
  tmp1=apply(tmp[,2,,,],c(1,2,4),mean)-apply(tmp[,1,,,],c(1,2,4),mean)
  tmp2=apply(tmp1,c(2,3),mean)
  tmp3=apply(tmp1,c(1,3),sum) # Change in total heavy rain days
  
  plot(NA,xlim=c(-90,90),ylim=ylims[k,],xlab="",ylab="",main=pnames[k],axes=F,cex.main=1.5)
  axis(1,at=seq(-80,80,20),labels=c(paste0(seq(80,20,-20),"S"),"0",paste0(seq(20,80,20),"N")))
  axis(2,at=yaxs[[k]])
  box()
  
  box()
  polygon(c(lat,rev(lat)),c(apply(tmp3,2,min),rev(apply(tmp3,2,max))),
          col=rgb(0,0,0,1/4),border=NA)
  abline(h=0,lty=2)
  lines(lat,apply(tmp3,2,mean),col="black",lwd=3)
  cols=viridis(5)
  for(i in 1:4) lines(lat,tmp2[i,],lwd=3,col=colseq[i])
  legend("topleft",c(catnames[c(2,3,4,1)],"All days"),
         bty="n",lwd=3,col=c(colseq[c(2,3,4,1)],"black"),ncol=3)
}
dev.off()


### Data for analysis based on latitude bands
## 1. Mean frequencies and contributions to rainfall means for different latitude bands

latreg=c("SH_polar","SH_extratropics","SH_midlat","Tropics",
         "NH_midlat","NH_extratropics","NH_polar","All_extratropics","All_midlat","Global")
latlist=list(which(lat<=-75),which(lat>-75 & lat<=-45),which(lat>-45 & lat<=-15),
             which(abs(lat)<15),which(lat>=15 & lat<45),which(lat>=45 & lat<75),which(lat>=75),
             which(abs(lat)>=45 & abs(lat)<75),which(abs(lat)>=15 & abs(lat)<45),which(abs(lat)<=90))
s=6

cycstats=array(0,c(length(latreg),length(cmiplist),5,6,3))
dimnames(cycstats)[[1]]=latreg
dimnames(cycstats)[[2]]=cmiplist
dimnames(cycstats)[[3]]=c("All days","No cyclone","Deep cyclone","Shallow surface","Shallow upper")
dimnames(cycstats)[[4]]=c("Days","Total rainfall","Rain days","Rain rate","P99 days","P99.7 days")
dimnames(cycstats)[[5]]=snames[5:7]

for(i in 1:length(cmiplist))
  for(j in 1:6)
    for(s in 5:7)
    {
      tmp1=switch(j,apply(cyc_freq[i,1,,,,mlist[[s]]],c(1,2,3),sum),
                  apply(cyc_rain[i,1,,,,mlist[[s]]],c(1,2,3),sum),
                  apply(cyc_raindays[i,1,,,,mlist[[s]]],c(1,2,3),sum),
                  apply(cyc_raindayrain[i,1,,,,mlist[[s]]],c(1,2,3),sum)/apply(cyc_raindays[i,1,,,,mlist[[s]]],c(1,2,3),sum),
                  apply(cyc_threshdays[i,1,,,,mlist[[s]]],c(1,2,3),sum),
                  apply(cyc_threshdays2[i,1,,,,mlist[[s]]],c(1,2,3),sum))
      
      for(k in 1:5)
      {
        if(k==1) tmp2=apply(tmp1,c(2,3),sum) else # Total frequency
          tmp2=tmp1[k-1,,]/apply(tmp1,c(2,3),sum) # Proportion for each type
        
        for(l in 1:length(latreg)) cycstats[l,i,k,j,s-4]=mean(tmp2[,latlist[[l]]],na.rm=T)
        
      }
    }

## 2. Projected changes in frequency 

cycchange=array(0,c(length(latreg),10,5,6,5,3))
dimnames(cycchange)[[1]]=latreg
dimnames(cycchange)[[2]]=cmiplist[2:11]
dimnames(cycchange)[[3]]=c("No cyclone","Deep cyclone","Shallow surface","Shallow upper","All")
dimnames(cycchange)[[4]]=c("Days","Total rainfall","Rain days","Rain rate","P99 days","P99.7 days")
dimnames(cycchange)[[5]]=c("1979-2005","2070-2099 rcp85","2070-2099 rcp45","change rcp85","change rcp45")
dimnames(cycchange)[[6]]=snames[5:7]

for(i in 1:10)
  for(j in 1:6)
    for(s in 5:7)
    {
      tmp1=switch(j,apply(cyc_freq[i+1,,,,,mlist[[s]]],c(1,2,3,4),sum),
                  apply(cyc_rain[i+1,,,,,mlist[[s]]],c(1,2,3,4),sum),
                  apply(cyc_raindays[i+1,,,,,mlist[[s]]],c(1,2,3,4),sum),
                  apply(cyc_raindayrain[i+1,,,,,mlist[[s]]],c(1,2,3,4),sum)/apply(cyc_raindays[i+1,,,,,mlist[[s]]],c(1,2,3,4),sum),
                  apply(cyc_threshdays[i+1,,,,,mlist[[s]]],c(1,2,3,4),sum),
                  apply(cyc_threshdays2[i+1,,,,,mlist[[s]]],c(1,2,3,4),sum))
      
      for(l in 1:length(latreg))
        for(r in 1:3)
        {
          for(k in 1:4) cycchange[l,i,k,j,r,s-4]=sum(tmp1[r,k,,latlist[[l]]])
          cycchange[l,i,5,j,r,s-4]=sum(tmp1[r,,,latlist[[l]]])
        }
    }


cycchange[,,,,4,]=100*((cycchange[,,,,2,]/cycchange[,,,,1,])-1)
cycchange[,,,,5,]=100*((cycchange[,,,,3,]/cycchange[,,,,1,])-1)