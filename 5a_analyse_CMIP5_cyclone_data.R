# THis has a set of code that I run interactively to analyse cyclone data

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

# A function that applies spatial smoothing to grids

makesmooth<-function(data,winwid=3,circ=T,lon=NaN,lat=NaN,spread=F,total=F)
{
  a=dim(data)
  
  if(is.na(lon[1]))
  {
    lon=seq(1,a[1])
    lat=seq(1,a[2])
  }
  
  if(lat[2]<lat[1]) ll=a[2]:1 else ll=1:a[2]
  
  ## Only make cyclic if it's actually global, -180 to 180 or 0 to 360
  if(lon[2]-lon[1]>=350) m1=abind(data[(a[1]-winwid+1):a[1],ll],data[,ll],data[1:winwid,ll],along=1) else
    m1=abind(matrix(0,winwid,a[2]),data[,ll],matrix(0,winwid,a[2]),along=1)
  
  m2=raster(t(m1),xmn=min(lon)-winwid,xmx=max(lon)+winwid,ymn=min(lat),ymx=max(lat))
  if(circ) w2=focalWeight(m2,winwid,type="circle") else w2=matrix(1,1+2*winwid,1+2*winwid)
  if(total) w2[w2>0]=1
  if(spread) m3=focal(m2,w2,pad=T,padValue=0) else m3=focal(m2,w2)
  tmp=t(as.matrix(m3)[ll,(winwid+1):(a[1]+winwid)])
  if(spread) tmp[tmp>0]=1
  return(tmp)
}

# Set up some of the base directories etc
basedir="/g/data/eg3/asp561/CycloneTracking/"

# Details about the names, member, thresholds etc for each model
cmiplist=c("ERA5","ACCESS1-0","ACCESS1-3","CCSM4","GFDL-ESM2M","MRI-CGCM3","CNRM-CM5","MIROC5","NorESM1-M","CanESM2","HadGEM2-CC")
memlist=c(NaN,"r1i1p1","r1i1p1","r6i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1","r1i1p1")
cmip=c(F,rep(T,10))
cvthresh=rbind(c(1.00,0.90,0.90,1.15,1.10,1.30,0.90,1.15,1.05,1.15,0.95), # Minimum Laplacian for surface lows
               c(8.1,6.2,5.9,6.0,5.7,6.8,5.3,7.0,5.2,6.1,6.4)) # Minimum laplacian for upper lows

#Common criteria o apply
dur=NaN # Minimum duration (timesteps)?
move=NaN # Minimum movement (km)?
closedS=F # Do surface lows need to be closed
closedU=F # Do upper lows need to be closed?

# These are some extra parameters I need to find the right cyclone directory - you will need different ones
thresh="rad2cv0.5"
threshD="rad2cv1"
proj=100
lev=500

## Years to group analysis over and make the file size more manageable
years=1979:2099
yeargroup=rep(0,length(years))
yeargroup[years%in%1979:2005]=1
yeargroup[years%in%2010:2039]=2
yeargroup[years%in%2040:2069]=3
yeargroup[years%in%2070:2099]=4
yeargrouplist=list(1979:2005,2010:2039,2040:2069,2070:2099)

# Details for each rco
rcp=c("historical","rcp45","rcp85")
ryears=c("19502005","20062100","20062100")

# Set up arrays
lat=seq(-89.5,89.5,1)
lon=seq(0.5,359.5,1)
surfreq=array(0,c(length(lon),length(lat),4,12,length(cmiplist),3))
dimnames(surfreq)[[3]]<-c("1979-2005","2010-2039","2040-2069","2070-2099")
dimnames(surfreq)[[5]]<-cmiplist
dimnames(surfreq)[[6]]<-rcp
upfreq<-surfreq

# Set up dates and seasons for use in plotting

datelist=seq.Date(as.Date(paste0(min(years),"-01-01")), as.Date(paste0(max(years),"-12-31")),by="1 day")
datelist=data.frame(Date=datelist,Year=as.numeric(format(datelist,"%Y")),Month=as.numeric(format(datelist,"%m")),YYYYMMDD=as.numeric(format(datelist,"%Y%m%d")))
snames=c("MAM","JJA","SON","DJF","MJJASO","NDJFMA","Annual","April-September","October-March")
mlist=list(3:5,6:8,9:11,c(12,1:2),5:10,c(11:12,1:4),1:12,4:9,c(10:12,1:3))


for(r in 1:3)
  for(i in 1:length(cmip))
  {
    
    ## Sets up the directories to find the surface and upper lows - this will need to change
    if(rcp[r]=="rcp45" & cmiplist[i]=="MIROC5") next
    if(!cmip[i] & r>1) next
    
    if(r==1 & !cmip[i])
    {
      udir=paste0(basedir,cmiplist[i],"/",lev,"hPa_z/daily_proj",proj,"_lows_",threshD,"_global")
      sdir=paste0(basedir,cmiplist[i],"/proj",proj,"_lows_",thresh,"_open")
      ylist=1979:2005
    } else {
      udir=paste0(basedir,"CMIP5/",cmiplist[i],"/",rcp[r],"/",memlist[i],"/",lev,"hPa_z/daily_proj",proj,"_lows_",threshD)
      sdir=paste0(basedir,"CMIP5/",cmiplist[i],"/",rcp[r],"/",memlist[i],"/proj",proj,"_lows_",thresh)
      if(r==1) ylist=1979:2005 else ylist=2010:2099
    }
    
    Y=which(years%in%ylist)
    
    # Loops over all years to load tracks and add to the correct yeargroup in the array
    
    for(y in Y)
    {
      print(years[y])
      fixesS=read.table(paste0(sdir,"/tracks_",years[y],".dat"), sep="",skip=0)
      fixesU=read.table(paste0(udir,"/tracks_",years[y],".dat"), sep="",skip=0)
      
      colnames(fixesS)<-colnames(fixesU)<-c("ID","Fix","Date","Time","Open",
                                            "Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
      
      ##Clean up dates & remove data for preceding December
      yy=floor(fixesS$Date/10000)
      yy2=unique(yy)
      if(length(yy2)>1) fixesS=fixesS[yy==yy2[2],]
      yy=floor(fixesU$Date/10000)
      yy2=unique(yy)
      if(length(yy2)>1) fixesU=fixesU[yy==yy2[2],]
      
      fixesS$Date=(fixesS$Date%%10000) + years[y]*10000
      fixesU$Date=(fixesU$Date%%10000) + years[y]*10000
      fixesS$Month=floor(fixesS$Date/100)%%100
      fixesU$Month=floor(fixesU$Date/100)%%100
      
      ##Clean up longitudes
      fixesS$Lon2=fixesS$Lon
      fixesS$Lon2[fixesS$Lon>180]=fixesS$Lon[fixesS$Lon>180]-360
      fixesU$Lon2=fixesU$Lon
      fixesU$Lon2[fixesU$Lon>180]=fixesU$Lon[fixesU$Lon>180]-360
      
      ##Optional - duration threshold
      if(!is.na(dur)) 
      {
        x<-rle(fixesS$ID)
        I=which(fixesS$ID%in%sort(unique(x$values[x$lengths>=dur])))
        fixesS=fixesS[I,]
        
        # Only affect curface cyclones - using daily data means upper already have an effective duration
        x<-rle(fixesU$ID)
        I=which(fixesU$ID%in%sort(unique(x$values[x$lengths>=dur/4])))
        fixesU=fixesU[I,]
      }
      
      # Optional: movement threshold
      
      if(!is.na(move))
      {
        x<-rle(fixesS$ID)
        events<-data.frame(ID=x$values,Length=x$lengths,Move=rep(0,length(x$values)))
        for(z in 1:length(events$ID))
        {
          I=which(fixesS$ID==events$ID[z])
          events$Move[z]=distGeo(cbind(fixesS$Lon2[min(I)],fixesS$Lat[min(I)]),cbind(fixesS$Lon2[max(I)],fixesS$Lat[max(I)]))/1000
        }
        
        I=which(fixesS$ID%in%sort(unique(events$ID[events$Move>=move])))
        fixesS=fixesS[I,]
        
        x<-rle(fixesU$ID)
        events<-data.frame(ID=x$values,Length=x$lengths,Move=rep(0,length(x$values)))
        for(z in 1:length(events$ID))
        {
          I=which(fixesU$ID==events$ID[z])
          events$Move[z]=distGeo(cbind(fixesU$Lon2[min(I)],fixesU$Lat[min(I)]),cbind(fixesU$Lon2[max(I)],fixesU$Lat[max(I)]))/1000
        }
        
        I=which(fixesU$ID%in%sort(unique(events$ID[events$Move>=move])))
        fixesU=fixesU[I,]
      }
      
      ##Optional: only closed systems
      if(closedS)
      {
        I=which(fixesS$Open%in%c(0,10))
        fixesS=fixesS[I,]
      }
      if(closedU)
      {
        I=which(fixesU$Open%in%c(0,10))
        fixesU=fixesU[I,]
      }
      
      ##Find the correct threshold
      fixesS=fixesS[fixesS$CV>=cvthresh[1,i],]
      fixesU=fixesU[fixesU$CV>=cvthresh[2,i],]
      
     
      #### Aaaand finally, we make the table of frequency!
      fixesS$Lat2=floor(fixesS$Lat)
      fixesS$Lon2=(floor(fixesS$Lon))%%360
      fixesU$Lat2=floor(fixesU$Lat)
      fixesU$Lon2=(floor(fixesU$Lon))%%360
      print(warnings()) 
      yg=yeargroup[y] 
      surfreq[,,yg,,i,r]=surfreq[,,yg,,i,r]+table(factor(fixesS$Lon2,levels=floor(lon)),factor(fixesS$Lat2,levels=floor(lat)),factor(fixesS$Month,levels=1:12))
      upfreq[,,yg,,i,r]=upfreq[,,yg,,i,r]+table(factor(fixesU$Lon2,levels=floor(lon)),factor(fixesU$Lat2,levels=floor(lat)),factor(fixesU$Month,levels=1:12))
    }
  }

# Convert to annual frequency

for(y in 1:4)
{
  upfreq[,,y,,,]=upfreq[,,y,,,]/length(yeargrouplist[[y]])
  surfreq[,,y,,,]=surfreq[,,y,,,]/length(yeargrouplist[[y]])
}

# Clean up data we expect to be missing, e.g. historical data for rcp4.5/rcp8.5
upfreq[,,2:4,,1,]<-surfreq[,,2:4,,1,]<-NaN
upfreq[,,2:4,,2:length(cmiplist),1]<-surfreq[,,2:4,,2:length(cmiplist),1]<-NaN
upfreq[,,1,,2:length(cmiplist),2:3]<-surfreq[,,1,,2:length(cmiplist),2:3]<-NaN
upfreq[,,,,8,2]<-surfreq[,,,,8,2]<-NaN # Missing MIROC5 rcp4.5


### Supplementary Figure S5: cyclone frequency (%) in 1979-2005 for ERAI and 10 model mean, both surface and upper

breaks=c(seq(0,0.3,0.025),100)
cols=rev(viridis(length(breaks)-1))
breaksa=c("","","","","0.1%","","","","0.2%","","","","0.3%","")

breaks2=c(-100,seq(-0.2,0.2,0.025),100)
breaks2a=c(-100,paste0(seq(-0.2,-0.025,0.025),"%"),paste0("+",seq(0,0.2,0.025),"%"),1000)
cols2=colorRampPalette(brewer.pal(11,"RdBu"))(length(breaks2)-1)

for(s in 7) # Annual
{

  S=length(which(datelist$Month%in%mlist[[s]] & datelist$Year%in%yeargrouplist[[2]]))/length(yeargrouplist[[2]])  # Average number of days p.a. in this season
  pdf(file=paste0("FigureS5_CMIP5_global_cycfreq_cv1var_dist500_19792005_surfvsupper_",snames[s],"_smooth5mean_biasabs_fixed_10models.pdf"),width=16,height=6,pointsize=16)
  par(mar=c(3,3,3,1))
  layout(rbind(c(1,2,7,3,8),c(4,5,7,6,8)),width=c(1,1,0.35,1,0.35))
  
  tmp1=apply(upfreq[,,1,mlist[[s]],,1],c(1,2,4),sum)
  for(i in 1:length(cmiplist)) tmp1[,,i]=makesmooth(tmp1[,,i],lon=lon,lat=lat,win=5,total=F)
  
  image(lon,lat,100*tmp1[,,1]/S,
        breaks=breaks,col=cols,main="a) Upper cyclones: ERA5",xlab="",ylab="",cex.main=1.5)
  map("world2",add=T)
  
  image(lon,lat,100*apply(tmp1[,,2:length(cmiplist)],c(1,2),mean)/S,
        breaks=breaks,col=cols,main="b) Upper cyclones: CMIP5",xlab="",ylab="",cex.main=1.5)
  map("world2",add=T)
  
  tmp2=tmp1[,,2:length(cmiplist)]
  for(i in 1:10) tmp2[,,i]=tmp2[,,i]-tmp1[,,1]
  image(lon,lat,100*apply(tmp2,c(1,2),mean)/S,
        breaks=breaks2,col=cols2,main="c) Upper cyclone bias",xlab="",ylab="",cex.main=1.5)
  map("world2",add=T)
  
  tmp1=apply(surfreq[,,1,mlist[[s]],,1],c(1,2,4),sum)
  for(i in 1:length(cmiplist)) tmp1[,,i]=makesmooth(tmp1[,,i],lon=lon,lat=lat,win=5,total=F)
  
  image(lon,lat,100*tmp1[,,1]/(4*S),
        breaks=breaks,col=cols,main="d) Surface cyclones: ERA5",xlab="",ylab="",cex.main=1.5)
  map("world2",add=T)
  
  image(lon,lat,100*apply(tmp1[,,2:length(cmiplist)],c(1,2),mean)/(4*S),
        breaks=breaks,col=cols,main="e) Surface cyclones: CMIP5",xlab="",ylab="",cex.main=1.5)
  map("world2",add=T)
  
  tmp2=tmp1[,,2:length(cmiplist)]
  for(i in 1:10) tmp2[,,i]=tmp2[,,i]-tmp1[,,1]
  image(lon,lat,100*apply(tmp2,c(1,2),mean)/(4*S),
        breaks=breaks2,col=cols2,main="f) Surface cyclone bias",xlab="",ylab="",cex.main=1.5)
  map("world2",add=T)
  
  ColorBar(brks=breaksa,cols=cols,vert=T)
  ColorBar(brks=breaks2a,cols=cols2,vert=T,subsampleg=4)
  dev.off()
}

### Supplementary Figure S3 - change between 1979-2005 and 2070-2099 - upper and lower

breaks2=c(-100000,seq(-50,50,10),1000000)
breaks2a=c(-100,paste0(seq(-50,-10,10),"%"),paste0("+",seq(0,50,10),"%"),1000)
cols2=colorRampPalette(brewer.pal(11,"RdBu"))(length(breaks2)-1)

for(s in 7) # Annual
  for(g in 4) # 2070-2099
{
  pdf(file=paste0("FigureS3_CMIP5_global_cycchangepc_cv1var_dist500_surfvsupper_",snames[s],"_smooth5_19792005_",ygname[g],"_rcp85_points_10models.pdf"),width=10,height=3)
  par(mar=c(2,2,3,1))
  layout(cbind(1,2,3),width=c(1,1,0.25))
  
  # Add 5 degrees of spatial smoothing so plots are less messy, and remove areas where historical frequency is very low
  tmp1=makesmooth(apply(surfreq[,,1,mlist[[s]],2:length(cmiplist),1],c(1,2),sum),lon=lon,lat=lat,win=5)/10 # Sum over months and models then divide by number of models
  tmp1[tmp1<0.01]=NaN
  tmp2=makesmooth(apply(surfreq[,,g,mlist[[s]],2:length(cmiplist),3],c(1,2),sum),lon=lon,lat=lat,win=5)/10
  
  image(lon,lat,100*((tmp2/tmp1)-1), # Percentage change
        breaks=breaks2,col=cols2,main="a) Change in surface cyclone frequency",xlab="",ylab="",cex.main=1.5)
  map("world2",add=T)

  # Count how many models have increasing or decreasing trends for mask
  postrend<-negtrend<-array(0,dim(tmp1))
  I=which(is.na(tmp1))
  postrend[I]=NaN
  negtrend[I]=NaN
  for(i in 2:length(cmiplist))
  {
    tmp1=makesmooth(apply(surfreq[,,1,mlist[[s]],i,1],c(1,2),sum),lon=lon,lat=lat,win=5)
    tmp2=makesmooth(apply(surfreq[,,g,mlist[[s]],i,3],c(1,2),sum),lon=lon,lat=lat,win=5)/tmp1
    I=which(!is.na(tmp2) & tmp2>1 & tmp1>=0.01)
    postrend[I]=postrend[I]+1
    I=which(!is.na(tmp2) & tmp2<1 & tmp1>=0.01)
    negtrend[I]=negtrend[I]+1
  }

  I=which(!is.na(tmp1) & postrend<8 & negtrend<8,arr.ind=T)
  points(lon[I[,1]],lat[I[,2]],cex=0.1,pch=16)

  # Repwat for upper cyclones
  tmp1=makesmooth(apply(upfreq[,,1,mlist[[s]],2:length(cmiplist),1],c(1,2),sum),lon=lon,lat=lat,win=5)/10
  tmp1[tmp1<0.01]=NaN
  tmp2=makesmooth(apply(upfreq[,,g,mlist[[s]],2:length(cmiplist),3],c(1,2),sum),lon=lon,lat=lat,win=5)/10
  
  image(lon,lat,100*((tmp2/tmp1)-1),
        breaks=breaks2,col=cols2,main="b) Change in upper cyclone frequency",xlab="",ylab="",cex.main=1.5)
  map("world2",add=T)

  postrend<-negtrend<-array(0,dim(tmp1))
  I=which(is.na(tmp1))
  postrend[I]=NaN
  negtrend[I]=NaN
  for(i in 2:length(cmiplist))
  {
    tmp1=makesmooth(apply(upfreq[,,1,mlist[[s]],i,1],c(1,2),sum),lon=lon,lat=lat,win=5)
    tmp2=makesmooth(apply(upfreq[,,g,mlist[[s]],i,3],c(1,2),sum),lon=lon,lat=lat,win=5)/tmp1
    I=which(!is.na(tmp2) & tmp2>1 & tmp1>=0.01)
    postrend[I]=postrend[I]+1
    I=which(!is.na(tmp2) & tmp2<1 & tmp1>=0.01)
    negtrend[I]=negtrend[I]+1
  }
  I=which(!is.na(tmp1) & postrend<8 & negtrend<8,arr.ind=T)
  points(lon[I[,1]],lat[I[,2]],cex=0.1,pch=16)

  ColorBar(brks=breaks2a,cols=cols2,vert=T)
  dev.off()
  }
