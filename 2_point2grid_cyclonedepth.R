# This version has three categories
# Any type of deep cyclone, shallow surface only, shallow upper only

library(raster)
library(sp)
library(ncdf4)
library(abind)
library(geosphere)

spreadeffect<-function(data,winwid=3,circ=T,lon=NaN,lat=NaN)
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
  m3=focal(m2,w2,pad=T,padValue=0)
  tmp=t(as.matrix(m3)[ll,(winwid+1):(a[1]+winwid)])
  tmp[tmp>0]=1
  return(tmp)
}


storm_combo<-function(year,outf=NA,sdir,udir,outdir=sdir,cvthresh=c(1,8),winwid=12,move=NaN,dur=NaN,closedS=F,closedU=F,event=T,dist=1000,lats=seq(-89.5,89.5,1),lons=seq(0.5,359.5,1))
{
  
  
  datelist=seq.Date(as.Date(paste0(min(year),"-01-01")), as.Date(paste0(max(year),"-12-31")),by="1 day")
  datelist=data.frame(Date=datelist,Year=as.numeric(format(datelist,"%Y")),Month=as.numeric(format(datelist,"%m")),YYYYMMDD=as.numeric(format(datelist,"%Y%m%d")))
  
  fixesS=read.table(paste0(sdir,"/tracks_",year,".dat"), sep="",skip=0)
  fixesU=read.table(paste0(udir,"/tracks_",year,".dat"), sep="",skip=0)
  
  
  
  print("Loaded cyclones") 
  colnames(fixesS)<-colnames(fixesU)<-c("ID","Fix","Date","Time","Open",
                                        "Lon","Lat","MSLP","CV","Depth","Radius","Up","Vp")
  ##Clean up dates
  yy=floor(fixesS$Date/10000)
  yy2=unique(yy)
  if(length(yy2)>1) fixesS=fixesS[yy==yy2[2],]
  yy=floor(fixesU$Date/10000)
  yy2=unique(yy)
  if(length(yy2)>1) fixesU=fixesU[yy==yy2[2],]
  
  fixesS$Date=(fixesS$Date%%10000) + year*10000
  fixesU$Date=(fixesU$Date%%10000) + year*10000
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
    
    x<-rle(fixesU$ID)
    I=which(fixesU$ID%in%sort(unique(x$values[x$lengths>=dur/4])))
    fixesU=fixesU[I,]
  }
  
  # Optional: movement threshold
  
  if(!is.na(move))
  {
    fixesS$Move<-NaN
    I=which(fixesS$Fix>1)
    if(I[1]==1) I=I[-1]
    for(i in 1:length(I)) fixesS$Move[I[i]]=spDistsN1(cbind(fixesS$Lon[I[i]],fixes$SLat[I[i]]),cbind(fixesS$Lon[I[i]-1],fixesS$Lat[I[i]-1]),longlat=T)
    
    x<-rle(fixesS$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Move=rep(0,length(x$values)))
    for(i in 1:length(events$ID)) events$Move[i]=spDistsN1(cbind(fixesS$Lon[min(I)],fixesS$Lat[min(I)]),
                                                           cbind(fixesS$Lon[max(I)],fixesS$Lat[max(I)]),longlat=T)
    I=which(fixesS$ID%in%sort(unique(x$values[events$Move>=move])))
    fixesS=fixesS[I,]
    
    fixesU$Move<-NaN
    I=which(fixesU$Fix>1)
    if(I[1]==1) I=I[-1]
    for(i in 1:length(I)) fixesU$Move[I[i]]=spDistsN1(cbind(fixesU$Lon[I[i]],fixes$SLat[I[i]]),cbind(fixesU$Lon[I[i]-1],fixesU$Lat[I[i]-1]),longlat=T)
    
    x<-rle(fixesU$ID)
    events<-data.frame(ID=x$values,Length=x$lengths,Move=rep(0,length(x$values)))
    for(i in 1:length(events$ID)) events$Move[i]=spDistsN1(cbind(fixesU$Lon[min(I)],fixesU$Lat[min(I)]),
                                                           cbind(fixesU$Lon[max(I)],fixesU$Lat[max(I)]),longlat=T)
    I=which(fixesU$ID%in%sort(unique(x$values[events$Move>=move])))
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
  fixesS=fixesS[fixesS$CV>=cvthresh[1],]
  fixesU=fixesU[fixesU$CV>=cvthresh[2],]
  
  ##### OKAY! Time to do some matching
  #### As before - match by same day in UTC for version 1
  #### Then look at the event-max CV for version2
  
  fixesS$UpperCV2<-fixesS$UpperCV<-0
  
  for(j in 1:length(fixesS[,1]))
  {
    
    J=which(fixesU$Date==fixesS$Date[j]  & sign(fixesU$Lat)==sign(fixesS$Lat[j])) ## Remove all that are further N/S than possible, with a leeway
    
    if(length(J)>0) 
    {
      tmp=distGeo(cbind(fixesU$Lon2[J],fixesU$Lat[J]),cbind(fixesS$Lon2[j],fixesS$Lat[j]))/1000
      K=which(tmp<dist)
      if(length(K)>0) fixesS$UpperCV[j]=max(fixesU$CV[J[K]])
    }
  }
  
  fixesU$LowerCV2<-fixesU$LowerCV<-0
  
  for(j in 1:length(fixesU[,1]))
  {
    J=which(fixesS$Date==fixesU$Date[j]  & sign(fixesS$Lat)==sign(fixesU$Lat[j]))
    if(length(J)>0) {
      tmp=distGeo(cbind(fixesS$Lon2[J],fixesS$Lat[J]),cbind(fixesU$Lon2[j],fixesU$Lat[j]))/1000
      K=which(tmp<dist)
      if(length(K)>0) fixesU$LowerCV[j]=max(fixesS$CV[J[K]])
    }
  }
  
  ## Event-based approach
  if(event)
  {
    IDs=unique(fixesS$ID)
    for(id in IDs)
    {
      I=which(fixesS$ID==id)
      J=which(fixesS$ID==id & !is.na(fixesS$UpperCV))
      if(length(J)>0) fixesS$UpperCV2[I]=max(fixesS$UpperCV[J],na.rm=T)
    }
    
    IDs=unique(fixesU$ID)
    for(id in IDs)
    {
      I=which(fixesU$ID==id)
      J=which(fixesU$ID==id & !is.na(fixesU$LowerCV))
      if(length(J)>0) fixesU$LowerCV2[I]=max(fixesU$LowerCV[J],na.rm=T)
    }
  }
  print("Completed mtaching")
  ### Now I need it to make the grid, but first I need to set the output lat and lon resolutionn
  ## So I'm gonna make a surface and an upper grid - no cyclone is 0, deep is 1, shallow is 2 - so deep overrides shallow.
  ## Or is it better to make deep+shallow as 3, so I can check frequency? Maybe that's safer
  
  #dlat=abs(lats[2]-lats[1])
  #dlon=abs(lons[2]-lons[1])
  #lat2=lats/dlat
  #lon2=lons/dlon
  
  lat2=seq(1,length(lats))
  lon2=seq(1,length(lons))
  
  fixesS$Lat2=sapply(fixesS$Lat,function(x)which.min(abs(x - lats)))
  fixesS$Lon2=sapply(fixesS$Lon,function(x)which.min(abs(x - lons)))
  fixesU$Lat2=sapply(fixesU$Lat,function(x)which.min(abs(x - lats)))
  fixesU$Lon2=sapply(fixesU$Lon,function(x)which.min(abs(x - lons)))
  
  #fixesS$Lat2=round(fixesS$Lat/dlat,0)
  #fixesS$Lon2=fixesS$Lon%%360
  #if(min(lons)<0) fixesS$Lon2[fixesS$Lon2>180]=fixesS$Lon2[fixesS$Lon2>180]-360
  #fixesS$Lon2=round(fixesS$Lon2/dlon,0)
  #fixesU$Lat2=round(fixesU$Lat/dlat,0)
  #fixesU$Lon2=fixesU$Lon%%360
  #if(min(lons)<0) fixesU$Lon2[fixesU$Lon2>180]=fixesU$Lon2[fixesU$Lon2>180]-360
  #fixesU$Lon2=round(fixesU$Lon2/dlon,0)
  
  if(event) I=which(fixesS$UpperCV2>=cvthresh[2]) else I=which(fixesS$UpperCV>=cvthresh[2])
  grid1=table(factor(fixesS$Lon2[I],levels=lon2),factor(fixesS$Lat2[I],levels=lat2),factor(fixesS$Date[I],levels=datelist$YYYYMMDD))
  for(i in 1:length(datelist[,1])) grid1[,,i]=spreadeffect(grid1[,,i],lon=lons,lat=lats,circ=T,winwid=winwid)
  
  grid2=table(factor(fixesS$Lon2[-I],levels=lon2),factor(fixesS$Lat2[-I],levels=lat2),factor(fixesS$Date[-I],levels=datelist$YYYYMMDD))
  for(i in 1:length(datelist[,1])) grid2[,,i]=spreadeffect(grid2[,,i],lon=lons,lat=lats,circ=T,winwid=winwid)
  
  if(event) I=which(fixesU$LowerCV2>=cvthresh[1]) else I=which(fixesU$LowerCV>=cvthresh[1])
  grid3=table(factor(fixesU$Lon2[I],levels=lon2),factor(fixesU$Lat2[I],levels=lat2),factor(fixesU$Date[I],levels=datelist$YYYYMMDD))
  for(i in 1:length(datelist[,1])) grid3[,,i]=spreadeffect(grid3[,,i],lon=lons,lat=lats,circ=T,winwid=winwid)
  
  grid4=table(factor(fixesU$Lon2[-I],levels=lon2),factor(fixesU$Lat2[-I],levels=lat2),factor(fixesU$Date[-I],levels=datelist$YYYYMMDD))
  for(i in 1:length(datelist[,1])) grid4[,,i]=spreadeffect(grid4[,,i],lon=lons,lat=lats,circ=T,winwid=winwid)
  
  cycfreq=array(0,dim(grid1))
  cycfreq[grid1==1 | grid3==1]=1
  cycfreq[cycfreq==0 & grid2==1]=2
  cycfreq[cycfreq==0 & grid4==1]=3
  
  print("Completed gridding")
  
  ##
  dimX<-ncdim_def("lon","degrees_E",lons)
  dimY<-ncdim_def("lat","degrees_N",lats)
  
  tmp=seq.POSIXt(as.POSIXct(paste0(year,"0101 09:00"),format="%Y%m%d %H:%M",tz="GMT"),as.POSIXct(paste0(year,"1231 09:00"),format="%Y%m%d %H:%M",tz="GMT"),by="24 hours")
  dimT<-ncdim_def("time","hours since 1970-1-1 00:00:00",as.numeric(tmp)/(60*60))  
  fillvalue <- 9999
  cyc_def <- ncvar_def("cycfreq","count",list(dimX,dimY,dimT),fillvalue,"Cyclones from ERA5 - Deep (1), Shallow surface (2) or Shallow upper (3)",prec="integer")
  
  ncout <- nc_create(paste0(outdir,"/",outf),cyc_def) #force_v4=T)
  
  # put variables
  ncvar_put(ncout,cyc_def,cycfreq)
  
  # put additional attributes into dimension and data variables
  ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
  ncatt_put(ncout,"lat","axis","Y")
  ncatt_put(ncout,"time","axis","T")
  nc_close(ncout)
}

#Run
storm_combo(year, #Year to run code over
            winwid=10, # Radius of influence for cyclones (degrees)
            dist=500,event=T, # Distance (km) between surface and upper lows to be "deep", and whether to use an event approach
            cvthresh=c(1,8), # Surface and upper minimum Laplacian thresholds
            dur=NaN,move=NaN, #Minimum suration (timesteps) or movement (km) for the event
            closedS=F,closedU=F, # Are Surface or Upper lows required to be closed?
            lats=elat,lons=elon, # Supply your lat/lon grid of choice
            sdir=sdir, # Directory where surface tracks are stored
            udir=udir, # DIrectory where upper tracks are stored
            outf=outf) # Name for output netCDF file
print(warnings())





