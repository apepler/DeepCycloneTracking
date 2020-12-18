load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; Code to edit - set up source and output directory, thresholds and years of interest

model="HadGEM2-CC"
member="r1i1p1"
rcp="historical"
syear=1979
eyear=2005

pdir="/g/data/al33/replicas/CMIP5/combined/MOHC/"+model+"/historical/day/atmos/day/"+member+"/v20110930/pr/" ; CMIP5 rain data
edir="/g/data/eg3/asp561/CycloneTracking/CMIP5/"+model+"/"+rcp+"/"member+"/proj100_lows_rad2cv0.5/" ; Directory with gridded cyclone data, generated in step 2
threshdir="/g/data/eg3/asp561/CMIP5/" ; Directory with percentile values, generated in step 3
thresh=99 ; Percentile threshold for extreme rain days
thresh2=1 ; mm threshold for rain days

; Load the percentile values 

Data=addfile(threshdir+model+"_pr_percentile_"+thresh+".nc","r")
lat=Data->lat
lon=Data->lon
pctile=Data->percentile_values
delete(Data)

; Loop over years

do year=syear,eyear

; Define input and output files
efile=edir+model+"_"+rcp+"_"+rlist(dd)+"_deepcyclones_dailygrid_cv1var_dist500event_"+year+"_10deg_v2.nc"
ofile=edir+model+"_"+rcp+"_"+rlist(dd)+"_deepcyclones_dailygrid_cv1var_dist500event_"+year+"_rain_Q"+thresh+"_10deg_monthly.nc"
ofile2=edir+model+"_"+rcp+"_"+rlist(dd)+"_deepcyclones_dailygrid_cv1var_dist500event_"+year+"_rain_Above"+thresh2+"mm_10deg_monthly.nc"

; Define data arrays

surf_rain=new((/12,dimsizes(lat),dimsizes(lon),4/),"float")
surf_rain=0
surf_rain!0="time"
surf_rain!1="lat"
surf_rain&lat=lat
surf_rain!2="lon"
surf_rain&lon=lon
surf_rain!3="type"

surf_threshdays=surf_rain
cyc_freq=surf_rain

surf_rain@standard_name="precipitation"
surf_rain@units="mm"
surf_threshdays@standard_name="count"
surf_threshdays@units="days"
cyc_freq@standard_name="count"
cyc_freq@units="days"

surf_rain2=surf_rain
surf_threshdays2=surf_threshdays

; Load all precip data

Data=addfiles(systemfunc (" ls -1 " + pdir+"pr_day*.nc"),"r")
time=cd_calendar(Data[:]->time,0)

; Load all cyclone data

Data2=addfile(efile,"r")
time2=cd_calendar(Data2->time,0)
len=dimsizes(time2)
surf=Data2->cycfreq


; Loop over the 12 months and the 4 cyclone categories to calculate the total number of cyclone days

do m=1,12
do t=0,3
 M=ind(time2(:,1).eq.m)
 cyc_freq(m-1,:,:,t)=(/dim_sum_n_Wrap(where(surf(M,:,:).eq.t,1,0),0)/)
 delete(M)
end do
end do

; Loop over all days to load daily rainfall, compare it to thresholds & cyclone category data, and add to the appropriate array

print("Loaded")
do i=0,len(0)-1
  I=ind(time(:,0).eq.time2(i,0).and.time(:,1).eq.time2(i,1).and.time(:,2).eq.time2(i,2))

; Only run if the data exists - to deal with leap years and other oddities
 if(any(ismissing(I))) then
  print("Missing day")
  print(time2(i,0:2))
 else
  m=tointeger(time2(i,1)-1)
  tmp1=Data[:]->pr(I,:,:)*60*60*24
  tmp2=where(tmp1.ge.thresh2,1,0)

  do t=0,3
    surf_rain(m,:,:,t)=surf_rain(m,:,:,t)+tmp1*where(surf(i,:,:).eq.t,1,0) ; Total rain
    surf_threshdays(m,:,:,t)=surf_threshdays(m,:,:,t)+where(tmp1.ge.pctile,1,0)*where(surf(i,:,:).eq.t,1,0) ; Days above percentile
    surf_threshdays2(m,:,:,t)=surf_threshdays2(m,:,:,t)+tmp2*where(surf(i,:,:).eq.t,1,0) ; Rain days above thresh2
	surf_rain2(m,:,:,t)=surf_rain2(m,:,:,t)+tmp1*tmp2*where(surf(i,:,:).eq.t,1,0) ; Total rain from days above thresh2
  end do
  delete(tmp1)
 end if
 delete(I)
end do

; Output netCDF data

cyc_freq@long_name="Annual "+model+" days with cyclones present within a 10 degree radius. None (1), Deep (2), Shallow surface (3) or Shallow upper (4)"
 surf_rain@long_name="Annual "+model+" rainfall based on cyclones present within a 10 degree radius. None (1), Deep (2), Shallow surface (3) or Shallow upper (4)"
 surf_threshdays@long_name="Annual "+model+" days above "+thresh+"th percentile based on cyclones present within a 10 degree radius. None (1), Deep (2), Shallow surface (3) or Shallow upper (4)"

system("/bin/rm -f "+ofile)   ; remove any pre-existing file
ncdf = addfile(ofile ,"c")  ; open output netCDF file
ncdf->cyc_rain=surf_rain
ncdf->cyc_threshdays=surf_threshdays
ncdf->cyc_freq=cyc_freq

 surf_rain2@long_name="Monthly "+model+" rainfall on rain days >="+thresh2+"mm based on cyclones present within a 10 degree radius. None (1), Deep (2), Shallow surface (3) or Shallow upper (4)"
 surf_threshdays2@long_name="Monthly "+model+" days above "+thresh2+"mm based on cyclones present within a 10 degree radius. None (1), Deep (2), Shallow surface (3) or Shallow upper (4)"

system("/bin/rm -f "+ofile2)   ; remove any pre-existing file
ncdf2 = addfile(ofile2 ,"c")  ; open output netCDF file
ncdf2->cyc_rain=surf_rain2
ncdf2->cyc_threshdays=surf_threshdays2
ncdf2->cyc_freq=cyc_freq

delete([/ncdf,surf_rain,surf_threshdays,cyc_freq,Data,Data2,time,time2,surf,tmp2,surf_rain2,surf_threshdays2,ncdf2/])

end do
