load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; Code to edit - set up source and output directory, thresholds and years of interest

basedir="/g/data/al33/replicas/CMIP5/combined/"
odir="/g/data/eg3/asp561/CMIP5/"
syear=1979
eyear=2005
thresh=99.7
model="HadGEM2-CC"
member="r1i1p1"
cdir="/g/data/al33/replicas/CMIP5/combined/MOHC/"+model+"/historical/day/atmos/day/"+member+"/v20110930/pr/"

; Input and output files - assumes the easiest way is loading all daily precip data to handle different name structures
ifile=addfiles(systemfunc (" ls -1 " + cdir+"pr_day*.nc"),"r")
ofile=odir+model+"_pr_percentile_"+thresh+".nc"

lat=ifile[0]->lat
lon=ifile[0]->lon
time=cd_calendar(ifile[:]->time,0)
I=ind(time(:,0).ge.syear.and.time(:,0).le.eyear)

; Set up output array
tp_thresh=new((/dimsizes(lat),dimsizes(lon)/),"float")
tp_thresh!0="lat"
tp_thresh&lat=lat
tp_thresh!1="lon"
tp_thresh&lon=lon
tp_thresh@standard_name="percentile"
tp_thresh@units="mm"
tp_thresh@long_name = "Value of "+thresh+"th percentile of daily pr"

; Loop over all grid cells to calculate percentile value
do i=0,dimsizes(lat)-1
 do j=0,dimsizes(lon)-1
  tmp=ifile[:]->pr(I,i,j)
  tmp2=tmp(ind(.not.ismissing(tmp)))*60*60*24
  qsort(tmp2)
  tp_thresh(i,j)=(/tmp2(round(thresh*dimsizes(tmp2)/100,3))/)
  delete([/tmp,tmp2/])
 end do
end do

; Output the data as a netCDF

system("/bin/rm -f "+ofile)
ncdf=addfile(ofile,"c")
fAtt               = True            ; assign file attributes
fAtt@title         = model+" "+thresh+"th percentile for daily pr, "+syear+"-"+eyear
fAtt@source        = "CMIP5 daily pr data on NCI in "+cdir
fAtt@author        = "Acacia Pepler <acacia.pepler@bom.gov.au>"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

ncdf->percentile_values=tp_thresh

delete([/ncdf,ofile,tp_thresh,ifile,lat,lon,time,I/])


