; ################################################################
;
; The purpose of this program is to read in the WRF RCM data and
; create composites of fields based on SOM analysis.
;
; ################################################################
pro composite,num_rows,num_cols,expID
  
; Configuration
;expID      = 'WRF_hadgem_50km'
SOMexpID   = 'ONDJFM'
dirIN      = '/Projects/HydroMet/dswales/NA-CORDEX/'+expID+'/'
;num_rows   = 2
;num_cols   = 2

; ################################################################
; Read in SOM BMU information
; ################################################################
somBMUfile = '/data/dswales/NA-CORDEX/SOMs/bmus/BMU.'+string(num_rows,format='(i1)')+$
             'x'+string(num_cols,format='(i1)')+'.'+SOMexpID+'.'+expID+'.txt'
dumb=''
openr,101,somBMUfile
readf,101,dumb
readf,101,ntime
yearSOM  = intarr(ntime)
monthSOM = intarr(ntime)
daySOM   = intarr(ntime)
bmusSOM  = intarr(ntime)
readf,101,yearSOM
readf,101,dumb
readf,101,ntime
readf,101,monthSOM
readf,101,dumb
readf,101,ntime
readf,101,daySOM
readf,101,dumb
readf,101,ntime
readf,101,bmusSOM
close,101

; ################################################################
; Read in WRF climatology data.
; ################################################################
fileClim = '/home/dswales/Projects/NA-CORDEX/data/clim/'+expID+'/ivt.mon.clim.nc'
fileID   = ncdf_open(fileClim)
ncdf_varget,fileID,ncdf_varid(fileID,'mean'),ivtMclim
ncdf_varget,fileID,ncdf_varid(fileID,'stdev'),ivtSclim
ncdf_close,fileID
fileClim = '/home/dswales/Projects/NA-CORDEX/data/clim/'+expID+'/precip.mon.clim.nc'
fileID   = ncdf_open(fileClim)
ncdf_varget,fileID,ncdf_varid(fileID,'mean'),precipMclim
ncdf_close,fileID

; ################################################################
; Read in WRF RCM data
; ################################################################
files01 = file_search(dirIN,'wrfout_*01.nc')
files02 = file_search(dirIN,'wrfout_*02.nc')
files03 = file_search(dirIN,'wrfout_*03.nc')
files10 = file_search(dirIN,'wrfout_*10.nc')
files11 = file_search(dirIN,'wrfout_*11.nc')
files12 = file_search(dirIN,'wrfout_*12.nc')
files   = [files01,files02,files03,files10,files11,files12]

nDay = 0
for ij=0,n_elements(files)-1 do begin
   print,ij,' of ~',n_elements(files)

   fileID = ncdf_open(files(ij))
   if (ij eq 0) then begin
      ncdf_varget,fileID,ncdf_varid(fileID,'XLONG'),lon
      ncdf_varget,fileID,ncdf_varid(fileID,'XLAT'),lat
      nlon = n_elements(lon(*,0))
      nlat = n_elements(lat(0,*))
      ivtSOM     = fltarr(nLon,nLat,num_rows*num_cols)
      ivtASOM    = fltarr(nLon,nLat,num_rows*num_cols)
      precipSOM  = fltarr(nLon,nLat,num_rows*num_cols)
      precipASOM = fltarr(nLon,nLat,num_rows*num_cols)
      countSOM   = fltarr(num_rows*num_cols)
   endif
   ncdf_varget,fileID,ncdf_varid(fileID,'IVTU'),ivtU
   ncdf_varget,fileID,ncdf_varid(fileID,'IVTV'),ivtV
   ncdf_varget,fileID,ncdf_varid(fileID,'RAINNC'),rain
   ncdf_varget,fileID,ncdf_varid(fileID,'Year'),y
   ncdf_varget,fileID,ncdf_varid(fileID,'Month'),m
   ncdf_varget,fileID,ncdf_varid(fileID,'Day'),d
   ncdf_varget,fileID,ncdf_varid(fileID,'Hour'),h
   ncdf_close,fileID
   
   ; Compute daily values from 3-hour values. Remove monthly climotology
   ; and divide by sigma for IVT
   init = 1
   for ik=min(d),max(d) do begin
      di  = where(d eq ik and m eq m(1))

      if (di(0) ne -1 and length(di) gt 1) then begin
         ; Compute ivt
         ivt = total(sqrt(ivtU(*,*,di)*ivtU(*,*,di)+ivtV(*,*,di)*ivtV(*,*,di)),3)/n_elements(di)

         ; Now remove climatology and standardize.
         ivtA = (ivt - ivtMclim(*,*,m(0)-1) ) / ivtSclim(*,*,m(0)-1)
   
         ; Store data
         if (init) then begin
            precip = rain(*,*,di(length(di)-1))-rain(*,*,di(0))
            if (nDay eq 0) then begin
               year   = y(di(0))
               month  = m(di(0))
               day    = d(di(0))
            endif
            if (nDay gt 0) then begin
               year   = [year,y(di(0))]
               month  = [month,m(di(0))]
               day    = [day,d(di(0))]
            endif
         endif
         if (not init) then begin
            precip = rain(*,*,di(length(di)-1))-rain(*,*,di(0)-1)
            year   = [year,y(di(0))]
            month  = [month,m(di(0))]
            day    = [day,d(di(0))]
         endif

         if (min(precip) lt 0) then begin
            goto,badPrecip
         endif

         testFail = 0
         test = where(finite(ivt) eq 0)
         if (test(0) ne -1) then testFail=1
      
         ; Compute precipitation anomaly.
         precipA = precip - precipMclim(*,*,m(0)-1)

         ; For this day, what BMU (SOM node) does it belong to?
         bmu  = bmusSOM(where(year(nDay) eq yearSOM and month(nDay) eq monthSOM and day(nDay) eq daySOM))
         bmui = bmu-1
      
         ; Add map to correct node.
         if (bmui(0) ge 0 and not testFail) then begin
            precipSOM(*,*,bmui)  = precipSOM(*,*,bmui)  + precip
            precipASOM(*,*,bmui) = precipASOM(*,*,bmui) + precipA
            ivtSOM(*,*,bmui)     = ivtSOM(*,*,bmui)     + ivt
            ivtASOM(*,*,bmui)    = ivtASOM(*,*,bmui)    + ivtA
            countSOM(bmui)      = countSOM(bmui)+1
         endif
         init = 0
         nDay = nDay + 1

         badPrecip:
      endif
   end
end

; Divide each SOM node by number of matching BMUs
for ij=0,num_rows*num_cols-1 do begin
   precipSOM(*,*,ij)  = precipSOM(*,*,ij)/countSOM(ij)
   precipASOM(*,*,ij) = precipASOM(*,*,ij)/countSOM(ij)
   ivtSOM(*,*,ij)     = ivtSOM(*,*,ij)/countSOM(ij)
   ivtASOM(*,*,ij)    = ivtASOM(*,*,ij)/countSOM(ij)
end
pat_freq = countSOM/total(countSOM)


; Write to netCDF file
fileOUT = 'SOMcomposites/comp.'+string(num_rows,format='(i1)')+$
          'x'+string(num_cols,format='(i1)')+'.'+SOMexpID+'.'+expID+'.nc'
fileID = ncdf_create(fileOUT,/clobber)
dimID1 = ncdf_dimdef(fileID,'lon',nlon)
dimID2 = ncdf_dimdef(fileID,'lat',nlat)
dimID3 = ncdf_dimdef(fileID,'BMU',num_rows*num_cols)
varID1 = ncdf_vardef(fileID,'lon',        [dimID1,dimID2       ],/float)
varID2 = ncdf_vardef(fileID,'lat',        [dimID1,dimID2       ],/float)
varID3 = ncdf_vardef(fileID,'BMU',        [              dimID3],/float)
varID4 = ncdf_vardef(fileID,'precip',     [dimID1,dimID2,dimID3],/float)
varID5 = ncdf_vardef(fileID,'IVT',        [dimID1,dimID2,dimID3],/float)
varID6 = ncdf_vardef(fileID,'precip_anom',[dimID1,dimID2,dimID3],/float)
varID7 = ncdf_vardef(fileID,'IVT_anom',   [dimID1,dimID2,dimID3],/float)
varID8 = ncdf_vardef(fileID,'pfreq',      [              dimID3],/float)
ncdf_attput,fileID,varID4,"comment","SOM composite: Precipitation (mm)"
ncdf_attput,fileID,varID5,"comment","SOM composite: IVT (kg/m/s)"
ncdf_attput,fileID,varID6,"comment","SOM composite: Precipitation Anomaly(mm)"
ncdf_attput,fileID,varID7,"comment","SOM composite: Standardized IVT Anomaly (kg/m/s)"

ncdf_attput,fileID,varID6,"comment","SOM Pattern Frequency (1)"
ncdf_control,fileID,/endef
ncdf_varput,fileID,varID1,lon
ncdf_varput,fileID,varID2,lat
ncdf_varput,fileID,varID3,indgen(num_rows*num_cols)+1
ncdf_varput,fileID,varID4,precipSOM
ncdf_varput,fileID,varID5,ivtSOM
ncdf_varput,fileID,varID6,precipASOM
ncdf_varput,fileID,varID7,ivtASOM
ncdf_varput,fileID,varID8,pat_freq
ncdf_close,fileID


; END PROGRAM
end
