%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of this program is to create SOMs of the standardized IVT
% anomalies. Only land-points during the cool-season (ONDJFM) are used.
%
%
% The dates of the BMUs are saved in .mat files located in bmus/runID.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Housekeeping
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input data control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Which GCMs to use?
expIDs = {'WRF_hadgem_50km','WRF_mpi_50km','WRF_erain_50km','WRF_mpi_25km',...
    'WRF_erain_25km','WRF_hadgem_25km','WRF_gfdl_50km','WRF_gfdl_25km'};
lsMaskFile = {'wrf.xland.50km.nc','wrf.xland.50km.nc','wrf.xland.50km.nc',...
    'wrf.xland.25km.nc','wrf.xland.25km.nc','wrf.xland.25km.nc',...
    'wrf.xland.50km.nc','wrf.xland.25km.nc'};

% Where to store SOM output (BMUS)?
dirOUT     = '/data/dswales/NA-CORDEX/SOMs/';

% Experiment ID?
SOMexpID   = 'ONDJFM';

% Establish the latitude and longitude max and mins for the SOM analysis.
latmax = 50;
latmin = 30;
lonmax = -130;
lonmin = -110;

% SOM parameters
num_rowsi         = [2,3,4]; % Number of rows in SOM
num_colsi         = [2,3,4]; % Number of columns in SOM
init              = 1;       % 1 = linear initialization, 0 = random
lattice           = 'rect';  % 'rect' for rectangular or 'hexa' for hexagonal
shape             = 'sheet'; % 'sheet', 'cyl', or 'toroid'
neighborhood_fct  = 'ep';    % 'gaussian', 'cutgauss', 'bubble', or 'ep'
rad_ini           = 3;       % Initial neighborhood radius used in the training
rad_fin           = 0;       % Final neighborhood radius at the end of SOM
trainlen_rough    = 50;      % Iterations during the rough initial phase
trainlen_finetune = 50;      % Iterations during the finetuning phase

for iMod = 1:length(expIDs)
    expID      = string(expIDs(iMod));
    dirIN      = strcat('/Projects/HydroMet/dswales/NA-CORDEX/',expID,'/');
    fileLSMASK = strcat('/home/dswales/Projects/NA-CORDEX/data/',string(lsMaskFile(iMod)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read in WRF IVT data.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Read in files for cool-season months (ONDJFM).
    files01 = dir(fullfile(dirIN, 'wrfout_*01.nc'));
    files02 = dir(fullfile(dirIN, 'wrfout_*02.nc'));
    files03 = dir(fullfile(dirIN, 'wrfout_*03.nc'));
    files10 = dir(fullfile(dirIN, 'wrfout_*10.nc'));
    files11 = dir(fullfile(dirIN, 'wrfout_*11.nc'));
    files12 = dir(fullfile(dirIN, 'wrfout_*12.nc'));
    files   = [files01;files02;files03;files10;files11;files12];
    
    % Read in longitude and latitude
    lon  = ncread(strcat(dirIN,files(1).name),'XLONG');
    lat  = ncread(strcat(dirIN,files(1).name),'XLAT');
    xi   = find(lon(:,1) >= lonmax & lon(:,1) <= lonmin);
    yi   = find(lat(1,:) >= latmin & lat(1,:) <= latmax);
    lon  = lon(min(xi):max(xi),min(yi):max(yi));
    lat  = lat(min(xi):max(xi),min(yi):max(yi));
    nlon = size(lon,1);
    nlat = size(lat,2);
    
    % Read in IVT
    init = 0;
    for ij=1:size(files,1)
        ij
        y      = ncread(strcat(dirIN,files(ij).name),'Year');
        m      = ncread(strcat(dirIN,files(ij).name),'Month');
        d      = ncread(strcat(dirIN,files(ij).name),'Day');
        h      = ncread(strcat(dirIN,files(ij).name),'Hour');
        nTime  = length(y);
        ivtU   = ncread(strcat(dirIN,files(ij).name),'IVTU',[min(xi),min(yi),1],[nlon,nlat,nTime],[1,1,1]);
        ivtV   = ncread(strcat(dirIN,files(ij).name),'IVTV',[min(xi),min(yi),1],[nlon,nlat,nTime],[1,1,1]);
        
        % Compute daily values
        for ik=min(d):max(d)
            di=find(d == ik);
            % Store data
            if (isempty(di) == 0)
                if (init == 0)
                    ivtTot = mean(sqrt(ivtU(:,:,di).*ivtU(:,:,di)+ivtV(:,:,di).*ivtV(:,:,di)),3);
                    year  = y(di(1));
                    month = m(di(1));
                    day   = d(di(1));
                    init  = 1;
                else
                    ivtTot = cat(3,ivtTot,mean(sqrt(ivtU(:,:,di).*ivtU(:,:,di)+ivtV(:,:,di).*ivtV(:,:,di)),3));
                    year  = [year;y(di(1))];
                    month = [month;m(di(1))];
                    day   = [day;d(di(1))];
                end
            end
        end
    end
    nTime = size(year,1);
    
    % Standardize: Remove climatology and divide by standard deviation
    fileClim = strcat('/home/dswales/Projects/NA-CORDEX/data/clim/',expID,'/ivt.mon.clim.nc');
    ivtClim  = ncread(fileClim,'mean',[min(xi),min(yi),1],[nlon,nlat,12],[1,1,1]);
    ivtSTD   = ncread(fileClim,'stdev',[min(xi),min(yi),1],[nlon,nlat,12],[1,1,1]);
    ivtTotN = zeros(nlon,nlat,nTime);
    for ij=1:nTime
        ivtTotN(:,:,ij) = (ivtTot(:,:,ij)-ivtClim(:,:,month(ij)))./ivtSTD(:,:,month(ij));
    end
    
    % Reshape array (time,xyPts).
    ivtTotNb = reshape(ivtTotN,nlon*nlat,nTime)';
    
    % Guard against bad IVT data
    [ti,xii] = find(isinf(ivtTotNb)==1);
    if (isempty(ti)==0)
        badTimeIndices = unique(ti);
        nBadTimes = length(badTimeIndices);
        disp('Not Including the following days in SOM analysis:')
        for ij=1:nBadTimes
            disp(strcat(num2str(month(badTimeIndices(ij))),'/',num2str(day(badTimeIndices(ij))),...
                '/',num2str(year(badTimeIndices(ij)))))
        end
        tii = boolean([1:length(year)]);
        tii(badTimeIndices)=0;
        ivtTotNb = ivtTotNb(tii,:);
        year     = year(tii);
        month    = month(tii);
        day      = day(tii);
    end
    
    % Mask out oceanic scences.
    lsmask     = ncread(fileLSMASK,'XLAND',[min(xi),min(yi),1],[nlon,nlat,1],[1,1,1]);
    lsmask     = reshape(lsmask,nlon*nlat,1);
    land       = find(lsmask == 1);
    ivtTotNb   = ivtTotNb(:,land);
    
    % Guard against any NaN's
    badPts = find(isinf(ivtTotNb)==1);
    if (isempty(badPts)==0)
        ivtTotNb(badPts) = 0.0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOM Analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for isom = 1:length(num_rowsi)
        close all
        num_rows = num_rowsi(isom);
        num_cols = num_colsi(isom);
        runID = strcat(num2str(num_rows),'x',num2str(num_cols),'.',SOMexpID);
        disp(strcat('Creating SOM for ',num2str(num_rows),' rows and ',num2str(num_cols),' columns'))
        
        % Call SOM subroutine and get pattern frequencies and BMUs
        [pat_freq,bmus] = get_som_bmus(num_rows,num_cols,ivtTotNb,lattice,...
            shape,rad_ini,rad_fin,init,...
            neighborhood_fct,trainlen_rough,...
            trainlen_finetune,nTime);
        
        % Save BMU dates
        fileOUT = strcat(dirOUT,'bmus/BMU.',runID,'.',expID,'.mat');
        save(fileOUT,'year','month','day','pat_freq','bmus','ivtTotNb')
        bunk   = find(isnan(bmus(:)) == 1);
        nobunk = find(isnan(bmus(:)) == 0);
        fileOUT = strcat(dirOUT,'bmus/BMU.',runID,'.',expID,'.txt');
        fileID=fopen(fileOUT,'w');
        fprintf(fileID,'%s\n','year');
        fprintf(fileID,'%i\n',size(year(nobunk),1));
        fprintf(fileID,'%6i\n',year(nobunk));
        fprintf(fileID,'%s\n','month');
        fprintf(fileID,'%i\n',size(month(nobunk),1));
        fprintf(fileID,'%6i\n',month(nobunk));
        fprintf(fileID,'%s\n','day');
        fprintf(fileID,'%i\n',size(day(nobunk),1));
        fprintf(fileID,'%6i\n',day(nobunk));
        fprintf(fileID,'%s\n','bmus');
        fprintf(fileID,'%i\n',size(bmus(nobunk),1));
        fprintf(fileID,'%6i\n',bmus(nobunk));
        fclose(fileID);
    end
    
end

