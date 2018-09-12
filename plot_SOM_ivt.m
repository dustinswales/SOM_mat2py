%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of this program is to create composites for the 2x2 SOM runs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Housekeeping
clear all
close all

% Setup
expID      = 'WRF_gfdl_25km';
SOMexpID   = 'ONDJFM';
num_rows   = 3;
num_cols   = 3;
runID      = strcat(num2str(num_rows),'x',num2str(num_cols),'.',SOMexpID,'.',expID);

% WRF input data
dirIN      = 'SOMcomposites/';

% Output directory for plots
dirOUT     = '/home/dswales/Projects/NA-CORDEX/plots/SOMs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colortables
dirC  = '/home/dswales/tools/MATLAB/colortables/';
file1 = strcat(dirC,'ctable.70.CB-RdBu.txt');
file2 = strcat(dirC,'ctable.52.CB-GnBu.txt');
file3 = strcat(dirC,'ctable.58.CB-PuBuGn.txt');
file4 = strcat(dirC,'ctable.65.CB-YlOrBr.txt');
cmap  = textread(file1);
cmap2 = textread(file2);
cmap2(1:5,:)=1; %Make precip under 2.5mm white
cmap3 = textread(file3);
cmap4 = textread(file4);

% Composite plot domain
minLat = 20;
maxLat = 55;
minLon = 220;
maxLon = 260;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in WRF data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file        = strcat(dirIN,'comp.',runID,'.nc');
lon         = ncread(file,'lon');
lat         = ncread(file,'lat');
bmu         = ncread(file,'BMU');
ivt         = ncread(file,'IVT');
precip      = ncread(file,'precip');
precip_anom = ncread(file,'precip_anom');
ivt_anom    = ncread(file,'IVT_anom');
pat_freq    = ncread(file,'pfreq');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standardized IVT anomaly.
varName  = 'IVTanom';
varRange = [-2,2];
dVar     = 0.1;
varLName = 'Standardized IVT Anomaly';
plot_SOMcomposites(double(ivt_anom), double(lon), double(lat), pat_freq, varRange(1), varRange(2), dVar, minLon, minLat,...
            maxLon, maxLat, cmap, num_rows, num_cols, varLName)
saveas(1,strcat(dirOUT, 'SOM.', runID, varName, '.jpg'), 'jpg')
close all 

% IVT
varName  = 'IVT';
varRange = [10,200];
dVar     = 5.;
varLName = 'IVT';
plot_SOMcomposites(double(ivt), double(lon), double(lat), pat_freq, varRange(1), varRange(2), dVar, minLon, minLat,...
            maxLon, maxLat, cmap2, num_rows, num_cols, varLName)
saveas(1,strcat(dirOUT, 'SOM.', runID, varName, '.jpg'), 'jpg')
close all 

% Precipitation
varName  = 'Precip';
varRange = [0,20];
dVar     = 1;
varLName = 'Precipitation';
plot_SOMcomposites(double(precip), double(lon), double(lat), pat_freq, varRange(1), varRange(2), dVar, minLon, minLat,...
            maxLon, maxLat, cmap2, num_rows, num_cols, varLName)
saveas(1,strcat(dirOUT, 'SOM.', runID, varName, '.jpg'), 'jpg')
close all  




