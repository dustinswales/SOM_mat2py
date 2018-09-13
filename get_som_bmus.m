%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The purpose of this program is to compute the SOM for the given inputs
% and return the pattern-frequency and best-matching units (BMUS)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pat_freq,bmus] = get_som_bmus(num_rows,num_cols,var_anom,lattice,...
    shape,rad_ini,rad_fin,init,neighborhood_fct,trainlen_rough,...
    trainlen_finetune,num_obs)

    % Establish the data struct, as described in the SOM toolbox manual
    sD = som_data_struct(var_anom);
    
    % Initialize the map
    if init == 1
        sMap = som_lininit(sD, 'msize', [num_cols num_rows], lattice, shape);
    elseif init == 0
        sMap = som_randinit(sD, 'msize', [num_cols num_rows], lattice, shape);
    else
        disp('Improper SOM initialization assignment')
    end
    disp('Finished with SOM initialization')
    
    % Train the map with the batch version.
    sTrain            = som_train_struct('train', sD);
    sTrain.neigh      = neighborhood_fct;
    sTrain.radius_ini = rad_ini;
    sMap              = som_batchtrain(sMap, sD, sTrain, 'trainlen',...
        trainlen_rough);
    disp('Finished with first phase of SOM training')
    
    % Finetune in the second phase of training
    sTrain.radius_ini = rad_fin + 2;
    sTrain.radius_fin = rad_fin;
    [sMap, sTopol]    = som_batchtrain(sMap,sD, sTrain, 'trainlen',...
        trainlen_finetune);
    disp('Finished with second phase of SOM training')
    
    % Get the average quantization (rms) error and topographic error
    % as measures of quality
    [aqe, te] = som_quality(sMap,sD);
    
    % Get the best-matching unit (BMU) time series and quantization
    % error for each map.
    [bmus, qerrs] = som_bmus(sMap, sD);
    
    % Determine the frequency of occurrence of each pattern
    pat_freq = zeros(1,num_rows*num_cols);
    for p = 1:num_rows*num_cols
        ind = find(bmus == p);
        pat_freq(p) = length(ind)/num_obs;
    end
    
    % Get the SOM pattern matrix from the codebook
    pats = sMap.codebook;

end