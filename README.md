This repository contains matlab code to compute the Self Organizing Maps (SOMs) of the
integrated vapor transport (IVT) over the intermountain west (IMW).

The code has been recently modified to use 3-hour WRF-RCM data from NA-CORDEX.

This is a two-step process:
1) Generate the BMUs (see tranSOM_vt.m). This will make text files, and .mat files, containing
   the BMU for each day. Forewarning, there's quite a bit of pre-processing going on in here to
   compute daily values from the 3-hour data.
2) a) Create synoptic composites from the BMUs. I do this in IDL (composite.pro), which...
      i)   Read in the BMUs from step 1.
      ii)  Read in WRF climatology. For simplicity I pre-computed the monthly climotologies
           and ingest them here.
      iii) Read in WRF data (only IVT and precipitation here), compute daily values, assign daily
           map to correct SOM node.
      iv)  Output to netCDF.
    b) Plot SOMs (see plot_SOM_ivt.m and [plot_SOMcomposites.m) for examples.