This directory includes the code used to post-process cyclone track data and ERA5/CMIP5 rainfall data for Pepler & Dowdy (2021)

Cyclone tracks were generated using the Melbourne Univeristy cyclone tracking scheme (https://cyclonetracker.earthsci.unimelb.edu.au/) 
and are available for research use at: XXXX

Final post-processed netCDF files used in the generation of figures and data for the paper are available on figshare at YYYY:

The code here uses a mix of R and NCL scripts, with some minor CDO used for additional postprocessing (e.g. regridding). 
It was set up based on the directory structures of data stored at NCI, so will likely need additional editing for your own system

ANALYSIS STEPS:

1. Extract/download cyclone tracks (see above). 

These files have one row for every cyclone centre, giving the date, time, location, intensity, etc.

2. Convert cyclone data to daily grids of cyclone location using "2_point2grid_cyclonedepth.R"

This code requires both surface and upper cyclone data files for a given model as well as vectors of latitude and longitude, 
and outputs a netCDF which has one grid per day showing areas influenced by no cyclone (0), a deep cyclone (1), a shallow surface
cyclone (2) or a shallow upper cyclone (3)

You can choose a number of parameters, including:
- The minimum intensity thresholds for each level
- Other driteria such as requiring lows to be closed, to persist a minimum number of time steps or minimum track movement (in km)
- The radius (in degrees) used for a "cyclone region"
- The maximum distance between upper and surface cyclones required for a cyclone to be "deep", 
and whether this should be applied to each low individually or to the event as a whole

3. Calculate the relevant percentile threshold for each grid point from model data using "3_calculate_percentile.ncl"

Calculates percentiles from all daily data for a specified threshold and time period

4. Calculate the total rainfall and days above percentile thresholds on a monthly basis
for the three cyclone categories (and other days) using "4_gridrain_cyclonedepth.ncl"

4a. Use CDO to combine the annual files into a single file for the whole period and regrid to a common resolution,
producing the final netCDF files used for analysis (and vailable at YYYY)

cdo enssum $FILELIST tmp.nc
cdo remapbil,mygrid tmp.nc $FOUT

5. Analyse rainfall results and produce figures using "5_analyse_CMIP5_rainfall_data.R"
This includes figures 1-4, supplementary figures S2-4 and S6 and the mean changes in latitutindal ranges reported in the paper 
[and associated .RData files at ZZZZ?]

5a. Additional analysis of raw surface and upper cyclone tracks using "5a_analyse_CMIP5_cyclone_data.R"
This produces supplementary figures S1 and S5

