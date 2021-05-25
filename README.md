# RTC_GreenInfrastructure
SWMM input files and Python files from a comparative simulation analysis of phosphorus removal of a real-time controlled bioretention cell to an uncontrolled bioretention cell (underdrain always open) with and without soil amendments under several influent concentrations and storm conditions.

## Water Quality Model
* Phosphorus water quality model implemented using [StormReactor](https://github.com/kLabUM/StormReactor), a Python package for modeling any pollutant generation or removal method in SWMM.
* Model parameters were experimentally derived and calibrated by Li and Davis (2016) and are provided in the Python files below. 

## Hydrological Model
* The hydrological model was built based on a bioretention cell in Toledo, Ohio was retrofitted for real-time control, including soil characteristics and infiltration rates. 
* The parameters that could not be obtained from the Toledo site were based on design standards set forth by the Ohio Department of Natural Resources (Mathews 2006). See SWMM files below.

### Python files
* TP_Sims_6hr0.5in_0.2mgL.py
* TP_Sims_6hr0.5in_0.6mgL.py
* TP_Sims_6hr0.5in_1.0mgL.py
* TP_Sims_6hr0.5in_1.4mgL.py
* TP_Sims_6hr0.5in_1.8mgL.py
* TP_Sims_6hr1in_0.2mgL.py
* TP_Sims_6hr1in_0.6mgL.py
* TP_Sims_6hr1in_1.0mgL.py
* TP_Sims_6hr1in_1.4mgL.py
* TP_Sims_6hr1in_1.8mgL.py
* TP_Sims_6hr2in_0.2mgL.py
* TP_Sims_6hr2in_0.6mgL.py
* TP_Sims_6hr2in_1.0mgL.py
* TP_Sims_6hr2in_1.4mgL.py
* TP_Sims_6hr2in_1.8mgL.py
* TP_Sims_RealWeather_0.38mgL.py
* heatmap.py

### SWMM files
* BRC_6hr0.5in_0.2mgL.inp
* BRC_6hr0.5in_0.6mgL.inp
* BRC_6hr0.5in_1.0mgL.inp
* BRC_6hr0.5in_1.4mgL.inp
* BRC_6hr0.5in_1.8mgL.inp
* BRC_6hr1in_0.2mgL.inp
* BRC_6hr1in_0.38mgL.inp
* BRC_6hr1in_0.6mgL.inp
* BRC_6hr1in_1.0mgL.inp
* BRC_6hr1in_1.4mgL.inp
* BRC_6hr1in_1.8mgL.inp
* BRC_6hr2in_0.2mgL.inp
* BRC_6hr2in_0.6mgL.inp
* BRC_6hr2in_1.0mgL.inp
* BRC_6hr2in_1.4mgL.inp
* BRC_6hr2in_1.8mgL.inp
* BRC_RealWeather_0.38mgL.inp
* Rain_Law_6-11-2014.dat
