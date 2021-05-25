# RTC_GreenInfrastructure
SWMM input files and Python files from a comparative simulation analysis of phosphorus removal of a real-time controlled bioretention cell to an uncontrolled bioretention cell (underdrain always open) with and without soil amendments under several influent concentrations and storm conditions.

# Model Parameters
* The hydrological model was built based on a bioretention cell in Toledo, Ohio was retrofitted for real-time control, including soil characteristics and infiltration rates. The parameters that could not be obtained from the Toledo site were based on design standards set forth by the Ohio Department of Natural Resources (Mathews 2006). 
* The water quality model parameters were experimentally derived and calibrated by Li and Davis (2016). 
* For more information, see SWMM input and Python files below.

# Phosphorus model implemented through StormReactor
* https://github.com/kLabUM/StormReactor

# SWMM input files and rain data file for all simulations
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

# Python files to run compartive simulation analysis, real weather simulation, and create heatmap of results
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
