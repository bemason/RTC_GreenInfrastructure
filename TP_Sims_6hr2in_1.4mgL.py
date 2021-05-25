# Import modules
from pyswmm import Simulation, Nodes, Links
import numpy as np
import statistics as st
import matplotlib.pyplot as plt

# Traditional soil parameters
B1_pp = 0.000217
Ceq0_pp = 0.376
k_pp = 0.000533
B1_dp = 0.000167
Ceq0_dp = 0.368
k_dp = 0.000867
L = 3.0
A = 4357
E = 0.44
t = 0

# Lists to store results
inflow = []         # Total inflow to the BRC (ft^3/s)
runoff = []         # Surface runoff from the BRC (ft^3/s)
pond_height = []    # Ponding height in BRC (ft)
soil_height = []    # Height of water in BRC soil (ft)
drain_flow = []     # Drain outflow from the BRC (ft^3/s)
infil_rate = []     # Infiltration rate through the BRC (ft^3/s)
exfil_rate = []     # Exfiltration rate from the BRC (ft^3/s)

inflow_dp = []      # Influent DP concentration (mg/L)
inflow_pp = []      # Influent PP concentration (mg/L)
runoff_dp = []      # Runoff DP concentration (mg/L)
runoff_pp = []      # Runoff PP concentration (mg/L)
infil_dp = []       # Infiltrated DP concentration (mg/L)
infil_pp = []       # Infiltrated PP concentration (mg/L)
drained_dp = []     # Drained DP concentration (mg/L)
drained_pp = []     # Drained PP concentration (mg/L)
exfil_dp = []       # Exfiltrated DP concentration (mg/L)
exfil_pp = []       # Exfiltrated PP concentration (mg/L)


# Setup toolbox simulation for traditional soils
with Simulation("./BRC_6hr2in_1.4mgL.inp") as sim:
    
    # Get asset information
    soil = Links(sim)["IntoMedia"]
    pond = Nodes(sim)["Pond"]
    drained = Nodes(sim)["Drained"]
    exfiltrated = Nodes(sim)["Infiltrated"]
    overflow = Nodes(sim)["Runoff"]

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Get phosphorus concentrations        
        Cin_dp = soil.pollut_quality['DP']
        inflow_dp.append(Cin_dp)
        Cin_pp = soil.pollut_quality['PP']
        inflow_pp.append(Cin_pp)

        # Get depth data
        pond_height.append(pond.depth)
        soil_height.append(soil.depth)

        # Getflow data
        inflow.append(pond.total_inflow)
        runoff.append(overflow.total_inflow)
        Qin = soil.flow
        infil_rate.append(Qin)
        drain_flow.append(drained.total_inflow)
        exfil_rate.append(exfiltrated.total_inflow)

        if index < 10000:
            # Phosphorus Model
            if Qin >= 0.0001:
                # Accumulate time elapsed since water entered node
                t = t + 1
                # Calculate new concentration
                Cnew_dp = (Cin_dp*np.exp((-k_dp*L*A*E)/Qin))+\
                    (Ceq0_dp*np.exp(B1_dp*t))*(1-(np.exp((-k_dp*L*A*E)/Qin)))
                Cnew_pp = (Cin_pp*np.exp((-k_pp*L*A*E)/Qin))+\
                    (Ceq0_pp*np.exp(B1_pp*t))*(1-(np.exp((-k_pp*L*A*E)/Qin)))

                # Set new concentration
                sim._model.setNodePollutant("Media", 0, Cnew_pp)
                sim._model.setNodePollutant("Media", 1, Cnew_dp)
                #print("Cnew DP & PP:", Cnew_dp, Cnew_pp)
            else:
                Cnew_pp = Ceq0_pp
                Cnew_dp = Ceq0_dp  
                sim._model.setNodePollutant("Media", 0, Cnew_pp)
                sim._model.setNodePollutant("Media", 1, Cnew_dp)
                #print("Cnew2 DP & PP:", Cnew_dp, Cnew_pp)
        else:
            Cnew_pp = Ceq0_pp
            Cnew_dp = Ceq0_dp  
            sim._model.setNodePollutant("Media", 0, Cnew_pp)
            sim._model.setNodePollutant("Media", 1, Cnew_dp)

        # Save phosphorus concentrations for the rest of the system 
        runoff_dp.append(overflow.pollut_quality['DP'])
        runoff_pp.append(overflow.pollut_quality['PP'])    
        infil_dp.append(Cnew_dp)
        infil_pp.append(Cnew_pp)
        drained_dp.append(Cnew_dp)
        drained_pp.append(Cnew_pp)
        exfil_dp.append(Cnew_dp)
        exfil_pp.append(Cnew_pp)

    sim._model.swmm_end()
    #print("Runoff Error:", sim.runoff_error)
    #print("Routing Error:", sim.flow_routing_error)
    #print("Quality Error:", sim.quality_error)

# Convert flow rates from CFS to LPS
conv_cfs_lps = [28.3168]*len(inflow)
inflow_m = [a*b for a,b in zip(inflow,conv_cfs_lps)]
runoff_m = [a*b for a,b in zip(runoff,conv_cfs_lps)]
infil_m = [a*b for a,b in zip(infil_rate,conv_cfs_lps)]
drain_m = [a*b for a,b in zip(drain_flow,conv_cfs_lps)]
exfil_m = [a*b for a,b in zip(exfil_rate,conv_cfs_lps)]

# Convert height from ft to m
conv_ft_m = [0.3048]*len(inflow)
pond_height_m = [a*b for a,b in zip(pond_height,conv_ft_m)]
soil_height_m = [a*b for a,b in zip(soil_height,conv_ft_m)]

# Calculate TTP concentration
inflow_tp = [a+b for a,b in zip(inflow_dp, inflow_pp)]
runoff_tp = [a+b for a,b in zip(runoff_dp, runoff_pp)]
infil_tp = [a+b for a,b in zip(infil_dp, infil_pp)]
drain_tp = [a+b for a,b in zip(drained_dp, drained_pp)]
exfil_tp = [a+b for a,b in zip(exfil_dp, exfil_pp)]

# Calculate DP load (g) each timestep
conv_mgs_gs = [0.001]*len(inflow)
timestep = [5]*len(inflow)
inflow_load_dp = [a*b*c*d for a,b,c,d in zip(inflow_dp,inflow_m,conv_mgs_gs,timestep)]
runoff_load_dp = [a*b*c*d for a,b,c,d in zip(runoff_dp,runoff_m,conv_mgs_gs,timestep)]
infil_load_dp = [a*b*c*d for a,b,c,d in zip(infil_dp,infil_m,conv_mgs_gs,timestep)]
drain_load_dp = [a*b*c*d for a,b,c,d in zip(drained_dp,drain_m,conv_mgs_gs,timestep)]
exfil_load_dp = [a*b*c*d for a,b,c,d in zip(exfil_dp,exfil_m,conv_mgs_gs,timestep)]

# Calculate DP load (g) each timestep
inflow_load_pp = [a*b*c*d for a,b,c,d in zip(inflow_pp,inflow_m,conv_mgs_gs,timestep)]
runoff_load_pp = [a*b*c*d for a,b,c,d in zip(runoff_pp,runoff_m,conv_mgs_gs,timestep)]
infil_load_pp = [a*b*c*d for a,b,c,d in zip(infil_dp,infil_m,conv_mgs_gs,timestep)]
drain_load_pp = [a*b*c*d for a,b,c,d in zip(drained_pp,drain_m,conv_mgs_gs,timestep)]
exfil_load_pp = [a*b*c*d for a,b,c,d in zip(exfil_pp,exfil_m,conv_mgs_gs,timestep)]

# Calculate cumulative DP load (g)
inflow_cumload_dp = np.cumsum(inflow_load_dp)
runoff_cumload_dp = np.cumsum(runoff_load_dp)
infil_cumload_dp = np.cumsum(infil_load_dp)
drain_cumload_dp = np.cumsum(drain_load_dp)
exfil_cumload_dp = np.cumsum(exfil_load_dp)

# Calculate cumulative PP load (g)
inflow_cumload_pp = np.cumsum(inflow_load_pp)
runoff_cumload_pp = np.cumsum(runoff_load_pp)
infil_cumload_pp = np.cumsum(infil_load_pp)
drain_cumload_pp = np.cumsum(drain_load_pp)
exfil_cumload_pp = np.cumsum(exfil_load_pp)

# Calculate TP load (g)
influent_load = [a+b for a,b in zip(inflow_load_dp, inflow_load_pp)]
runoff_load = [a+b for a,b in zip(runoff_load_dp, runoff_load_pp)]
infil_load = [a+b for a,b in zip(infil_load_dp, infil_load_pp)]
drain_load = [a+b for a,b in zip(drain_load_dp, drain_load_pp)]
released_load = [a+b+c+d for a,b,c,d in zip(runoff_load_dp, runoff_load_pp, drain_load_dp, drain_load_pp)]
exfil_load = [a+b for a,b in zip(exfil_load_dp, exfil_load_pp)]
captured_load = [a-b for a,b in zip(influent_load, released_load)]

# Calculate cumulative TP load (g)
influent_cumload = np.cumsum(influent_load)
released_cumload = np.cumsum(released_load)
captured_cumload= np.cumsum(captured_load)

#_______________________________________________________________________

# Traditional soil parameters
t = 0

# Lists to store results
inflowC = []         # Total inflow to the BRC (ft^3/s)
runoffC = []         # Surface runoff from the BRC (ft^3/s)
pond_heightC = []    # Ponding height in BRC (ft)
soil_heightC = []    # Height of water in BRC soil (ft)
drain_flowC = []     # Drain outflow from the BRC (ft^3/s)
infil_rateC = []     # Infiltration rate through the BRC (ft^3/s)
exfil_rateC = []     # Exfiltration rate from the BRC (ft^3/s)
valve_setting = []  # Underdrain valve settings (%)

inflow_dpC = []      # Influent DP concentration (mg/L)
inflow_ppC = []      # Influent PP concentration (mg/L)
runoff_dpC = []      # Runoff DP concentration (mg/L)
runoff_ppC = []      # Runoff PP concentration (mg/L)
infil_dpC = []       # Infiltrated DP concentration (mg/L)
infil_ppC = []       # Infiltrated PP concentration (mg/L)
drained_dpC = []     # Drained DP concentration (mg/L)
drained_ppC = []     # Drained PP concentration (mg/L)
exfil_dpC = []       # Exfiltrated DP concentration (mg/L)
exfil_ppC = []       # Exfiltrated PP concentration (mg/L)

# Setup toolbox simulation for RTC of BRC with traditional soils
with Simulation("./BRC_6hr2in_1.4mgL.inp") as sim:
    
    # Get asset information
    soil = Links(sim)["IntoMedia"]
    pond = Nodes(sim)["Pond"]
    drain = Links(sim)["Underdrain"]
    drained = Nodes(sim)["Drained"]
    exfiltrated = Nodes(sim)["Infiltrated"]
    overflow = Nodes(sim)["Runoff"]

    # Tracking time for control actions every 15 minutes (5 sec time step)
    _tempcount = 180

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Get phosphorus concentrations        
        Cin_dpC = soil.pollut_quality['DP']
        inflow_dpC.append(Cin_dpC)
        Cin_ppC = soil.pollut_quality['PP']
        inflow_ppC.append(Cin_ppC)

        # Get depth data
        h = pond.depth
        pond_heightC.append(h)
        soil_heightC.append(soil.depth)

        # Getflow data
        inflowC.append(pond.total_inflow)
        runoffC.append(overflow.total_inflow)
        Qin = soil.flow
        infil_rateC.append(Qin)
        drain_flowC.append(drained.total_inflow)
        exfil_rateC.append(exfiltrated.total_inflow)

        if index < 10000:
            # Phosphorus Model
            if Qin >= 0.0001:
                # Accumulate time elapsed since water entered node
                t = t + 1
                # Calculate new concentration
                Cnew_dp = (Cin_dp*np.exp((-k_dp*L*A*E)/Qin))+\
                    (Ceq0_dp*np.exp(B1_dp*t))*(1-(np.exp((-k_dp*L*A*E)/Qin)))
                Cnew_pp = (Cin_pp*np.exp((-k_pp*L*A*E)/Qin))+\
                    (Ceq0_pp*np.exp(B1_pp*t))*(1-(np.exp((-k_pp*L*A*E)/Qin)))

                # Set new concentration
                sim._model.setNodePollutant("Media", 0, Cnew_pp)
                sim._model.setNodePollutant("Media", 1, Cnew_dp)
                #print("Cnew DP & PP:", Cnew_dp, Cnew_pp)
            else:
                Cnew_pp = Ceq0_pp
                Cnew_dp = Ceq0_dp  
                sim._model.setNodePollutant("Media", 0, Cnew_pp)
                sim._model.setNodePollutant("Media", 1, Cnew_dp)
                #print("Cnew2 DP & PP:", Cnew_dp, Cnew_pp)
        else:
            Cnew_pp = Ceq0_pp
            Cnew_dp = Ceq0_dp  
            sim._model.setNodePollutant("Media", 0, Cnew_pp)
            sim._model.setNodePollutant("Media", 1, Cnew_dp)

        # Save phosphorus concentrations for the rest of the system 
        runoff_dpC.append(overflow.pollut_quality['DP'])
        runoff_ppC.append(overflow.pollut_quality['PP'])    
        infil_dpC.append(Cnew_dp)
        infil_ppC.append(Cnew_pp)
        drained_dpC.append(Cnew_dp)
        drained_ppC.append(Cnew_pp)
        exfil_dpC.append(Cnew_dp)
        exfil_ppC.append(Cnew_pp)

        
        # Real-time control
        if _tempcount == 180:
            drain.target_setting = 0.1*h
            _tempcount= 0
        _tempcount+= 1
        
        valve_setting.append(drain.target_setting)
        
    sim._model.swmm_end()
    #print("Runoff Error:", sim.runoff_error)
    #print("Routing Error:", sim.flow_routing_error)
    #print("Quality Error:", sim.quality_error)

# Convert flow rates from CFS to LPS
conv_cfs_lps = [28.3168]*len(inflowC)
inflow_mC = [a*b for a,b in zip(inflowC,conv_cfs_lps)]
runoff_mC = [a*b for a,b in zip(runoffC,conv_cfs_lps)]
infil_mC = [a*b for a,b in zip(infil_rateC,conv_cfs_lps)]
drain_mC = [a*b for a,b in zip(drain_flowC,conv_cfs_lps)]
exfil_mC = [a*b for a,b in zip(exfil_rateC,conv_cfs_lps)]

# Convert height from ft to m
conv_ft_m = [0.3048]*len(inflow)
pond_height_mC = [a*b for a,b in zip(pond_heightC,conv_ft_m)]
soil_height_mC = [a*b for a,b in zip(soil_heightC,conv_ft_m)]

# Calculate TP concentration
inflow_tpC = [a+b for a,b in zip(inflow_dpC, inflow_ppC)]
runoff_tpC = [a+b for a,b in zip(runoff_dpC, runoff_ppC)]
infil_tpC = [a+b for a,b in zip(infil_dpC, infil_ppC)]
drain_tpC = [a+b for a,b in zip(drained_dpC, drained_ppC)]
exfil_tpC = [a+b for a,b in zip(exfil_dpC, exfil_ppC)]

# Calculate DP load (g) each timestep
conv_mgs_gs = [0.001]*len(inflow)
timestep = [5]*len(inflow)
inflow_load_dpC = [a*b*c*d for a,b,c,d in zip(inflow_dpC,inflow_mC,conv_mgs_gs,timestep)]
runoff_load_dpC = [a*b*c*d for a,b,c,d in zip(runoff_dpC,runoff_mC,conv_mgs_gs,timestep)]
infil_load_dpC = [a*b*c*d for a,b,c,d in zip(infil_dpC,infil_mC,conv_mgs_gs,timestep)]
drain_load_dpC = [a*b*c*d for a,b,c,d in zip(drained_dpC,drain_mC,conv_mgs_gs,timestep)]
exfil_load_dpC = [a*b*c*d for a,b,c,d in zip(exfil_dpC,exfil_mC,conv_mgs_gs,timestep)]

# Calculate DP load (g) each timestep
inflow_load_ppC = [a*b*c*d for a,b,c,d in zip(inflow_ppC,inflow_mC,conv_mgs_gs,timestep)]
runoff_load_ppC = [a*b*c*d for a,b,c,d in zip(runoff_ppC,runoff_mC,conv_mgs_gs,timestep)]
infil_load_ppC = [a*b*c*d for a,b,c,d in zip(infil_dpC,infil_mC,conv_mgs_gs,timestep)]
drain_load_ppC = [a*b*c*d for a,b,c,d in zip(drained_ppC,drain_mC,conv_mgs_gs,timestep)]
exfil_load_ppC = [a*b*c*d for a,b,c,d in zip(exfil_ppC,exfil_mC,conv_mgs_gs,timestep)]

# Calculate cumulative DP load (g)
inflow_cumload_dpC = np.cumsum(inflow_load_dpC)
runoff_cumload_dpC = np.cumsum(runoff_load_dpC)
infil_cumload_dpC = np.cumsum(infil_load_dpC)
drain_cumload_dpC = np.cumsum(drain_load_dpC)
exfil_cumload_dpC = np.cumsum(exfil_load_dpC)

# Calculate cumulative PP load (g)
inflow_cumload_ppC = np.cumsum(inflow_load_ppC)
runoff_cumload_ppC = np.cumsum(runoff_load_ppC)
infil_cumload_ppC = np.cumsum(infil_load_ppC)
drain_cumload_ppC = np.cumsum(drain_load_ppC)
exfil_cumload_ppC = np.cumsum(exfil_load_ppC)

# Calculate TP load (g)
influent_loadC = [a+b for a,b in zip(inflow_load_dpC, inflow_load_ppC)]
runoff_loadC = [a+b for a,b in zip(runoff_load_dpC, runoff_load_ppC)]
infil_loadC = [a+b for a,b in zip(infil_load_dpC, infil_load_ppC)]
drain_loadC = [a+b for a,b in zip(drain_load_dpC, drain_load_ppC)]
released_loadC = [a+b+c+d for a,b,c,d in zip(runoff_load_dpC, runoff_load_ppC, drain_load_dpC, drain_load_ppC)]
exfil_loadC = [a+b for a,b in zip(exfil_load_dpC, exfil_load_ppC)]
captured_loadC = [a-b for a,b in zip(influent_loadC, released_loadC)]

# Calculate cumulative TP load (g)
influent_cumloadC = np.cumsum(influent_loadC)
released_cumloadC = np.cumsum(released_loadC)
captured_cumloadC = np.cumsum(captured_loadC)

#_______________________________________________________________________

# Amended Soil parameters
B1_ppA = 0.0000833
Ceq0_ppA = 0.0226
k_ppA = 0.00143
B1_dpA = 0.0000333
Ceq0_dpA = 0.0081
k_dpA = 00.00320
t = 0

# Lists to store results
inflowA = []         # Total inflow to the BRC (ft^3/s)
runoffA = []         # Surface runoff from the BRC (ft^3/s)
pond_heightA = []    # Ponding height in BRC (ft)
soil_heightA = []    # Height of water in BRC soil (ft)
drain_flowA = []     # Drain outflow from the BRC (ft^3/s)
infil_rateA = []     # Infiltration rate through the BRC (ft^3/s)
exfil_rateA = []     # Exfiltration rate from the BRC (ft^3/s)

inflow_dpA = []      # Influent DP concentration (mg/L)
inflow_ppA = []      # Influent PP concentration (mg/L)
runoff_dpA = []      # Runoff DP concentration (mg/L)
runoff_ppA = []      # Runoff PP concentration (mg/L)
infil_dpA = []       # Infiltrated DP concentration (mg/L)
infil_ppA = []       # Infiltrated PP concentration (mg/L)
drained_dpA = []     # Drained DP concentration (mg/L)
drained_ppA = []     # Drained PP concentration (mg/L)
exfil_dpA = []       # Exfiltrated DP concentration (mg/L)
exfil_ppA = []       # Exfiltrated PP concentration (mg/L)

# Setup toolbox simulation for amended soils
with Simulation("./BRC_6hr2in_1.4mgL.inp") as sim:
    
    # Get asset information
    soil = Links(sim)["IntoMedia"]
    pond = Nodes(sim)["Pond"]
    drained = Nodes(sim)["Drained"]
    exfiltrated = Nodes(sim)["Infiltrated"]
    overflow = Nodes(sim)["Runoff"]

    # Step through the simulation    
    for index,step in enumerate(sim):

        # Get phosphorus concentrations        
        Cin_dpA = soil.pollut_quality['DP']
        inflow_dpA.append(Cin_dpA)
        Cin_ppA = soil.pollut_quality['PP']
        inflow_ppA.append(Cin_ppA)

        # Get depth data
        pond_heightC.append(pond.depth)
        soil_heightC.append(soil.depth)

        # Getflow data
        inflowA.append(pond.total_inflow)
        runoffA.append(overflow.total_inflow)
        Qin = soil.flow
        infil_rateA.append(Qin)
        drain_flowA.append(drained.total_inflow)
        exfil_rateA.append(exfiltrated.total_inflow)

        if index < 10000:
            # Phosphorus Model
            if Qin >= 0.0001:
                # Accumulate time elapsed since water entered node
                t = t + 1
                # Calculate new concentration
                Cnew_dpA = (Cin_dpA*np.exp((-k_dpA*L*A*E)/Qin))+\
                    (Ceq0_dpA*np.exp(B1_dpA*t))*(1-(np.exp((-k_dpA*L*A*E)/Qin)))
                Cnew_ppA = (Cin_ppA*np.exp((-k_ppA*L*A*E)/Qin))+\
                    (Ceq0_ppA*np.exp(B1_ppA*t))*(1-(np.exp((-k_ppA*L*A*E)/Qin)))

                # Set new concentration
                sim._model.setNodePollutant("Media", 0, Cnew_ppA)
                sim._model.setNodePollutant("Media", 1, Cnew_dpA)
                #print("Cnew DP & PP:", Cnew_dpA, Cnew_ppA)
            else:
                Cnew_ppA = Ceq0_ppA
                Cnew_dpA = Ceq0_dpA   
                sim._model.setNodePollutant("Media", 0, Cnew_ppA)
                sim._model.setNodePollutant("Media", 1, Cnew_dpA)
                #print("Cnew2 DP & PP:", Cnew_dpA, Cnew_ppA)
        else:
            Cnew_ppA = Ceq0_ppA
            Cnew_dpA = Ceq0_dpA  
            sim._model.setNodePollutant("Media", 0, Cnew_ppA)
            sim._model.setNodePollutant("Media", 1, Cnew_dpA)

        # Save phosphorus concentrations for the rest of the system 
        runoff_dpA.append(overflow.pollut_quality['DP'])
        runoff_ppA.append(overflow.pollut_quality['PP'])    
        infil_dpA.append(Cnew_dpA)
        infil_ppA.append(Cnew_ppA)
        drained_dpA.append(Cnew_dpA)
        drained_ppA.append(Cnew_ppA)
        exfil_dpA.append(Cnew_dpA)
        exfil_ppA.append(Cnew_ppA)

    sim._model.swmm_end()
    #print("Runoff Error:", sim.runoff_error)
    #print("Routing Error:", sim.flow_routing_error)
    #print("Quality Error:", sim.quality_error)

# Convert flow rates from CFS to LPS
conv_cfs_lps = [28.3168]*len(inflow)
inflow_mA = [a*b for a,b in zip(inflowA,conv_cfs_lps)]
runoff_mA = [a*b for a,b in zip(runoffA,conv_cfs_lps)]
infil_mA = [a*b for a,b in zip(infil_rateA,conv_cfs_lps)]
drain_mA = [a*b for a,b in zip(drain_flowA,conv_cfs_lps)]
exfil_mA = [a*b for a,b in zip(exfil_rateA,conv_cfs_lps)]

# Convert height from ft to m
conv_ft_m = [0.3048]*len(inflow)
pond_height_mA = [a*b for a,b in zip(pond_heightA,conv_ft_m)]
soil_height_mA = [a*b for a,b in zip(soil_heightA,conv_ft_m)]

# Calculate TTP concentration
inflow_tpA = [a+b for a,b in zip(inflow_dpA, inflow_ppA)]
runoff_tpA = [a+b for a,b in zip(runoff_dpA, runoff_ppA)]
infil_tpA = [a+b for a,b in zip(infil_dpA, infil_ppA)]
drain_tpA = [a+b for a,b in zip(drained_dpA, drained_ppA)]
exfil_tpA = [a+b for a,b in zip(exfil_dpA, exfil_ppA)]

# Calculate DP load (g) each timestep
conv_mgs_gs = [0.001]*len(inflow)
timestep = [5]*len(inflow)
inflow_load_dpA = [a*b*c*d for a,b,c,d in zip(inflow_dpA,inflow_mA,conv_mgs_gs,timestep)]
runoff_load_dpA = [a*b*c*d for a,b,c,d in zip(runoff_dpA,runoff_mA,conv_mgs_gs,timestep)]
infil_load_dpA = [a*b*c*d for a,b,c,d in zip(infil_dpA,infil_mA,conv_mgs_gs,timestep)]
drain_load_dpA = [a*b*c*d for a,b,c,d in zip(drained_dpA,drain_mA,conv_mgs_gs,timestep)]
exfil_load_dpA = [a*b*c*d for a,b,c,d in zip(exfil_dpA,exfil_mA,conv_mgs_gs,timestep)]

# Calculate DP load (g) each timestep
inflow_load_ppA = [a*b*c*d for a,b,c,d in zip(inflow_ppA,inflow_mA,conv_mgs_gs,timestep)]
runoff_load_ppA = [a*b*c*d for a,b,c,d in zip(runoff_ppA,runoff_mA,conv_mgs_gs,timestep)]
infil_load_ppA = [a*b*c*d for a,b,c,d in zip(infil_dpA,infil_mA,conv_mgs_gs,timestep)]
drain_load_ppA = [a*b*c*d for a,b,c,d in zip(drained_ppA,drain_mA,conv_mgs_gs,timestep)]
exfil_load_ppA = [a*b*c*d for a,b,c,d in zip(exfil_ppA,exfil_mA,conv_mgs_gs,timestep)]

# Calculate cumulative DP load (g)
inflow_cumload_dpA = np.cumsum(inflow_load_dpA)
runoff_cumload_dpA = np.cumsum(runoff_load_dpA)
infil_cumload_dpA = np.cumsum(infil_load_dpA)
drain_cumload_dpA = np.cumsum(drain_load_dpA)
exfil_cumload_dpA = np.cumsum(exfil_load_dpA)

# Calculate cumulative PP load (g)
inflow_cumload_ppA = np.cumsum(inflow_load_ppA)
runoff_cumload_ppA = np.cumsum(runoff_load_ppA)
infil_cumload_ppA = np.cumsum(infil_load_ppA)
drain_cumload_ppA = np.cumsum(drain_load_ppA)
exfil_cumload_ppA = np.cumsum(exfil_load_ppA)

# Calculate TP load (g)
influent_loadA = [a+b for a,b in zip(inflow_load_dpA, inflow_load_ppA)]
runoff_loadA = [a+b for a,b in zip(runoff_load_dpA, runoff_load_ppA)]
infil_loadA = [a+b for a,b in zip(infil_load_dpA, infil_load_ppA)]
drain_loadA = [a+b for a,b in zip(drain_load_dpA, drain_load_ppA)]
released_loadA = [a+b+c+d for a,b,c,d in zip(runoff_load_dpA, runoff_load_ppA, drain_load_dpA, drain_load_ppA)]
exfil_loadA = [a+b for a,b in zip(exfil_load_dpA, exfil_load_ppA)]
captured_loadA = [a-b for a,b in zip(influent_loadA, released_loadA)]

# Calculate cumulative TP load (g)
influent_cumloadA = np.cumsum(influent_loadA)
released_cumloadA = np.cumsum(released_loadA)
captured_cumloadA = np.cumsum(captured_loadA)

#----------------------------------------------------------------------#

# Print final DP and PP loads (g)
print("Inflow:", influent_cumload[-1])
print("InflowC:", influent_cumload[-1])
print("InflowA:", influent_cumload[-1])
print("Released:", released_cumload[-1])
print("ReleasedC:", released_cumloadC[-1])
print("ReleasedA:", released_cumloadA[-1])
print("Captured:", captured_cumload[-1])
print("CapturedC:", captured_cumloadC[-1])
print("CapturedA:", captured_cumloadA[-1])

#----------------------------------------------------------------------#

# Plot Depths Results
fig, ax = plt.subplots(2, 1)
ax[0].plot(pond_height_mC, color='#695580', linewidth=2, label="RTC")
ax[0].plot(pond_height_m, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[0].plot(pond_height_mA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[0].set_xticks([])
ax[0].set_ylabel("Pond Height (m)")
ax[0].legend()

ax[1].plot(soil_height_mC, color='#695580', linewidth=2, label="RTC")
ax[1].plot(soil_height_m, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[1].plot(soil_height_mA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[1].set_xticks([])
ax[1].set_ylabel("Soil Height (m)")
ax[1].legend()
plt.show()

# Plot Flow and Concentration Results
fig, ax = plt.subplots(5, 2)
ax[0,0].plot(inflow_mC, color='#695580', linewidth=2, label="RTC")
ax[0,0].plot(inflow_m, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[0,0].plot(inflow_mA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[0,0].set_xticks([])
ax[0,0].set_ylabel("Inflow (L/s)")
ax[0,0].legend()

ax[1,0].plot(runoff_mC, color='#695580', linewidth=2, label="RTC")
ax[1,0].plot(runoff_m, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[1,0].plot(runoff_mA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[1,0].set_xticks([])
ax[1,0].set_ylabel("Runoff (L/s)")

ax[2,0].plot(infil_mC, color='#695580', linewidth=2, label="RTC")
ax[2,0].plot(infil_m, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[2,0].plot(infil_mA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[2,0].set_xticks([])
ax[2,0].set_ylabel("Infiltration (L/s)")

ax[3,0].plot(drain_mC, color='#695580', linewidth=2, label="RTC")
ax[3,0].plot(drain_m, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[3,0].plot(drain_mA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[3,0].plot(valve_setting, ':', color='#818282', linewidth=2, label="Valve Setting")
ax[3,0].set_xticks([])
ax[3,0].set_ylabel("Drain Flow (L/s)")
ax[3,0].legend()

ax[4,0].plot(exfil_mC, color='#695580', linewidth=2, label="RTC")
ax[4,0].plot(exfil_m, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[4,0].plot(exfil_mA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[4,0].set_ylabel("Exfiltration (L/s)")
#ax[4,0].set_xticks([0,7374,14747,22120,29494,36867,44241])
#ax[4,0].set_xticklabels(["0","6","12","18","24","30","36"])
ax[3,0].set_xlabel("Time (hours)")

ax[0,1].plot(influent_loadC, color='#695580', linewidth=2, label='RTC')
ax[0,1].plot(influent_load, '-.', color='#6CC6D1', linewidth=2, label='Standard')
ax[0,1].plot(influent_loadA, '--', color='#3B4D7A', linewidth=2, label='Amended')
ax[0,1].set_ylabel("Influent Load (g)")
ax[0,1].set_xticks([])
ax[0,1].legend()

ax[1,1].plot(runoff_tpC, color='#695580', linewidth=2, label='RTC')
ax[1,1].plot(runoff_tp, '-.', color='#6CC6D1', linewidth=2, label='Standard')
ax[1,1].plot(runoff_tpA, '--', color='#3B4D7A', linewidth=2, label='Amended')
ax[1,1].set_ylabel("Soil/Drain Conc (mg/L)")
ax[1,1].set_xticks([])

ax[2,1].plot(infil_tpC, color='#695580', linewidth=2, label='RTC')
ax[2,1].plot(infil_tp, '-.', color='#6CC6D1', linewidth=2, label='Standard')
ax[2,1].plot(infil_tpA, '--', color='#3B4D7A', linewidth=2, label='Amended')
ax[2,1].set_ylabel("Soil/Drain Conc (mg/L)")
ax[2,1].set_xticks([])

ax[3,1].plot(drain_loadC, color='#695580', linewidth=2, label='RTC')
ax[3,1].plot(drain_load, '-.', color='#6CC6D1', linewidth=2, label='Standard')
ax[3,1].plot(drain_loadA, '--', color='#3B4D7A', linewidth=2, label='Amended')
ax[3,1].set_ylabel("Drain Load (g)")
ax[3,1].set_xticks([])

ax[4,1].plot(exfil_loadC, color='#695580', linewidth=2, label='RTC')
ax[4,1].plot(exfil_load, '-.', color='#6CC6D1', linewidth=2, label='Standard')
ax[4,1].plot(exfil_loadA, '--', color='#3B4D7A', linewidth=2, label='Amended')
ax[4,1].set_ylabel("Caputed Load (g)")
#ax[4,1].set_xticks([0,7373,14747,22120,29494,36867,44241])
#ax[4,1].set_xticklabels(["0","6","12","18","24","30","36"])
ax[4,1].set_xlabel("Time (hours)")
plt.show()

# Plot Cumulative Load Result
fig, ax = plt.subplots(1,3)
ax[0].plot(influent_cumloadC, color='#695580', linewidth=2, label="RTC")
ax[0].plot(influent_cumload, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[0].plot(influent_cumloadA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[0].set_ylabel("TP Cum Load Recieved (g)")
#ax[0].set_xticks([0,7373,14747,22120,29494,36867,44241])
#ax[0].set_xticklabels(["0","6","12","18","24","30","36"])
ax[0].set_xlabel("Time (hours)")
ax[0].legend()

ax[1].plot(released_cumloadC, color='#695580', linewidth=2, label="RTC")
ax[1].plot(released_cumload, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[1].plot(released_cumloadA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[1].set_ylabel("TP Cum Load Released(g)")
#ax[1].set_xticks([0,7373,14747,22120,29494,36867,44241])
#ax[1].set_xticklabels(["0","6","12","18","24","30","36"])
ax[1].set_xlabel("Time (hours)")

ax[2].plot(captured_cumloadC, color='#695580', linewidth=2, label="RTC")
ax[2].plot(captured_cumload, '-.', color='#6CC6D1', linewidth=2, label="Standard")
ax[2].plot(captured_cumloadA, '--', color='#3B4D7A', linewidth=2, label="Amended")
ax[2].set_ylabel("TP Cum Load Captured(g)")
#ax[2].set_xticks([0,7373,14747,22120,29494,36867,44241])
#ax[2].set_xticklabels(["0","6","12","18","24","30","36"])
ax[2].set_xlabel("Time (hours)")
plt.show()