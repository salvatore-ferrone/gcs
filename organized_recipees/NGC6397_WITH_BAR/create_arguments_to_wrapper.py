"""
This script will create a txt file
each row of that txt file will be a different set of input parameters for the script execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py
the sbatch array will read from that txt file and execute the script for each row
therefore, all the jobs will be independent and run in parallel on different processors
"""
import numpy as np 
import wrapper_bar_experiment as wbe



MAX_SIMULATION_CALL = 999 # guard for not exceeding the maximum number of simulations that can be run in parallel
fname = "arguments.txt"


# MOST probable values
bar_mass_most_prob = 1.6e10 # kipper et al 2020
bar_length_most_prob = 3.5 # lucey et al 
bar_angle_most_prob = 30*np.pi/180 # 
bar_axis_ratio_most_prob = 1/4 # not sure
pattern_speed_most_prob = 39.5 # Partail 2017 and 2016


MASS_INDEX_MOST_PROB = np.argmin(np.abs(wbe.BAR_MASSES-bar_mass_most_prob))
LENGTH_INDEX_MOST_PROB = np.argmin(np.abs(wbe.BARLENGTHS-bar_length_most_prob))
ANGLE_INDEX_MOST_PROB = np.argmin(np.abs(wbe.BAR_ANGLES-bar_angle_most_prob))
AXIS_RATIO_INDEX_MOST_PROB = np.argmin(np.abs(wbe.AXIS_RATIOS-bar_axis_ratio_most_prob))
PATTERN_SPEED_INDEX_MOST_PROB = np.argmin(np.abs(wbe.PATTERN_SPEEDS-pattern_speed_most_prob))


# number of particles
NPs_start,NPs_nsims,NPs_step = 1e5/5,5,1
NPs = np.arange(NPs_start,NPs_start+NPs_nsims+NPs_step,NPs_step,dtype=int)

# monte carlo 
montecarloindex = 0 # we will only run one monte carlo index right now 

GCname = "NGC6397"


N_SIMULATIONS = len(NPs)*wbe.npatternspeeds

if N_SIMULATIONS > MAX_SIMULATION_CALL:
    print("The number of simulations is too high, {:d} > {:d}".format(N_SIMULATIONS,MAX_SIMULATION_CALL))
    print("Please reduce the number of simulations")
    exit()
print("We will run {:d} simulations".format(N_SIMULATIONS))
print("Up date the sbatch array with the number of simulations")


# clear the file
with open(fname,"w") as f:
    f.write("")


for i in range(len(NPs)):
    for j in range(wbe.npatternspeeds):
        call_string = "{:s} {:d} {:d} {:d} {:d} {:d} {:d} {:d}".format(GCname,montecarloindex, NPs[i], ANGLE_INDEX_MOST_PROB, j, MASS_INDEX_MOST_PROB, LENGTH_INDEX_MOST_PROB, AXIS_RATIO_INDEX_MOST_PROB)
        with open(fname,"a") as f:
            f.write(call_string + "\n")


print(fname, "written")


