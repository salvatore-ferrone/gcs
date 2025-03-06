"""
This script will create a txt file
each row of that txt file will be a different set of input parameters for the script execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py
the sbatch array will read from that txt file and execute the script for each row
therefore, all the jobs will be independent and run in parallel on different processors
"""
import numpy as np 

MAX_SIMULATION_CALL = 999 # guard for not exceeding the maximum number of simulations that can be run in parallel
fname = "arguments_for_constant_cluster_initial_conditions_vary_bar_pattern_speed.txt"
# clear the file
with open(fname,"w") as f:
    f.write("")

# the pattern speeds to explore
# this needs to be consistent with what is in the script execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py
import execute_constant_cluster_initial_conditions_vary_bar_pattern_speed as ECCICVBPS
N_bar_pattern_speeds = ECCICVBPS.bar_pattern_speeds.shape[0]

# we pick the number of initial conditions we want to explore
N_min_monte_carlo_indices = 0
N_max_monte_carlo_indices = 1
step = 1  
monte_carlo_indices = np.arange(N_min_monte_carlo_indices,N_max_monte_carlo_indices+1,step)

# we pick the particle count we want to explore
N_particles = [int(32e3), int(33e3),int(35e3)] # i'm afraid to hit 1e5 particles, so I divide it into a few batches

N_SIMULATIONS = N_bar_pattern_speeds*len(monte_carlo_indices)*len(N_particles)

if N_SIMULATIONS > MAX_SIMULATION_CALL:
    print("The number of simulations is too high, {:d} > {:d}".format(N_SIMULATIONS,MAX_SIMULATION_CALL))
    print("Please reduce the number of simulations")
    exit()
print("We will run {:d} simulations".format(N_SIMULATIONS))

for monte_carlo_index in monte_carlo_indices:
    for j in range(N_bar_pattern_speeds):
        for NP in N_particles:
            with open(fname,"a") as f:
                f.write("{:d} {:d} {:d}\n".format(monte_carlo_index,j,NP))
print(fname, "written")


