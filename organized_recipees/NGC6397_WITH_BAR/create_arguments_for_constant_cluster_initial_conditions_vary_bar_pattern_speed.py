"""
This script will create a txt file
each row of that txt file will be a different set of input parameters for the script execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py
the sbatch array will read from that txt file and execute the script for each row
therefore, all the jobs will be independent and run in parallel on different processors
"""
import numpy as np 

MAX_SIMULATION_CALL = 999
fname = "arguments_for_constant_cluster_initial_conditions_vary_bar_pattern_speed.txt"
# clear the file
with open(fname,"w") as f:
    f.write("")

# the pattern speeds to explore
### THIS NEEDS TO BE CONSISTENT WITH THE SCRIPT execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py
omega_min, omega_max = 25,66
omega_step = 0.25
bar_pattern_speeds = np.arange(omega_min,omega_max+omega_step,omega_step)

# we pick the number of initial conditions we want to explore
N_min_monte_carlo_indices = 0
N_max_monte_carlo_indices = 1
step = 1  
monte_carlo_indices = np.arange(N_min_monte_carlo_indices,N_max_monte_carlo_indices+1,step)

# we pick the particle count we want to explore
N_particles = [int(32e3), int(33e3),int(35e3)] # i'm afraid to hit 1e5 particles, so I divide it into a few batches

N_SIMULATIONS = len(bar_pattern_speeds)*len(monte_carlo_indices)*len(N_particles)

if N_SIMULATIONS > MAX_SIMULATION_CALL:
    print("The number of simulations is too high, {:d} > {:d}".format(N_SIMULATIONS,MAX_SIMULATION_CALL))
    print("Please reduce the number of simulations")
    exit()
print("We will run {:d} simulations".format(N_SIMULATIONS))
# now we will write the file of all the input arguments 
for monte_carlo_index in monte_carlo_indices:
    for j in range(len(bar_pattern_speeds)):
        for NP in N_particles:
            with open(fname,"a") as f:
                f.write("{:d} {:d} {:d}\n".format(monte_carlo_index,j,NP))
print(fname, "written")


