#!/bin/bash



##### This scrip submits multiple jobs to the slurm queue with dependencies.
####      First, we integrate one set of GC initial conditions and parallelize over the bar pattern speeds
####          once a set of initial conditions are done, we move onto the next set of initial conditions. 
####        We can only queue 5 at a time

### made the script executable with:
###### chmod +x /obs/sferrone/gcs/recipees/slurm_loop_monte_carlo_over_constant_cluster_initial_conditions_vary_bar_pattern_speed.sh
#### Run the script with this: 

#################  EXECUTE ################# EXECUTE ################# EXECUTE ################# EXECUTE ################# EXECUTE #################
## nohup ./slurm_loop_monte_carlo_over_constant_cluster_initial_conditions_vary_bar_pattern_speed.sh > slurm_loop_monte_carlo.log 2>&1 &
#################  EXECUTE ################# EXECUTE ################# EXECUTE ################# EXECUTE ################# EXECUTE #################


#### MOVE TO NHAMPI 
#### scp -r /scratch2/sferrone/simulations/StreamEvolutionInBarredPotential/pouliasis2017pii/longmuralibar/Pal5/5000/ sferrone@nhampi:/home-filer/sferrone/tychoOUTPUTS/StreamEvolutionInBarredPotential/pouliasis2017pii/longmuralibar/Pal5/


# Number of iterations over the initial conditions
OUTER_LOOP=5

# Number of bar pattern speeds. THIS IS NOT arbitrary but has to be consistent with what's inside the python script: `execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.py`
#   This is consistent with the number of pattern speeds from 
#       omega_min, omega_max = 25,66
#       omega_step = 0.25
#       bar_pattern_speeds = np.arange(omega_min,omega_max+omega_step,omega_step)
INNER_LOOP=165



# Loop over the outer loop
for i in $(seq 0 $((OUTER_LOOP-1)))
do
    # Submit the inner loop as an array job
    jobid=$(sbatch --array=0-$((INNER_LOOP-1)) --export=ALL,OUTER_INDEX=$i slurm_execute_constant_cluster_initial_conditions_vary_bar_pattern_speed.sh | awk '{print $4}')
    echo "Submitted array job $jobid with OUTER_INDEX=$i"

done