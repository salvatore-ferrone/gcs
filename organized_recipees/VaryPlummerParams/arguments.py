# making the input parameters for plummer_pal5_mass_radius_grid.py so that I can launch the jobs in parallel
# I am looping over the particle number, mass, and radius 
import numpy as np 
import varyPlummerParams




# WE ARE ALSO ITERATING OVER THE PARTICLES
TOTAL = 1e5
DIVISOR = 30
STEP=1
mysequence = np.arange(TOTAL//DIVISOR - DIVISOR//2 - 1 , 
                       TOTAL//DIVISOR + DIVISOR//2 ,
                       STEP,dtype=int)

# OUTPUT FILE
output_file = "mass_radius_grid_params.txt"
# reset the file
with open(output_file,"w") as f:
    f.write("# mass radius n_particles\n")
# LOOP OVER MASS AND RADIUS
for i in range(varyPlummerParams.SIZE_GRID):
    for j in range(varyPlummerParams.SIZE_GRID):
        for n_particles in mysequence:
            with open(output_file,"a") as f:
                f.write(f"{i} {j} {n_particles}\n")
