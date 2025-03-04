# making the input parameters for plummer_pal5_mass_radius_grid.py so that I can launch the jobs in parallel
# I am looping over the particle number, mass, and radius 
import numpy as np 



# MASS GRID AND RADIUS GRID
SIZE_GRID = 5
N_MASS_SAMPLING = SIZE_GRID
N_RADIUS_SAMPLING = SIZE_GRID # square grid
MASS_GRID = np.logspace(4,5.2,N_MASS_SAMPLING) # in Msun
RADIUS_GRID = np.linspace(5,30,N_RADIUS_SAMPLING)/1000 # in kpc


# WE ARE ALSO ITERATING OVER THE PARTICLES
start = 9000
stop = 9900
step = 100
mysequence = np.arange(start,stop+step,step,dtype=int)
mysequence[0]=mysequence[0]//2

# OUTPUT FILE
output_file = "mass_radius_grid_params.txt"
for i in range(SIZE_GRID):
    for j in range(SIZE_GRID):
        for n_particles in mysequence:
            with open(output_file,"a") as f:
                f.write(f"{i} {j} {n_particles}\n")
