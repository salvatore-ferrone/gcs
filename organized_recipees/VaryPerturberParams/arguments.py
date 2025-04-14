import numpy as np 
import varyPerturberParams

GAURD = 999 # cannot run more than 999 simulations at once 

# number of simulations for the internal grid
N_internal=varyPerturberParams.SIZE_GRID*varyPerturberParams.SIZE_GRID

# number of simulations for the perturber grid
mass_radius=np.loadtxt("sub_halo_mass_radius.txt",dtype=float)
N_perturbers = mass_radius.shape[0]


# WE ARE ALSO ITERATING OVER THE PARTICLES
# compute the number of particle batches 
DIVISOR = GAURD//(N_perturbers*N_internal) 
TOTAL = 1e5
STEP=1
NPs = np.arange(TOTAL//DIVISOR - DIVISOR//2  ,
                TOTAL//DIVISOR + DIVISOR//2 ,
                STEP,dtype=int)

N_particle_batches = len(NPs)

outputfile="arguments.txt"

with open(outputfile,"w") as f:
    f.write("# n_particles internal_index perturber_index\n")


N_batches = len(NPs)*N_internal*N_perturbers
for i in range(N_particle_batches):
    for j in range(N_internal):
        for k in range(N_perturbers):
            with open(outputfile,"a") as f:
                f.write( "{:d} {:d} {:d} \n".format(NPs[i],j,k) )

print("can run ",N_batches," simulations at once")
# loop over perturber index 
