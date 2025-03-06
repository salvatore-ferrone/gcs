# the scientific question
The stream associated with NGC6397 is not as prominent as Boldrini \& Virtrol 2019 nor Arnold 
\& Baumgardt 2025 suggest it should be. Ibata et al 2019 observe a portion of a stream, however, it is quite limited. I hypothesize that the galactic bar is reponsible for destroying the stream that should originate from this cluster. In this folder, we test that hypothesis. Specifically, we probe the parameter space of the galactic bar: Mass, length, and pattern speed, as well as varying the initial conditions of NGC6397 based on the uncertanties in the baumgardt catalog. 


The folder essentially has the base script, and then other python codes that are wrappers over this code, fixing some parameters and iterating over others, depending on the experiment. Then there are the `slurm` wrappers that call the python code when submitting to the jobs to the cluster. There are also functions that generate the input parameters for the slurm script to facilitate launching the independent jobs. 

## The base code
`stream_evolution_in_barred_potential.py`
The input parameters are many: 
    - cluster_initial_conditions
        - RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m of the globular cluster

    - particle_initial_conditions
        - the x,y,z,vx,vy,vz positions of each particle
    
    - montecarlokey
        - the monte-carlo-sampling of the uncertainties for the host cluster

    - GCname
        - the name of the host clusters

    - MWpotential
        - the name of the MW. only `pouliasis2017pii` is supported

    - internal_dynamics
        - only `isotropic-plummer` is currently supported

    - barname
        - currently only `longmurali` is supported

    - barparams
        - the [Mass, x-axis, y-axis, z-axis]

    - barpoly
        - a polynomial describing the rotation of the bar. the 0th argument is the initial inclination, the 1st is the pattern speed, the 2nd is the acceleration, and so on. 

    - integrationtime
        - how long to go backward given as a number with astropy units
    
    - dt
        - the time step 
    
    - NSKIP
        - the number of time steps to skip when saving snap shots
    
    - temp_base_name
        - the name of the snapshots for the temporary files 
    
    - description
        - a string that will be added to the attributes of the output file describing the experiment 

    - writestream
        - a boolean dictating if intermediate snapshots will be saved or not. 

The code does the following: 

    1. Computes the orbit of the host globular cluster as a point mass in a barred potential 

    2. places a plummer sphere at the posiiton 5 gyr ago 

    3. integrates the plummer sphere to create the stream

    4. saves the output data in a file that holds the stream snapshots, host orbit, and the final stream snapshot

