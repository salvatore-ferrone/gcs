# Launching tstrippy from scratch 

**Created**: March 2025. 

**Author**: Salvatore Ferrone

**Target Audience**: Initially for the boss, Paola Di Matteo, and thereafter anyone who wishes to subsequently use this code. 

**Purpose**: Instructions on how install `tstrippy` and `gcs` and then launch simulations. `tstrippy` is the code that creates tidal streams in the Milky Way from globular clusters in the restricted three body problem approximation. `gcs` is my repository with *recipees* that gather the input data, launch the simulations, and store the outputs in an organized manner.

## Instructions 
### Step -2
Have conda configured 

### Step -1
Create a folder where you want to store the code. 
```
mkdir etidal
cd etidal
```

### Step 0 

Clone both `tstrippy` and `gcs` to this directory 
```
git clone https://github.com/salvatore-ferrone/tstrippy.git
git clone https://github.com/salvatore-ferrone/gcs.git
```

You have presumably made it this far if you are reading this readme. 

### Step 1
Configure the output paths. This must be done before installing the conda environment to ensure that the `paths.yaml` goes to the environment directory and doesn't stay local. 
```
cd ./gcs/gcs/
touch paths.yaml
```
Then, edith `paths.yaml` to be
```yaml
simulations: /your/target/directory/
temporary: /where/the/temp/files/will/be/stored/
```
When I work on `tycho`, the cluster at the paris observatory, I place `temporary` right in the poubelle. The files in this directory gets deleted every if not updated after 120 days. This is ok, because none of the simulations that I run have been longer than 20 days. The temporary directory is where snapshots are saved during the integration. This is only used if you are deciding to you save them. Once the simulation is completed, the snapshots are compiled together to make one final hdf5 output, which is placed in the `simulations` directory. 


### Step 2
Install the current environment which handles the proper python version packages. This step is key because unfortunately newer versions of numpy stopped supporting f2py. I haven't figured out how to work around this to still use the fortran code and may never do so. Therefore, it is crucial to have numpy==1.22 installed. 

```
conda env create environment.yml
conda activate gcs
```
This will properly install the dependencies and the two python packages that I created, `tstrippy` and `gcs`. when working in the `gcs` environment, you will be able to import both `tstrippy` and `gcs` in python. 

### Step 3: 
Perform a monte-carlo sampling of the baumgardt catalog. Buamgardt's catalog reports uncertainties on the heliocentric distances, radial velocities, proper motions, the covariance between the proper motions, and the mass for all the globular clusters. The following code samples these uncertainties 1000 times for each globular cluster and stores the output in hdf5 files. These monte-carlo samples are used when computing orbits for the globular clusters and generating a particle distribution. 
``` 
cd initialization
python create_monte_carlo_observables.py
```

### Step 4: 
Launch the specific simulations. My personal recipees are in the `recipees` folder. This may be a bit cluttered for people other than myself. Some simulations require the globular cluster orbits to be precomputed. However, for this delivery, I am showing an example of computing streams in a barred potential. In this case, the globular cluster orbit is computed in a first step, and the stream computed in a second step, and both are saved to the same output file and not separately. I have prepared another directory with one organized recipee for launching many simulations of NGC6397 in a series of galactic potentials varying the bar. 

```
cd ./gcs/organized_recipees/NGC6397_WITH_BAR/
```

### Step 5:
Launch the jobs. This is experiment specific and thus this readme ends her and differs to the readme.md in the `./gcs/organized_recipees/NGC6397_WITH_BAR/` folder. 

