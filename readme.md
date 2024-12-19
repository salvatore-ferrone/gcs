# What's the big idea?

This repository is in-essence a series of recipees for calling tstrippy to produce simulation data. Additionally, I have the functions that handle I/O to keep the file format standard. 

This repository IS NOT INTENDED FOR ANALYSIS. The main reason for this is that f2py is depreciated. So more recent versions of scipy, matplotlib, etc rely on more recent versions of numpy, which dropped f2py after 1.22. 


There are two classe of simulations, orbits and snapshots. An orbit file will save the entire evolution of a particle into an HDF5 file. A snapshop will save a group of particles at the same instance into an HDF5 file. Orbits are intended for the globular clusters or dark matter sub haloes where the snapshots are intended to be for the streams with a high particle number. The advantage of separating is that I can integrate the streams until today without saving intermediate timestamps. Otherwise this would be a ton of data. Also, I need the orbit of the globular cluster during the integration of the streams. This is why I need to save its data at many timestamps.

## What does the user need to do?
You must modify the paths.yaml file depending on what machine you're on. Here's an example while I am on nhampi. This will need to change if you are running the simulations on another place.

```yaml
simulations: /home-filer/sferrone/tychoOUTPUTS/
temporary: /home-filer/sferrone/temp/
```
## Orbits
- intended to be for globular cluster orbits. In essence it will be a time series of positions for one object

## Streams
- intended to an ensemble of objects, i.e. a stream. This saves snapshots of particles together. 

