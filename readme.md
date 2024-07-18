# What's the big idea?

I am re-doing this folder so I may cleanly reproduce my results from the first and second papers. As simple as possible and as easy reproduce by following a file and then hitting go along the way.


I intend this folder to produce simulations ONLY. No analysis will be provided here! 

There are two classe of simulations, orbits and snapshots. An orbit file will save the entire evolution of a particle into an HDF5 file. A snapshop will save a group of particles at the same instance into an HDF5 file. Orbits are intended for the globular clusters or dark matter sub haloes where the snapshots are intended to be for the streams with a high particle number. The advantage of separating is that I can integrate the streams until today without saving intermediate timestamps. Otherwise this would be a ton of data. Also, I need the orbit of the globular cluster during the integration of the streams. This is why I need to save its data at many timestamps.

## Orbits
- intended to be for globular cluster orbits

## Streams
- intended to 

