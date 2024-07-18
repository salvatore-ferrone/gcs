from gcs import path_handler as ph 
import h5py 
import numpy as np 


def initial_simulation_coordinates(filename,montecarlokey):
    with h5py.File(filename,'r') as filehdf5:
        x0 = filehdf5[montecarlokey]['xt'][0]
        y0 = filehdf5[montecarlokey]['yt'][0]
        z0 = filehdf5[montecarlokey]['zt'][0]
        vx0 = filehdf5[montecarlokey]['vxt'][0]
        vy0 = filehdf5[montecarlokey]['vyt'][0]
        vz0 = filehdf5[montecarlokey]['vzt'][0]
    return x0,y0,z0,vx0,vy0,vz0


def initial_simulation_coordinates_all(GCnames,potential,montecarlokey):
    assert(isinstance(GCnames,list))
    assert(isinstance(montecarlokey,str))
    assert(isinstance(potential,str))
    phase_space = np.zeros((len(GCnames),6))
    for i in range(len(GCnames)):
        filename = ph.GC_orbits(potential,GCnames[0])
        x0,y0,z0,vx0,vy0,vz0=initial_simulation_coordinates(filename,montecarlokey)
        phase_space[i] = np.array([x0,y0,z0,vx0,vy0,vz0])
    return phase_space