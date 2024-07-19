import h5py

def load_particles(fname,internal_dynamics,montecarlokey,NP):
    assert(type(fname) == str)
    assert(type(internal_dynamics) == str)
    assert(type(montecarlokey) == str)
    assert(type(NP) == int)
    with h5py.File(fname,'r') as myfile:
        x=myfile[internal_dynamics][str(int(NP))][montecarlokey]["x"][:]
        y=myfile[internal_dynamics][str(int(NP))][montecarlokey]["y"][:]
        z=myfile[internal_dynamics][str(int(NP))][montecarlokey]["z"][:]
        vx=myfile[internal_dynamics][str(int(NP))][montecarlokey]["vx"][:]
        vy=myfile[internal_dynamics][str(int(NP))][montecarlokey]["vy"][:]
        vz=myfile[internal_dynamics][str(int(NP))][montecarlokey]["vz"][:]
    return x,y,z,vx,vy,vz