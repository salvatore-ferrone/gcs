from gcs import path_handler as ph 
import h5py 
import numpy as np 


def extract_whole_orbit(filename,montecarlokey):
    """
    t,x,y,z,vx,vy,vz=extract_whole_orbit(filename,montecarlokey)
    """
    with h5py.File(filename,'r') as filehdf5:
        t = filehdf5[montecarlokey]['t'][:]
        x = filehdf5[montecarlokey]['xt'][:]
        y = filehdf5[montecarlokey]['yt'][:]
        z = filehdf5[montecarlokey]['zt'][:]
        vx = filehdf5[montecarlokey]['vxt'][:]
        vy = filehdf5[montecarlokey]['vyt'][:]
        vz = filehdf5[montecarlokey]['vzt'][:]
    return t,x,y,z,vx,vy,vz

def initial_simulation_coordinates(filename,montecarlokey):
    t,x,y,z,vx,vy,vz = extract_whole_orbit(filename,montecarlokey)
    return t[0],x[0],y[0],z[0],vx[0],vy[0],vz[0]


def extract_orbits_from_all_GCS(GCnames,potential,montecarlokey):
    assert(isinstance(GCnames,(list,tuple,np.ndarray))), "GCnames should be a list, tuple or numpy array but was "+str(type(GCnames))
    assert(isinstance(montecarlokey,str))
    assert(isinstance(potential,str))
    
    filename = ph.GC_orbits(potential,GCnames[0])
    
    t,x,y,z,vx,vy,vz = extract_whole_orbit(filename,montecarlokey)
    
    xs = np.zeros((len(GCnames),len(t)))
    ys = np.zeros((len(GCnames),len(t)))
    zs = np.zeros((len(GCnames),len(t)))
    vxs = np.zeros((len(GCnames),len(t)))
    vys = np.zeros((len(GCnames),len(t)))
    vzs = np.zeros((len(GCnames),len(t)))
    
    xs[0],ys[0],zs[0],vxs[0],vys[0],vzs[0] = x,y,z,vx,vy,vz
    for i in range(1,len(GCnames)):
        filename = ph.GC_orbits(potential,GCnames[i])
        _,x0,y0,z0,vx0,vy0,vz0=extract_whole_orbit(filename,montecarlokey)
        xs[i],ys[i],zs[i],vxs[i],vys[i],vzs[i] = x0,y0,z0,vx0,vy0,vz0
    return t,xs,ys,zs,vxs,vys,vzs