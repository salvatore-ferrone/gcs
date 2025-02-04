import gcs
import h5py
import tstrippy
import numpy as np 
from astropy import units as u
import matplotlib.pyplot as plt
import datetime

author = "Salvatore Ferrone"
author_affiliation = "Sapienza University of Rome"
author_email = "salvatore.ferrone@uniroma1.it"
description = "Integrating bouldrini subhaloes in the MW of pouliasis2017pii"


def main():
    """ 
    Take the subhaloes from Bouldrini, integrate them and store the output. 
    """
    # name this population of sub haloes 
    PopulationName = "Bouldrini2025_SET0"

    #### SET THE INTEGRATION TIME AND TIME STEP ####
    T0 = 0
    integration_time = 5e9 * u.yr
    dt = 1e4 * u.yr
    NSTEP = int(integration_time/dt)
    # convert to integration units 
    unitT = u.s * u.kpc/u.km
    unitV = u.km/u.s
    unitM = u.Msun
    timestep = dt.to(unitT).value
    integration_time = integration_time.to(unitT).value
    time_stamps = -np.linspace(T0,integration_time,NSTEP+1)


    # load the sub halo data
    DMpathname = gcs.path_handler.BouldriniDarkmatterSubloes()
    with h5py.File(DMpathname, 'r') as DMhaloes:
        xh,yh,zh=DMhaloes['xh'][:],DMhaloes['yh'][:],DMhaloes['zh'][:]
        vxh,vyh,vzh=DMhaloes['vxh'][:],DMhaloes['vyh'][:],DMhaloes['vzh'][:]
        vxh,vyh,vzh=vxh/100,vyh/100,vzh/100 # came in units of 100 km/s
        mh,rs=DMhaloes['mh'][:],DMhaloes['rs'][:]
    Nparticles=len(xh)
    
    # load the MW potential
    MWpotential =   "pouliasis2017pii"
    MWparams    =   tstrippy.Parsers.potential_parameters.pouliasis2017pii()

    # pick the halo mass model 
    mass_profile = m_enc_martos
    MF_params = MWparams[1:5] # the parameters of the MW halo in our potential

    # unnormalize the velocities 
    vxh,vyh,vzh = get_velocities_for_a_dm_halo(mass_profile,MF_params,[vxh,vyh,vzh],[xh,yh,zh])

    starttime=datetime.datetime.now()
    # initialize the integrator
    tstrippy.integrator.setinitialkinematics(xh,yh,zh,vxh,vyh,vzh)
    tstrippy.integrator.setstaticgalaxy("pouliasis2017pii", MWparams)
    tstrippy.integrator.setintegrationparameters(T0,timestep,NSTEP)
    tstrippy.integrator.setbackwardorbit()
    # perform the integration 
    xt,yt,zt,vxt,vyt,vzt = tstrippy.integrator.leapfrogintime(NSTEP,Nparticles)
    tstrippy.integrator.deallocate()
    endtime=datetime.datetime.now()
    print("Elapsed time: ", endtime-starttime)

    # flip the orbits and the sign of the velocities
    xt,yt,zt,vxt,vyt,vzt = xt[::-1],yt[::-1],zt[::-1],-vxt[::-1],-vyt[::-1],-vzt[::-1]
    # flip the time stamps 
    time_stamps = time_stamps[::-1]

    ## pack up attributes
    attributes = dict(
        author = author,
        author_affiliation = author_affiliation,
        author_email = author_email,
        description = description,
        integration_time = integration_time,
        dt = timestep,
        NSTEP = NSTEP,
        timestep = timestep,
        Nparticles = Nparticles,
        MWpotential = MWpotential,
    )
    # save the output
    outname = gcs.path_handler.Orbits_DarkMatterSubhalo(PopulationName,MWpotential)
    with h5py.File(outname, 'w') as f:
        f.create_dataset("xt", data=xt,dtype=np.float32)
        f.create_dataset("yt", data=yt,dtype=np.float32)
        f.create_dataset("zt", data=zt,dtype=np.float32)
        f.create_dataset("vxt", data=vxt,dtype=np.float32)
        f.create_dataset("vyt", data=vyt,dtype=np.float32)
        f.create_dataset("vzt", data=vzt,dtype=np.float32)
        f.create_dataset("time_stamps", data=time_stamps,dtype=np.float32)
        f.create_dataset("Mass", data=mh,dtype=np.float32)
        f.create_dataset("rs", data=rs,dtype=np.float32)
        for key in attributes:
            f.attrs[key] = attributes[key]
    print("Done with", outname)
    return None




# the halo from martos et al 2016
def m_enc_martos(params,x,y,z):
    Mstar,a,gamma,rcut = params
    r = np.sqrt(x**2 + y**2 + z**2)
    Menc = Mstar * (r/a)**gamma / (1 + (r/a)**(gamma-1))
    beyond = r > rcut
    Menc[beyond] = Mstar * (rcut/a)**gamma / (1 + (rcut/a)**(gamma-1))
    return Menc

def get_velocities_for_a_dm_halo(mass_profile,MF_params,velocities,positions):
    """
        This function will take the positions and velocities of the subhaloes and put them in the MW potential
    """
    x,y,z = positions
    vx,vy,vz = velocities # givein in 100 k/s / sqrt(Menc)
    menc = mass_profile(MF_params,x,y,z)
    vx = vx * np.sqrt(menc)
    vy = vy * np.sqrt(menc)
    vz = vz * np.sqrt(menc)
    return vx,vy,vz





if __name__ == "__main__":
    main()