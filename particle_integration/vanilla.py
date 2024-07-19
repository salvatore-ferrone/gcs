import tstrippy
from gcs import path_handler as ph
import gcs
import os 
import numpy as np 
import astropy.units as u


def load_arguments(GCname,montecarlokey,internal_dynamics,GCorbits_potential,MWpotential,NP):
    
    
    # the time stuff
    T,dt,Nstep,tnew=get_time_sampling()

    
    # load the particle positions
    fname           =   ph.ParticleInitialConditions(GCname)
    x,y,z,vx,vy,vz  =   gcs.extractors.ParticleInitialConditions.load_particles(fname,internal_dynamics,montecarlokey,NP)
    
    # now get the positions at the begining of the orbit 
    orbit_file_name             =   ph.GC_orbits(GCorbits_potential,GCname)
    tH,xH,yH,zH,vxH,vyH,vzH     =   gcs.extractors.GCOrbits.extract_whole_orbit(orbit_file_name,montecarlokey)
    # interpolate the host orbit to the time resolution of the particle
    xHost,yHost,zHost,vxHost,vyHost,vzHost        =   gcs.misc.interpolate_finer_grid(tnew,tH,xH,yH,zH,vxH,vyH,vzH)
    
    ## load the host mass and size
    Mass,rh_m,_,_,_,_,_,_=gcs.extractors.MonteCarloObservables.extract_all_GC_observables([GCname],montecarlokey)
    rplummer=gcs.misc.half_mass_to_plummer(rh_m[0]).value
    mass_host = Mass[0].value
    
    # get milky way params
    MWparams = tstrippy.Parsers.potential_parameters.pouliasis2017pii()
    
    # organize output 
    
    staticgalaxy = (MWpotential,MWparams)
    integrationparameters = (T.value,dt.value,Nstep)
    initialkinematics = (x+xHost[0],y+yHost[0],z+zHost[0],vx+vxHost[0],vy+vyHost[0],vz+vzHost[0])
    inithostperturber = (tnew,xHost,yHost,zHost,vxHost,vyHost,vzHost,mass_host,rplummer)
    
    return staticgalaxy,integrationparameters,initialkinematics,inithostperturber
    




def get_time_sampling(T=5e9*u.yr,dt=1e4*u.yr):
    unitT               =   u.s * u.kpc / u.km
    tnew                =   np.arange(0,T.value+dt.value,dt.value)*u.yr
    tnew                =   tnew - T
    tnew                =   tnew.to(unitT).value
    Nstep               =   len(tnew) -1 
    # put all in right time units
    T                   =   T.to(unitT)
    dt                  =   dt.to(unitT)    
    return T,dt,Nstep,tnew


if __name__ == "__main__" : 
    GCname              =   "NGC104"
    montecarlokey       =   "monte-carlo-000"
    internal_dynamics   =   "isotropic-plummer"
    GCorbits_potential  =   "pouliasis2017pii"
    MWpotential         =   "pouliasis2017pii"
    NP                  =   int(1e2)    
    main()