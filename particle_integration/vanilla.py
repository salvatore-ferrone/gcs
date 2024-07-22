import tstrippy
from gcs import path_handler as ph
import gcs
import numpy as np 
import astropy.units as u


def load_arguments(GCname,montecarlokey,internal_dynamics,GCorbits_potential,MWpotential,NP,T,dt):
    
    
    assert(type(GCname) == str)
    assert(type(montecarlokey) == str)
    assert(type(internal_dynamics) == str)
    assert(type(GCorbits_potential) == str)
    assert(type(MWpotential) == str)
    assert(type(NP) == int)
    assert(type(dt) == u.quantity.Quantity)
    assert(dt.unit == u.yr)
    
    # the time stuff
    T,dt,Nstep,tsampling =gcs.misc.get_time_sampling_from_years_to_integration_units(T=T,dt=dt)

    
    # load the particle positions
    fname           =   ph.ParticleInitialConditions(GCname)
    x,y,z,vx,vy,vz  =   gcs.extractors.ParticleInitialConditions.load_particles(fname,internal_dynamics,montecarlokey,NP)
    
    # Extract the orbit  
    orbit_file_name             =   ph.GC_orbits(GCorbits_potential,GCname)
    tH,xH,yH,zH,vxH,vyH,vzH     =   gcs.extractors.GCOrbits.extract_whole_orbit(orbit_file_name,montecarlokey)
    # interpolate the host orbit to the time resolution of the particle
    xHost,yHost,zHost,vxHost,vyHost,vzHost        =   gcs.misc.interpolate_finer_grid(tsampling,tH,xH,yH,zH,vxH,vyH,vzH)
    
    ## load the host mass and size
    Mass,rh_m,_,_,_,_,_,_=gcs.extractors.MonteCarloObservables.extract_all_GC_observables([GCname],montecarlokey)
    rplummer=gcs.misc.half_mass_to_plummer(rh_m[0]).value
    mass_host = Mass[0].value
    
    # get milky way params
    MWparams = tstrippy.Parsers.potential_parameters.pouliasis2017pii()
    
    # organize the input arguments for the integrator
    staticgalaxy = (MWpotential,MWparams)
    integrationparameters = (T.value,dt.value,Nstep)
    initialkinematics = (x+xHost[0],y+yHost[0],z+zHost[0],vx+vxHost[0],vy+vyHost[0],vz+vzHost[0])
    inithostperturber = (tsampling,xHost,yHost,zHost,vxHost,vyHost,vzHost,mass_host,rplummer)
    
    return staticgalaxy,integrationparameters,initialkinematics,inithostperturber



def initialize_vanilla_integrator(staticgalaxy,integrationparameters,initialkinematics,inithostperturber):
    integrator = tstrippy.integrator
    integrator.setstaticgalaxy(*staticgalaxy)
    integrator.setintegrationparameters(*integrationparameters)
    integrator.setinitialkinematics(*initialkinematics)
    integrator.inithostperturber(*inithostperturber)
    return integrator



def initialize_write_snapshots(integrator,NSKIP,GCname,temporary_directory):
    assert type(NSKIP) == int
    assert type(GCname) == str
    assert type(temporary_directory) == str
    
    integrator.initwritestream(NSKIP,GCname,temporary_directory)
    
    return integrator
    


def leapfrogtofinalpositions(integrator):
    
    integrator.leapfrogtofinalpositions()
    # integrate the particle
    xf  = integrator.xf
    yf  = integrator.yf
    zf  = integrator.zf
    vxf = integrator.vxf
    vyf = integrator.vyf
    vzf = integrator.vzf
    tesc= integrator.tesc
    phase_space = np.array([xf,yf,zf,vxf,vyf,vzf])
    integrator.deallocate()
    return phase_space,tesc




if __name__ == "__main__" : 
    import sys 
    import os
    GCname              =   "NGC104"
    montecarlokey       =   "monte-carlo-000"
    internal_dynamics   =   "isotropic-plummer"
    GCorbits_potential  =   "pouliasis2017pii"
    MWpotential         =   "pouliasis2017pii"
    NP                  =   int(1e2)
    T                   =   5e9*u.yr
    dt                  =   1e5*u.yr
    NSKIP               =   int(100)
    
    attributes = {
        "NP":NP,
        "T":T,
        "dt":dt,
        "NSKIP":NSKIP,
        "GCname":GCname,
        "montecarlokey":montecarlokey,
        "internal_dynamics":internal_dynamics,
        "GCorbits_potential":GCorbits_potential,
        "MWpotential":MWpotential
    }
    
    tempdir=ph._TemporaryStreamSnapShots(MWpotential,GCname,NP,internal_dynamics,montecarlokey)
    filename=ph.Stream(GCname,NP,MWpotential,internal_dynamics,montecarlokey)
    
    if os.path.exists(filename):
        print(filename, "Already exists. \n Skipping!")
        sys.exit(0)

    staticgalaxy,integrationparameters,initialkinematics,inithostperturber = \
        load_arguments(GCname,montecarlokey,internal_dynamics,GCorbits_potential,MWpotential,NP,T,dt) 


    ###############################################
    ########### INITIALIZE THE INTEGRATOR #########
    ################################
    integrator=initialize_vanilla_integrator(staticgalaxy,integrationparameters,initialkinematics,inithostperturber)
    integrator=initialize_write_snapshots(integrator,NSKIP,GCname,tempdir)
    
    ###############################################
    ########### PERFORM THE INTEGRATION ###########
    ###############################################
    phase_space,tesc = leapfrogtofinalpositions(integrator)
    ##### SAVE THE FINAL SNAP SHOT 
    gcs.writers.Stream.stream(filename,phase_space,tesc,attributes)
    print(filename, "saved")
    
    ### ASSEMBLE THE INTERMEDIATE SNAPSHOTS INTO ONE FILE
    T,dt,Nstep,tsampling=gcs.misc.get_time_sampling_from_years_to_integration_units(T=T,dt=dt)
    snapshottimesampling = tsampling[::NSKIP]
    outfilename = ph.StreamSnapShots(GCname,NP,MWpotential,internal_dynamics,montecarlokey)
    gcs.writers.Stream.StreamSnapShots(outfilename,snapshottimesampling,attributes,tempdir)
    print(outfilename, "saved")
    