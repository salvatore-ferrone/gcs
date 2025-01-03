"""
Functions that are commonly used for organizing data 
"""

import tstrippy
from gcs import path_handler as ph
import gcs
import numpy as np 
import astropy.units as u
from astropy import coordinates as coord
import copy
import yaml

def load_vanilla_arguments(GCname,montecarlokey,internal_dynamics,GCorbits_potential,MWpotential,NP,integrationtime,dt):
    
    assert(type(GCname) == str)
    assert(type(montecarlokey) == str)
    assert(type(internal_dynamics) == str)
    assert(type(GCorbits_potential) == str)
    assert(type(MWpotential) == str)
    assert(type(NP) == int)
    assert(type(dt) == u.quantity.Quantity)
    assert(dt.unit == u.yr)
    assert(type(integrationtime) == u.quantity.Quantity)
    assert(integrationtime.unit == u.yr)
    
    # the time stuff
    NSTEP = int(integrationtime.value/dt.value)
    unitT = u.s*(u.kpc/u.km)
    integrationtime = integrationtime.to(unitT)
    dt = dt.to(unitT)
    T0 = -integrationtime
    tsampling = np.linspace(T0.value,0,NSTEP+1)


    
    # load the particle positions
    fname                       =   ph.ParticleInitialConditions(GCname)
    x,y,z,vx,vy,vz              =   gcs.extractors.ParticleInitialConditions.load_particles(fname,internal_dynamics,montecarlokey,NP)
    
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
    integrationparameters = (T0.value,dt.value,NSTEP)
    initialkinematics = (x+xHost[0],y+yHost[0],z+zHost[0],vx+vxHost[0],vy+vyHost[0],vz+vzHost[0])
    inithostperturber = (tsampling,xHost,yHost,zHost,vxHost,vyHost,vzHost,mass_host,rplummer)
    
    return staticgalaxy,integrationparameters,initialkinematics,inithostperturber


def load_bar_yaml(barfile):
    """
    Load in the bar parameters and append the gravitational constant
    returns: barname,barparams,barpoly
    """
    with open(barfile, 'r') as stream:
        try:
            data_loaded = yaml.safe_load(stream)
            # append the gravitational constant
            G = tstrippy.Parsers.potential_parameters.G
            barparams = data_loaded['bar_parameters']
            barparams.insert(0,G)
            barname=data_loaded['name']
            barpoly = data_loaded['bar_polynomial']
            return barname,barparams,barpoly
        except yaml.YAMLError as exc:
            print(exc)


def adjust_disc_mass(staticgalaxy_,barparameters_):
    """
    Remove some mass from the discs to compensate for the bar mass
    """

    staticgalaxy_bar_adjusted = copy.deepcopy(staticgalaxy_)
    M_disc1=staticgalaxy_[1][5]
    M_disc2=staticgalaxy_[1][8]
    Mbar = barparameters_[0]

    M_disc1_new = M_disc1 - Mbar/2
    M_disc2_new = M_disc2 - Mbar/2

    staticgalaxy_bar_adjusted[1][5] = M_disc1_new
    staticgalaxy_bar_adjusted[1][8] = M_disc2_new
    return staticgalaxy_bar_adjusted


def initialize_vanilla_integrator(staticgalaxy,integrationparameters,initialkinematics,inithostperturber):
    integrator = tstrippy.integrator
    integrator.setstaticgalaxy(*staticgalaxy)
    integrator.setintegrationparameters(*integrationparameters)
    integrator.setinitialkinematics(*initialkinematics)
    integrator.inithostperturber(*inithostperturber)
    return integrator


def initialize_write_snapshots(integrator,NSKIP,GCname,temporary_directory):
    """ 
    if we are goign to save out snapshots, this is the call
    """
    assert type(NSKIP) == int
    assert type(GCname) == str
    assert type(temporary_directory) == str
    
    integrator.initwritestream(NSKIP,GCname,temporary_directory)
    
    return integrator
  

def leapfrogtofinalpositions(integrator):
    """
    Obtains the final stream and escape time and deallocates the integrator
    """
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


def get_random_GC_initial_conditions(GCname):
    """
    Sample the initial conditions of a globular cluster from a multivariate normal distribution
    """

    # get random set of initial conditions 
    GCdata=tstrippy.Parsers.baumgardtMWGCs()
    means,cov=GCdata.getGCCovarianceMatrix(GCname)
    initial_conditions=np.random.multivariate_normal(means,cov,1)
    RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=initial_conditions[0]
    return RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m

def sky_to_galactocentric(RA,DEC,Rsun,RV,mualpha,mu_delta):
    skycoords = coord.SkyCoord(ra=RA*u.deg, dec=DEC*u.deg, distance=Rsun*u.pc, pm_ra_cosdec=mualpha*u.mas/u.yr, pm_dec=mu_delta*u.mas/u.yr, radial_velocity=RV*u.km/u.s)
    refframe=tstrippy.Parsers.potential_parameters.MWreferenceframe()
    galcentric=skycoords.transform_to(refframe)
    x,y,z,vx,vy,vz=galcentric.cartesian.x.to(u.kpc).value,galcentric.cartesian.y.to(u.kpc).value,galcentric.cartesian.z.to(u.kpc).value,galcentric.velocity.d_x.to(u.km/u.s).value,galcentric.velocity.d_y.to(u.km/u.s).value,galcentric.velocity.d_z.to(u.km/u.s).value
    return x,y,z,vx,vy,vz