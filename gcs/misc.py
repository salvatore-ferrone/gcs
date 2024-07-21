import numpy as np
from astropy import units as u

def half_mass_to_plummer(half_mass):
    """Convert half-mass radius to Plummer radius."""
    factor = ( (1/2)**(-2/3) - 1 )**0.5
    return half_mass*factor


def interpolate_finer_grid(tnew,t,x,y,z,vx,vy,vz):
    """perform a linear interpolation of x(t) at tnew"""
    xx = np.interp(tnew,t,x)
    yy  = np.interp(tnew,t,y)
    zz = np.interp(tnew,t,z)
    vxx = np.interp(tnew,t,vx)
    vyy = np.interp(tnew,t,vy)
    vzz = np.interp(tnew,t,vz)
    return xx,yy,zz,vxx,vyy,vzz


def get_time_sampling_from_years_to_integration_units(T=5e9*u.yr,dt=1e4*u.yr):
    unitT               =   u.s * u.kpc / u.km
    tsampling           =   np.arange(0,T.value+dt.value,dt.value)*u.yr
    tsampling           =   tsampling - T
    tsampling           =   tsampling.to(unitT).value
    Nstep               =   len(tsampling) -1 
    T                   =   T.to(unitT)
    dt                  =   dt.to(unitT)    
    return T,dt,Nstep,tsampling
