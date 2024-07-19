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


def dt_in_years_to_integration_units(dt, T = 5e9*u.yr, unitT = (u.s * u.kpc / u.km)):
    assert(isinstance(dt,u.Quantity))
    assert(isinstance(T,u.Quantity))

    dt = dt.to(u.yr)
    T = T.to(u.yr)
    NSTEP = int((T/dt).value)
    T = T.to(unitT)
    return np.linspace(0,T.value,NSTEP)
    