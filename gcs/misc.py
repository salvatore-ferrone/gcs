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


def unwrap_angles(theta):
    unwrapped_theta = theta.copy()
    correction = 0
    for i in range(1, len(theta)):
        delta = theta[i] - theta[i-1]
        if delta < -np.pi:
            correction += 2*np.pi
        elif delta > np.pi:
            correction -= 2*np.pi
        unwrapped_theta[i] = theta[i] + correction
    return unwrapped_theta




def cartesian_to_spherical(x, y, z, vx, vy, vz):
    # Convert position from Cartesian to spherical coordinates
    r = np.sqrt(x**2 + y**2 + z**2)

    iszero = (r == 0)
    theta= np.zeros_like(r)
    theta[iszero] = 0
    theta[~iszero] = np.arccos(z[~iszero] / r[~iszero])
    phi = np.arctan2(y, x)

    # Compute the spherical velocity components
    rdot = (x * vx + y * vy + z * vz) / r
    
    thetadot= np.zeros_like(r)
    thetadot[iszero] = 0
    thetadot[~iszero] = (x[~iszero] * vx[~iszero] + y[~iszero] * vy[~iszero] - (x[~iszero]**2 + y[~iszero]**2) * vz[~iszero]) / (r[~iszero]**2 * np.sqrt(x[~iszero]**2 + y[~iszero]**2))
    phidot = (x * vy - y * vx) / (x**2 + y**2)



    return r, theta, phi, rdot, thetadot, phidot


def spherical_to_cartesian(r, theta, phi, rdot, thetadot, phidot):
    # Convert position from spherical to Cartesian coordinates
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    # Convert velocity from spherical to Cartesian coordinates
    vx = (rdot * np.sin(theta) * np.cos(phi) +
          r * thetadot * np.cos(theta) * np.cos(phi) -
          r * np.sin(theta) * np.sin(phi) * phidot)
    
    vy = (rdot * np.sin(theta) * np.sin(phi) +
          r * thetadot * np.cos(theta) * np.sin(phi) +
          r * np.sin(theta) * np.cos(phi) * phidot)
    
    vz = rdot * np.cos(theta) - r * thetadot * np.sin(theta)
    
    return x, y, z, vx, vy, vz


def cartesian_to_cylindrical(x, y, z, vx, vy, vz):
    # Convert position from Cartesian to cylindrical coordinates
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    z_cyl = z

    # Compute the cylindrical velocity components
    rhodot = (x * vx + y * vy) / rho 
    phidot = (x * vy - y * vx) / (x**2 + y**2)
    zdot_cyl = vz

    return rho, phi, z_cyl, rhodot, phidot, zdot_cyl


def cylindrical_to_cartesian(rho, phi, z_cyl, rhodot, phidot, zdot_cyl):
    # Convert position from cylindrical to Cartesian coordinates
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    z = z_cyl

    # Convert velocity from cylindrical to Cartesian coordinates
    vx = rhodot * np.cos(phi) - rho * phidot * np.sin(phi)
    vy = rhodot * np.sin(phi) + rho * phidot * np.cos(phi)
    vz = zdot_cyl

    return x, y, z, vx, vy, vz