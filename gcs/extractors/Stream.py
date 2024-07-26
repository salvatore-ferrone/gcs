from gcs import path_handler as ph 
import h5py 


def extract_old_streams(streampath,internal_dynamics,montecarlokey,NP):
    """This is the format before the july 2024 refactoring. It is kept for compatibility reasons.

    Parameters
    ----------
    streampath : str
        path to the stream file
    internal_dynamics : str
        the identifier of the internal dynamics
    montecarlokey : str
        the index of the monte-carlo key
    NP : int
        number of patciles

    Returns
    -------
    tesc,xp,yp,zp,vx,vy,vz
    np.ndarray
        the escape time and phase space coordinates
    """
    with  h5py.File(streampath, 'r') as myfile:
        tesc=myfile[internal_dynamics][str(NP)][montecarlokey]['tesc'][:]
        xp=myfile[internal_dynamics][str(NP)][montecarlokey]['x'][:]
        yp=myfile[internal_dynamics][str(NP)][montecarlokey]['y'][:]
        zp=myfile[internal_dynamics][str(NP)][montecarlokey]['z'][:]
        vx=myfile[internal_dynamics][str(NP)][montecarlokey]['vx'][:]
        vy=myfile[internal_dynamics][str(NP)][montecarlokey]['vy'][:]
        vz=myfile[internal_dynamics][str(NP)][montecarlokey]['vz'][:]
    return tesc,xp,yp,zp,vx,vy,vz