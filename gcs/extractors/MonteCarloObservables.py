import h5py
import gcs.path_handler as ph
from astropy import units as u


def extract_GC_observables(filename,montecarlokey,):
    with h5py.File(filename,'r') as filehdf5:
        RA          =   filehdf5[montecarlokey]['RA'][0]        *u.Unit(filehdf5.attrs['RA'])
        DEC         =   filehdf5[montecarlokey]['DEC'][0]       *u.Unit(filehdf5.attrs['DEC'])
        rh_m        =   filehdf5[montecarlokey]['rh_m'][0]      *u.Unit(filehdf5.attrs['rh_m'])
        Mass        =   filehdf5[montecarlokey]['Mass'][0]      *u.Unit(filehdf5.attrs['Mass'])
        Rsun        =   filehdf5[montecarlokey]['Rsun'][0]      *u.Unit(filehdf5.attrs['Rsun'])
        RV          =   filehdf5[montecarlokey]['RV'][0]        *u.Unit(filehdf5.attrs['RV'])
        mualpha     =   filehdf5[montecarlokey]['mualpha'][0]   *u.Unit(filehdf5.attrs['mualpha'])
        mu_delta    =   filehdf5[montecarlokey]['mu_delta'][0]  *u.Unit(filehdf5.attrs['mu_delta'])
    return  Mass,rh_m,RA,DEC,Rsun,RV,mualpha,mu_delta

def extract_all_GC_observables(GCnames,montecarlokey):
    assert(isinstance(GCnames,list))
    assert(isinstance(montecarlokey,str))
    Masses,rh_mes,RAes,DECes,Rsunes,RVes,mualphaes,mu_deltaes=[],[],[],[],[],[],[],[]
    for i in range(len(GCnames)):
        filename=ph.MonteCarloObservables(GCnames[i])
        

        Mass,rh_m,RA,DEC,Rsun,RV,mualpha,mu_delta=extract_GC_observables(filename,montecarlokey)
        Masses.append(Mass)
        rh_mes.append(rh_m)
        RAes.append(RA)
        DECes.append(DEC)
        Rsunes.append(Rsun)
        RVes.append(RV)
        mualphaes.append(mualpha)
        mu_deltaes.append(mu_delta)
        

    return Masses,rh_mes,RAes,DECes,Rsunes,RVes,mualphaes,mu_deltaes