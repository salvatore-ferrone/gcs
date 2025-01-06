import h5py
import gcs.path_handler as ph
from astropy import units as u


def extract_GC_observables(filename,montecarloindex):
    assert (isinstance(filename,str)), "filename must be a string"
    assert (isinstance(montecarloindex,int)), "montecarloindex must be an integer"
    with h5py.File(filename,'r') as filehdf5:
        RA,DEC,rh_m,Mass,Rsun,RV,mualpha,mu_delta=filehdf5['initial_conditions'][montecarloindex]
    return  RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m

def extract_all_GC_observables(GCnames,montecarlokey):
    assert(isinstance(GCnames,list))
    assert(isinstance(montecarlokey,str))
    Masses,rh_mes,RAes,DECes,Rsunes,RVes,mualphaes,mu_deltaes=[],[],[],[],[],[],[],[]
    for i in range(len(GCnames)):
        filename=ph.MonteCarloObservables(GCnames[i])
        

        RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=extract_GC_observables(filename,montecarlokey)
        Masses.append(Mass)
        rh_mes.append(rh_m)
        RAes.append(RA)
        DECes.append(DEC)
        Rsunes.append(Rsun)
        RVes.append(RV)
        mualphaes.append(mualpha)
        mu_deltaes.append(mu_delta)
        

    return RAes,DECes,Rsunes,RVes,mualphaes,mu_deltaes,Masses,rh_mes