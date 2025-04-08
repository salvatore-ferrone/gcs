import h5py
import gcs.path_handler as ph
from astropy import units as u
import numpy as np 



def extract_GC_observablesPre2025(filename,montecarloindex):
    assert isinstance(filename,str), "filename must be a string but instead was: "+str(type(filename))
    assert isinstance(montecarloindex,int), "montecarloindex must be an integer but instead was: "+str(type(montecarloindex))
    montecarlokey = "monte-carlo-"+str(montecarloindex).zfill(3)
    with h5py.File(filename,'r') as filehdf5:
        RA =filehdf5[montecarlokey]['RA'][0]
        DEC =filehdf5[montecarlokey]['DEC'][0]
        Rsun =filehdf5[montecarlokey]['Rsun'][0]
        RV =filehdf5[montecarlokey]['RV'][0]
        mualpha =filehdf5[montecarlokey]['mualpha'][0]
        mu_delta =filehdf5[montecarlokey]['mu_delta'][0]
        Mass =filehdf5[montecarlokey]['Mass'][0]
        rh_m =filehdf5[montecarlokey]['rh_m'][0]
    return  RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m

def extract_all_GC_observablesPre2025(GCnames,montecarloindex):
    assert isinstance(GCnames,(list,tuple,np.ndarray)), "GCnames must be a list but instead was: "+str(type(GCnames))
    assert isinstance(montecarloindex,int), "montecarlokey must be a string but instead was: "+str(type(montecarloindex))
    Masses,rh_mes,RAes,DECes,Rsunes,RVes,mualphaes,mu_deltaes=[],[],[],[],[],[],[],[]
    for i in range(len(GCnames)):
        filename=ph.MonteCarloObservablesPre2025(GCnames[i])
        RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=extract_GC_observablesPre2025(filename,montecarloindex)
        Masses.append(Mass)
        rh_mes.append(rh_m)
        RAes.append(RA)
        DECes.append(DEC)
        Rsunes.append(Rsun)
        RVes.append(RV)
        mualphaes.append(mualpha)
        mu_deltaes.append(mu_delta)
    return RAes,DECes,Rsunes,RVes,mualphaes,mu_deltaes,Masses,rh_mes

def extract_GC_observables(filename,montecarloindex):
    assert isinstance(filename,str), "filename must be a string but instead was: "+str(type(filename))
    assert isinstance(montecarloindex,int), "montecarloindex must be an integer but instead was: "+str(type(montecarloindex))
    with h5py.File(filename,'r') as filehdf5:
        RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=filehdf5['initial_conditions'][montecarloindex]
    return  RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m

def extract_all_GC_observables(GCnames,montecarloindex):
    assert isinstance(GCnames,list), "GCnames must be a list but instead was: "+str(type(GCnames))
    assert isinstance(montecarloindex,int), "montecarlokey must be a string but instead was: "+str(type(montecarloindex))
    Masses,rh_mes,RAes,DECes,Rsunes,RVes,mualphaes,mu_deltaes=[],[],[],[],[],[],[],[]
    for i in range(len(GCnames)):
        filename=ph.MonteCarloObservables(GCnames[i])
        RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=extract_GC_observables(filename,montecarloindex)
        Masses.append(Mass)
        rh_mes.append(rh_m)
        RAes.append(RA)
        DECes.append(DEC)
        Rsunes.append(Rsun)
        RVes.append(RV)
        mualphaes.append(mualpha)
        mu_deltaes.append(mu_delta)
    return RAes,DECes,Rsunes,RVes,mualphaes,mu_deltaes,Masses,rh_mes