import gcs 
import tstrippy

def load_GCnames_except_for_the_target(GCname):
    GCnames =list(tstrippy.Parsers.baumgardtMWGCs().data['Cluster'][:])
    GCnames.remove(GCname)
    return GCnames


def load_arguments(GCnames,GCorbits_potential,montecarlokey):
    assert isinstance(GCnames,list)
    ts,xs,ys,zs,_,_,_=gcs.extractors.GCOrbits.extract_orbits_from_all_GCS(GCnames,GCorbits_potential,montecarlokey)
    Masses,rh_mes,_,_,_,_,_,_=gcs.extractors.MonteCarloObservables.extract_all_GC_observables(GCnames,montecarlokey)
    r_plums = [gcs.misc.half_mass_to_plummer(rh_m).value for rh_m in rh_mes]
    Masses = [Mass.value for Mass in Masses]
    args=ts,xs,ys,zs,Masses,r_plums
    return args


def initperturbers(integrator,perturberargs):
    assert len(perturberargs) == 6
    integrator.initperturbers(*perturberargs)
    return integrator



