import os
import yaml
from pathlib import Path

# Get the directory of the current script
current_script_directory = Path(__file__).parent

# Specify the path to `paths.yaml` within the gcs package directory
paths_file_path = current_script_directory / 'paths.yaml'

# Now you can use `paths_file_path` to open and read `paths.yaml`
with open(paths_file_path, 'r') as file:
    paths = yaml.safe_load(file)

        
def ParticleDistribution(GCname,Distribution_Name,NP,montecarlokey):
    assert(isinstance(GCname, str))
    assert(isinstance(NP, int))
    assert(isinstance(Distribution_Name, str))
    assert(isinstance(montecarlokey, str))
    outdir = os.path.join(_ParticleDistribution(), GCname)
    os.makedirs(outdir,exist_ok=True)
    filename = "{:s}_{:s}_{:s}_{:s}.hdf5".format(GCname, Distribution_Name, str(NP), montecarlokey)
    return os.path.join(outdir, filename)

def _ParticleDistribution():
    outpath = paths['simulations'] + "ParticleDistribution/"
    os.makedirs(outpath,exist_ok=True)
    return outpath

def MonteCarloObservables(GCname:str):
    assert(isinstance(GCname, str))
    return _MonteCarloObservables() + GCname + "-observables.hdf5"

        
def _MonteCarloObservables():
    outpath = paths['simulations'] + "MonteCarloObservables/"
    os.makedirs(outpath,exist_ok=True)
    return outpath


def _temporary_gc_orbits(MWpotential):
    mydir=paths['temporary'] + "gc_orbits/" + MWpotential
    os.makedirs(mydir,exist_ok=True)
    return mydir


def temporary_gc_orbits(MWpotential, montecarlokey, coordinate):
    assert(isinstance(MWpotential, str))
    assert(isinstance(montecarlokey, str))
    assert(isinstance(coordinate, str))
    return _temporary_gc_orbits(MWpotential) + "/" +montecarlokey + "-" + coordinate + ".npy"