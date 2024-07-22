import os
import yaml
from pathlib import Path

# Get the directory of the current script
current_script_directory = Path(__file__).parent

# Specify the path to `paths.yaml` within the gcs package directory
paths_file_path = current_script_directory / 'paths.yaml'

# if not paths_file_path.exists():
#     errmessage = f"File not found: {paths_file_path}\n \
#     Create a YAML file with a place to store the simulation results and the temporary files. \n \
#     For example, the file should look like this: \n \
#         simulations: /home/sferrone/tychoOUTPUTS/
#         temporary: /home/sferrone/temp/"
    
#     raise FileNotFoundError(f"File not found: {paths_file_path}")

# Now you can use `paths_file_path` to open and read `paths.yaml`
with open(paths_file_path, 'r') as file:
    paths = yaml.safe_load(file)


def Stream(GCname,NP,potential_env,internal_dynamics,montecarlokey):
    outname = GCname + "-Stream-" + montecarlokey + ".hdf5"
    return _Stream(GCname,NP,potential_env,internal_dynamics) + outname


def StreamSnapShots(GCname,NP,potential_env,internal_dynamics,montecarlokey):
    outname = GCname + "-StreamSnapShots-" + montecarlokey + ".hdf5"
    return _StreamSnapShots(GCname,NP,potential_env,internal_dynamics) + outname




def GC_orbits(MWpotential, GCname):
    assert(isinstance(GCname, str))
    assert(isinstance(MWpotential, str))
    return _GC_orbits(MWpotential) + GCname + "-orbits.hdf5"


def ParticleDistribution(GCname,Distribution_Name,NP,montecarlokey):
    assert(isinstance(GCname, str))
    assert(isinstance(NP, int))
    assert(isinstance(Distribution_Name, str))
    assert(isinstance(montecarlokey, str))
    outdir = os.path.join(_ParticleDistribution(), GCname)
    os.makedirs(outdir,exist_ok=True)
    filename = "{:s}_{:s}_{:s}_{:s}.hdf5".format(GCname, Distribution_Name, str(NP), montecarlokey)
    return os.path.join(outdir, filename)


def MonteCarloObservables(GCname:str):
    assert(isinstance(GCname, str))
    return _MonteCarloObservables() + GCname + "-observables.hdf5"


def ParticleInitialConditions(GCname:str):
    assert(isinstance(GCname, str))
    return _ParticleInitialConditions(GCname) + GCname + "-ParticleInitialConditions.hdf5"


def _Stream(GCname,NP,potential_env,internal_dynamics):
    outpath = paths['simulations'] + "Stream/" + potential_env + "/" + GCname + "/" + str(NP) + "/" +internal_dynamics + "/"
    os.makedirs(outpath,exist_ok=True)
    return outpath

def _StreamSnapShots(GCname,NP,potential_env,internal_dynamics):
    outpath = paths['simulations'] + "StreamSnapShots/" + potential_env + "/" + GCname + "/" + str(NP) + "/" +internal_dynamics + "/"
    os.makedirs(outpath,exist_ok=True)
    return outpath

def _TemporaryStreamSnapShots(MWpotential,GCname,NP,internal_dynamics,montecarlokey):
    outpath = paths['temporary'] + "StreamSnapShots/" + MWpotential + "/" + GCname + "/" + str(NP) + "/" + internal_dynamics + "/" + montecarlokey + "/"
    os.makedirs(outpath,exist_ok=True)
    return outpath

def _ParticleInitialConditions(GCname):
    outpath = paths['simulations'] + "ParticleInitialConditions/"
    os.makedirs(outpath,exist_ok=True)
    return outpath


def _ParticleDistribution():
    outpath = paths['simulations'] + "ParticleDistribution/"
    os.makedirs(outpath,exist_ok=True)
    return outpath
        

def _MonteCarloObservables():
    outpath = paths['simulations'] + "MonteCarloObservables/"
    os.makedirs(outpath,exist_ok=True)
    return outpath


def _GC_orbits(MWpotential):
    mydir=paths['simulations'] + "Orbits/"+ MWpotential +"/"
    os.makedirs(mydir,exist_ok=True)
    return mydir


def _temporary_gc_orbits(MWpotential):
    mydir=paths['temporary'] + "gc_orbits/" + MWpotential
    os.makedirs(mydir,exist_ok=True)
    return mydir


def temporary_gc_orbits(MWpotential, montecarlokey, coordinate):
    assert(isinstance(MWpotential, str))
    assert(isinstance(montecarlokey, str))
    assert(isinstance(coordinate, str))
    return _temporary_gc_orbits(MWpotential) + "/" +montecarlokey + "-" + coordinate + ".npy"