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

def StreamMassRadiusVaryPerturber(GCname,NP,potential_env,internal_dynamics,montecarlokey,HostMass,HostRadius,PerturberName,PerturberMass,PerturberRadius):
    """
    An experiment for varying one single perturber and the mass and radius of the host cluster
    """
    assert isinstance(HostMass,int), "HostMass should be an integer but was {:s}".format(type(HostMass))
    assert isinstance(HostRadius,int), "HostRadius should be an integer but was {:s}".format(type(HostRadius))
    assert isinstance(PerturberMass,int), "PerturberMass should be an integer but was {:s}".format(type(PerturberMass))
    assert isinstance(PerturberRadius,int), "PerturberRadius should be an integer but was {:s}".format(type(PerturberRadius))
    assert isinstance(PerturberName,str), "PerturberName should be a string but was {:s}".format(type(PerturberName))
    assert isinstance(GCname,str), "GCname should be a string but was {:s}".format(type(GCname))
    assert isinstance(NP,int), "NP should be an integer but was {:s}".format(type(NP))
    assert isinstance(potential_env,str), "potential_env should be a string but was {:s}".format(type(potential_env))
    assert isinstance(internal_dynamics,str), "internal_dynamics should be a string but was {:s}".format(type(internal_dynamics))
    assert isinstance(montecarlokey,str), "montecarlokey should be a string but was {:s}".format(type(montecarlokey))
    # RADIUS IN PC and MASS IN MSUN
    outname = GCname + "-Stream-" + montecarlokey + "_hostMass_{:s}_hostRadius_{:s}_{:s}_perturberMass_{:s}_perturberRadius_{:s}.hdf5".format(str(HostMass).zfill(3),str(HostRadius).zfill(3),PerturberName,str(PerturberMass).zfill(3),str(PerturberRadius).zfill(3))
    return _Stream(GCname,NP,potential_env,internal_dynamics) + outname

def StreamSnapShotsMassRadiusVaryPerturber(GCname,NP,potential_env,internal_dynamics,montecarlokey,HostMass,HostRadius,PerturberName,PerturberMass,PerturberRadius):
    """
    An experiment for varying one single perturber and the mass and radius of the host cluster
    """
    assert isinstance(HostMass,int), "HostMass should be an integer but was {:s}".format(type(HostMass))
    assert isinstance(HostRadius,int), "HostRadius should be an integer but was {:s}".format(type(HostRadius))
    assert isinstance(PerturberMass,int), "PerturberMass should be an integer but was {:s}".format(type(PerturberMass))
    assert isinstance(PerturberRadius,int), "PerturberRadius should be an integer but was {:s}".format(type(PerturberRadius))
    assert isinstance(PerturberName,str), "PerturberName should be a string but was {:s}".format(type(PerturberName))
    assert isinstance(GCname,str), "GCname should be a string but was {:s}".format(type(GCname))
    assert isinstance(NP,int), "NP should be an integer but was {:s}".format(type(NP))
    assert isinstance(potential_env,str), "potential_env should be a string but was {:s}".format(type(potential_env))
    assert isinstance(internal_dynamics,str), "internal_dynamics should be a string but was {:s}".format(type(internal_dynamics))
    assert isinstance(montecarlokey,str), "montecarlokey should be a string but was {:s}".format(type(montecarlokey))
    # RADIUS IN PC and MASS IN MSUN
    outname = GCname + "-StreamSnapShots-" + montecarlokey + "_hostMass_{:s}_hostRadius_{:s}_{:s}_perturberMass_{:s}_perturberRadius_{:s}.hdf5".format(str(HostMass).zfill(3),str(HostRadius).zfill(3),PerturberName,str(PerturberMass).zfill(3),str(PerturberRadius).zfill(3))
    return _StreamSnapShots(GCname,NP,potential_env,internal_dynamics) + outname


def StreamMassRadius(GCname,NP,potential_env,internal_dynamics,montecarlokey,Mass,radius):
    """
    An experiment for varying the mass and radius of the host cluster
    """
    assert isinstance(Mass,int)
    assert isinstance(radius,int)
    # RADIUS IN PC and MASS IN MSUN
    outname = GCname + "-Stream-" + montecarlokey + "_mass_{:s}_radius_{:s}.hdf5".format(str(Mass).zfill(3),str(radius).zfill(3))
    return _Stream(GCname,NP,potential_env,internal_dynamics) + outname

def StreamSnapShotsMassRadius(GCname,NP,potential_env,internal_dynamics,montecarlokey,Mass,radius):
    assert isinstance(Mass,int)
    assert isinstance(radius,int)
    # RADIUS IN PC and MASS IN MSUN
    outname = GCname + "-StreamSnapShots-" + montecarlokey + "_mass_{:s}_radius_{:s}.hdf5".format(str(Mass).zfill(3),str(radius).zfill(3))
    return _StreamSnapShots(GCname,NP,potential_env,internal_dynamics) + outname

def Orbits_DarkMatterSubhalo(PopulationName,MWpotential):
    outname = PopulationName + "-orbits.hdf5"
    dirname = _Orbits(MWpotential)
    return dirname + outname

def Stream(GCname,NP,potential_env,internal_dynamics,montecarlokey):
    outname = GCname + "-Stream-" + montecarlokey + ".hdf5"
    return _Stream(GCname,NP,potential_env,internal_dynamics) + outname


def StreamSnapShots(GCname,NP,potential_env,internal_dynamics,montecarlokey):
    outname = GCname + "-StreamSnapShots-" + montecarlokey + ".hdf5"
    return _StreamSnapShots(GCname,NP,potential_env,internal_dynamics) + outname


def old_streams(MWpotential,GCname,montecarlokey,NP):
    """This is before the new data storage system was introduced as of july 2024. It is kept for compatibility reasons.

    Parameters
    ----------
    MWpotential : str
        name of the MW potential
    GCname : STR
        name of the globular cluster
    montecarlokey : str
        the index of the sampling of the initial
    NP : integer
        number of particles

    Returns
    -------
    str
        path to data file
    """
    fname = GCname + "-stream-" + montecarlokey + ".hdf5"
    return paths['simulations'] + "Streams/" + MWpotential + "/" + GCname + "/" + str(NP) + "/" + fname


def BouldriniDarkmatterSubloes():
    """
        for now there is just one 
    """
    return "/scratch2/sferrone/simulations/DarkMatterSubhaloes/DMsubhaloMW_data.hdf5"

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


def MonteCarloObservablesPre2025(GCname:str):
    assert(isinstance(GCname, str))
    return _MonteCarloObservablesPre2025() + GCname + "-observables.hdf5"

def MonteCarloObservables(GCname:str):
    assert(isinstance(GCname, str))
    return _MonteCarloObservables() + GCname + "-observables.hdf5"


def ParticleInitialConditions(GCname:str):
    assert(isinstance(GCname, str))
    return _ParticleInitialConditions(GCname) + GCname + "-ParticleInitialConditions.hdf5"

def ForceOnOrbit(GCname:str,MWpotential:str,montecarlokey:str):
    assert(isinstance(GCname, str))
    assert(isinstance(MWpotential, str))
    assert(isinstance(montecarlokey, str))
    return _ForceOnOrbit(GCname,MWpotential) + GCname + "-" + montecarlokey + ".hdf5"

def tauDensityMaps(GCname,MWpotential,montecarlokey,NP,internal_dynamics):
    outfname = GCname+"-"+montecarlokey+"-tauDensityMap.h5"
    return _tauDensityMaps(GCname,MWpotential,NP,internal_dynamics) + outfname

def PerturberSuspects(GCname,MWpotential,montecarlokey):
    outfname = GCname + "-" + montecarlokey + "-suspects.csv"
    return _PerturberSuspects(GCname,MWpotential) + outfname

def ImpactGeometry(GCname,MWpotential,montecarlokey,suspect,targetnumber):
    """ added target number because one GC can impact a stream more than once """
    outfname = GCname + "-" + montecarlokey + "-" + suspect + "-" + str(targetnumber) + ".hdf5"
    return _ImpactGeometry(GCname,MWpotential,montecarlokey) + outfname


def _Orbits(MWpotential):
    outdir = paths['simulations'] + "Orbits/" + MWpotential + "/"
    os.makedirs(outdir,exist_ok=True)
    return outdir

def _ImpactGeometry(GCname,MWpotential,montecarlokey):
    outdir = paths['simulations'] + "ImpactGeometry/" + MWpotential + "/" + GCname + "/" + montecarlokey + "/"
    os.makedirs(outdir,exist_ok=True)
    return outdir


def _PerturberSuspects(GCname,MWpotential):
    outdir = paths['simulations'] + "PerturberSuspects/" + MWpotential + "/" + GCname + "/" 
    os.makedirs(outdir,exist_ok=True)
    return outdir 
    

def _tauDensityMaps(GCname,MWpotential,NP,internal_dynamics):
    outdir=paths['simulations'] + "tauDensityMaps/" + MWpotential + "/" + GCname + "/" + str(NP) + "/" + internal_dynamics+"/"
    os.makedirs(outdir,exist_ok=True)
    return outdir

def _ForceOnOrbit(GCname:str,MWpotential:str):
    assert(isinstance(GCname, str))
    assert(isinstance(MWpotential, str))
    outpath = paths['simulations'] + "ForceOnOrbit/" + MWpotential + "/" + GCname + "/"
    os.makedirs(outpath,exist_ok=True)
    return outpath

def _Stream(GCname,NP,potential_env,internal_dynamics):
    outpath = paths['simulations'] + "Stream/" + potential_env + "/" + GCname + "/" + str(NP) + "/" +internal_dynamics + "/"
    os.makedirs(outpath,exist_ok=True)
    return outpath

def _StreamSnapShots(GCname,NP,potential_env,internal_dynamics):
    outpath = paths['simulations'] + "StreamSnapShots/" + potential_env + "/" + GCname + "/" + str(NP) + "/" +internal_dynamics + "/"
    os.makedirs(outpath,exist_ok=True)
    return outpath

def _TemporaryStreamSnapShotsMassRadiusGrid(MWpotential,GCname,NP,internal_dynamics,montecarlokey,MassIndex,radiusIndex):
    assert(isinstance(MassIndex,int))
    assert(isinstance(radiusIndex,int))
    mass_radius = "mass_{:s}_radius_{:s}".format(str(MassIndex).zfill(3),str(radiusIndex).zfill(3))
    outpath = paths['temporary'] + "StreamSnapShots/" + MWpotential + "/" + GCname + "/" + str(NP) + "/" + internal_dynamics + "/" + montecarlokey + "/" + mass_radius + "/"
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

def _MonteCarloObservablesPre2025():
    outpath = paths['simulations'] + "MonteCarloObservablesPre2025/"
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