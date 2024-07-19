import tstrippy
import h5py 
import datetime 
import os 
import sys 
import multiprocessing as mp
from gcs import path_handler as ph


potential_name = "isotropic_plummer"

def main(GCname, NP, ncpu = 10):
    
    
    # get the file of the initial conditions 
    init_conds_path = ph.MonteCarloObservables(GCname)    
    with h5py.File(init_conds_path, 'r') as f:
        masses = f["initial_conditions"][:,6]
        rh_m = f["initial_conditions"][:,7]
        n_samples = len(masses)
    montecarlokeys=["monte-carlo-"+str(i).zfill(3) for i in range(n_samples)]    
    # create the pool
    pool = mp.Pool(ncpu)
    
    # create the arguments
    args = [(GCname, NP, montecarlokeys[i], masses[i], rh_m[i]) for i in range(n_samples)]
    # run the pool
    pool.starmap(inner_loop, args)
    # close the pool
    pool.close()
    pool.join()
    return None
           
def inner_loop(GCname, NP, montecarlokey, Mass, rh_m):    
    # make sure that the user entered index is valid
        
    outname=ph.ParticleDistribution(GCname, potential_name, NP, montecarlokey)
    if os.path.isfile(outname):
        print("{:s} already exists",outname,"\n Skipping")
        return
    # extract the data for making the initial conditions
    G=tstrippy.Parsers.potential_parameters.G
    attrs = set_attributes(GCname,NP,montecarlokey)
    # perform the sampling
    x,y,z,vx,vy,vz=tstrippy.ergodic.isotropicplummer(G,Mass,rh_m,NP)
    # write the file
    write_to_file(outname,x,y,z,vx,vy,vz,Mass,rh_m,attrs)
    return None
    
    

def set_attributes(GCname,NP,mcarlokey):
    return {
        "GCname": GCname,
        "Datasource": "https://people.smp.uq.edu.au/HolgerBaumgardt/globular",
        "Creationdate": datetime.datetime.now().strftime("%Y-%m-%d"),
        "Author": "Salvatore Ferrone",
        "Authoraffiliation": "Observatoire de Paris",
        "author-email": "salvatore.ferrone@obspm.fr",
        "Description": "Inverse transform sampling of the isotropic Plummer model",
        "potential_name":potential_name,
        "Nparticles": NP,
        "monte-carlo-key": mcarlokey,
    }     
    
    
def write_to_file(outname,x,y,z,vx,vy,vz,Mass,rh_m,attrs):
    with h5py.File(outname,"w") as myfile:
        myfile.create_dataset("x",data=x)
        myfile.create_dataset("y",data=y)
        myfile.create_dataset("z",data=z)
        myfile.create_dataset("vx",data=vx)
        myfile.create_dataset("vy",data=vy)
        myfile.create_dataset("vz",data=vz)
        myfile.create_dataset("Mass",data=Mass)
        myfile.create_dataset("rh_m",data=rh_m)
        myfile.attrs.update(attrs)
        
        
        
if __name__=="__main__":
    # GCname = sys.argv[1]
    # NP  = int(sys.argv[2])
    GCname = "Pal5"
    NP=int(1e5)
    main(GCname, NP)
