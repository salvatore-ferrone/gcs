"""
This file creates some random sampling of the observed structural and kinematic parameters of the Globular Cluster.
This is so that we can reproduce the random experiments that we will perform in the future.

"""

import tstrippy
import numpy as np 
import datetime
import sys 
import os 
from gcs import path_handler as ph
import h5py


def main(n_samples = 50,seed_number = 40):
    
    
    myGCS=tstrippy.Parsers.baumgardtMWGCs()
    GCnames = myGCS.data['Cluster'][:]
    
    mcarlokeys = ["monte-carlo-"+str(i).zfill(3) for i in range(n_samples)]
    NGCs = len(GCnames)
    columns = myGCS._unitkeys
    columns.pop(-1)
    units = [str(myGCS.units[key]) for key in columns]
    
    for jj in range(NGCs):
        GCname = GCnames[jj]
        inner_loop(GCname,myGCS,n_samples,mcarlokeys,columns,units,seed_number+jj)
    
    
    
def inner_loop(GCname,myGCS,n_samples,mcarlokeys,columns,units,seed_number):
    outfilename=ph.MonteCarloObservables(GCname)
    # check if it exists
    if os.path.isfile(outfilename):
        print("{:s} already exists",outfilename,"\n Skipping")
        return
    
    means,cov=myGCS.getGCCovarianceMatrix(GCname)
    np.random.seed(seed_number)
    initial_conditions=np.random.multivariate_normal(means,cov,n_samples)
    # make the first entry the same as the means 
    initial_conditions[0]=means
    
    attrs=set_attributes(GCname,myGCS._pathtoclusterdata,seed_number)
    with h5py.File(outfilename,"w") as myfile:
        myfile.attrs.update(attrs)
        myfile.create_dataset("initial_conditions",data=initial_conditions)
        myfile.create_dataset("montecarlokeys",data=mcarlokeys)
        myfile.create_dataset("columns",data=columns)
        myfile.create_dataset("units",data=units)    
    return None
    
def set_attributes(GCname,path_to_file,seed_number):
    attrs = {
        "GCname": GCname,
        "Datasource": "https://people.smp.uq.edu.au/HolgerBaumgardt/globular",
        "Creationdate": datetime.datetime.now().strftime("%Y-%m-%d"),
        "Author": "Salvatore Ferrone",
        "Authoraffiliation": "Observatoire de Paris",
        "author-email": "salvatore.ferrone@obspm.fr",
        "Description": "Monte Carlo sampling of the observed structural and kinematic parameters of the Globular Cluster",
        "Datafile": path_to_file,
        "seed":seed_number
    } 
    return attrs 


if __name__=="__main__":
    main()  