import tstrippy
from gcs import path_handler as ph
import gcs
import numpy as np 
import astropy.units as u
import sys 
import os
import copy
import yaml


def main(
    GCname = "Pal5",
    montecarlokey = "monte-carlo-000",
    internal_dynamics = "isotropic-plummer",
    GCorbits_potential = "pouliasis2017pii",
    MWpotential = "pouliasis2017pii",
    barname = "longmurali",
    barparams = [22968000000, 4, 1, 0.5],
    barpoly = [0.4363323129985824, 38],
    NP                  =   int(1e3),
    T                   =   5e9*u.yr,
    dt                  =   1e5*u.yr,
    NSKIP               =   int(100),
):



    GCname = "Pal5"
    montecarlokey = "monte-carlo-000"
    internal_dynamics = "isotropic-plummer"
    GCorbits_potential = "pouliasis2017pii"
    MWpotential = "pouliasis2017pii"
    barfile="../data/longmurali_ferrone_2023.yaml"
    NP = int(1e3)
    T = 1e9*u.yr
    dt = 1e4*u.yr

    # CHECK IF FILE EXISTS
    # tempdir=ph._TemporaryStreamSnapShots(MWpotential,GCname,NP,internal_dynamics,montecarlokey)
    # filename=ph.Stream(GCname,NP,MWpotential,internal_dynamics,montecarlokey)    

    staticgalaxy,integrationparameters,initialkinematics,inithostperturber  = load_vanilla_arguments(GCname,montecarlokey,internal_dynamics,GCorbits_potential,MWpotential,NP,T,dt)
    staticgalaxy_bar_adjusted = adjust_disc_mass(staticgalaxy,barparams)

    ###############################################
    ########### INITIALIZE THE INTEGRATOR #########
    ###############################################
    integrator = initialize_vanilla_integrator(staticgalaxy_bar_adjusted,integrationparameters,initialkinematics,inithostperturber)
    integrator.initgalacticbar(barname,barparams,barpoly)

    ###############################################
    ########### PERFORM THE INTEGRATION ###########
    ###############################################
    phase_space,tesc = leapfrogtofinalpositions(integrator)
    integrator.deallocate()
    # ### PSUDEO CODE
    
    attributes = {
        "NP":NP,
        "T":T,
        "dt":dt,
        "NSKIP":NSKIP,
        "GCname":GCname,
        "montecarlokey":montecarlokey,
        "internal_dynamics":internal_dynamics,
        "GCorbits_potential":GCorbits_potential,
        "MWpotential":MWpotential
    }
    
    tempdir=ph._TemporaryStreamSnapShots(MWpotential,GCname,NP,internal_dynamics,montecarlokey)
    filename=ph.Stream(GCname,NP,MWpotential,internal_dynamics,montecarlokey)
    
    # if os.path.exists(filename):
    #     print(filename, "Already exists. \n Skipping!")
    #     sys.exit(0)

    # staticgalaxy,integrationparameters,initialkinematics,inithostperturber = \
    #     load_arguments(GCname,montecarlokey,internal_dynamics,GCorbits_potential,MWpotential,NP,T,dt) 


    # ###############################################
    # ########### INITIALIZE THE INTEGRATOR #########
    # ################################
    # integrator=initialize_vanilla_integrator(staticgalaxy,integrationparameters,initialkinematics,inithostperturber)
    # integrator=initialize_write_snapshots(integrator,NSKIP,GCname,tempdir)
    
    # ###############################################
    # ########### PERFORM THE INTEGRATION ###########
    # ###############################################
    # phase_space,tesc = leapfrogtofinalpositions(integrator)
    # ##### SAVE THE FINAL SNAP SHOT 
    # gcs.writers.Stream.stream(filename,phase_space,tesc,attributes)
    # print(filename, "saved")
    
    # ### ASSEMBLE THE INTERMEDIATE SNAPSHOTS INTO ONE FILE
    # T,dt,Nstep,tsampling=gcs.misc.get_time_sampling_from_years_to_integration_units(T=T,dt=dt)
    # snapshottimesampling = tsampling[::NSKIP]
    # outfilename = ph.StreamSnapShots(GCname,NP,MWpotential,internal_dynamics,montecarlokey)
    # gcs.writers.Stream.StreamSnapShots(outfilename,snapshottimesampling,attributes,tempdir)
    # print(outfilename, "saved")
        








if __name__ == "__main__" : 
    main()