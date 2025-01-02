import gcs
from gcs import path_handler as ph
import GCNBody # type: ignore
import vanilla
from astropy import units as u
import datetime


if __name__ == "__main__" : 
    import sys 
    import os
    GCname              =   "Pal5"
    montecarlokey       =   "monte-carlo-"+str(int(sys.argv[1])).zfill(3)
    internal_dynamics   =   "isotropic-plummer"
    GCorbits_potential  =   "pouliasis2017pii-GCNBody"
    MWpotential         =   "pouliasis2017pii"
    NP                  =   int(1e5)
    T                   =   5e9*u.yr
    dt                  =   1e4*u.yr
    NSKIP               =   int(100)
    
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
    
    outfilename=ph.Stream(GCname,NP,GCorbits_potential,internal_dynamics,montecarlokey)
    snapshotfilename = ph.StreamSnapShots(GCname,NP,GCorbits_potential,internal_dynamics,montecarlokey)
    cond = False
    if os.path.exists(snapshotfilename):
        print(snapshotfilename, "Already exists. \n Skipping!")
        cond = True
    if os.path.exists(outfilename):
        print(outfilename, "Already exists. \n Skipping!")
        cond = True
    if cond:
        sys.exit(0)

    tempdir=ph._TemporaryStreamSnapShots(GCorbits_potential,GCname,NP,internal_dynamics,montecarlokey)
    ###############################################
    ########### LOAD THE ARGUMENTS ###############
    ###############################################
    staticgalaxy,integrationparameters,initialkinematics,inithostperturber = \
        vanilla.load_arguments(GCname,montecarlokey,internal_dynamics,GCorbits_potential,MWpotential,NP,T,dt) 
    GCnames     = GCNBody.load_GCnames_except_for_the_target(GCname)
    perturbers  = GCNBody.load_perturbers(GCnames,GCorbits_potential,montecarlokey)
    ###############################################
    ########### INITIALIZE THE INTEGRATOR #########
    ###############################################
    integrator=vanilla.initialize_vanilla_integrator(staticgalaxy,integrationparameters,initialkinematics,inithostperturber)
    integrator=GCNBody.initperturbers(integrator,perturbers) 
    integrator=vanilla.initialize_write_snapshots(integrator,NSKIP,GCname,tempdir)
    print("Writing binary files to", tempdir)
    ###############################################
    ########### PERFORM THE INTEGRATION ###########
    ###############################################
    starttime=datetime.datetime.now()
    phase_space,tesc = vanilla.leapfrogtofinalpositions(integrator)
    endtime=datetime.datetime.now()
    print("Integration took",endtime-starttime)
    ################################################
    ############ THE SAVIOR OF THE DATA ############
    ################################################
    attributes["GCnames"]=GCnames
    gcs.writers.Stream.stream(outfilename,phase_space,tesc,attributes)
    print(outfilename, "saved")
    ################################################
    ############## SAVE THE SNAP SHOTS #############
    ################################################
    T,dt,Nstep,tsampling=gcs.misc.get_time_sampling_from_years_to_integration_units(T=T,dt=dt)
    snapshottimesampling = tsampling[::NSKIP]
    gcs.writers.Stream.StreamSnapShots(snapshotfilename,snapshottimesampling,attributes,tempdir)
    print(snapshotfilename, "saved")
    