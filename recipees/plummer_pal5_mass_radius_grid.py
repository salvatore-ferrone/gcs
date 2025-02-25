"""
This script is in response to the reviewer comments for the paper Ferrone et al. 2025.

particularly in modeling the internal dynamics of the globular cluster in the Milky Way potential.

The reviewer wants to see a sensitivity analysis to the internal dynamics and for changing the 
 mass and radius of palomar 5 

"""

author = "Salvatore Ferrone"
author_affiliation = "Sapienza University of Rome"
author_email = "salvatore.ferrone@uniroma1.it"
description = "modeling a cluster with varying mass and size "
import gcs
from gcs import path_handler as ph
from astropy import units as u
import datetime
import sys 
import os
import tstrippy
import csv
import numpy as np 
README = "The reviewer wanted a more massive Palomar 5. This script is the answer to that request. The script is a copy of the one in recipees/execute_GCNBody_Palomar5.py with the only difference being the mass of the Palomar 5 cluster. The mass of the Palomar 5 cluster is set to XXX Msun."


def load_GCnames_except_for_the_target(GCname):
    GCnames =list(tstrippy.Parsers.baumgardtMWGCs().data['Cluster'][:])
    GCnames.remove(GCname)
    return GCnames


def load_perturbers(GCnames,GCorbits_potential,montecarloindex):
    assert isinstance(GCnames,list)
    montecarlokey = "monte-carlo-"+str(montecarloindex).zfill(3)
    ts,xs,ys,zs,_,_,_=gcs.extractors.GCOrbits.extract_orbits_from_all_GCS(GCnames,GCorbits_potential,montecarlokey)
    Masses,rh_mes,_,_,_,_,_,_=gcs.extractors.MonteCarloObservables.extract_all_GC_observables(GCnames,montecarloindex)
    r_plums = [gcs.misc.half_mass_to_plummer(rh_m) for rh_m in rh_mes]
    Masses = [Mass for Mass in Masses]
    perturbers=ts,xs,ys,zs,Masses,r_plums
    return perturbers


def initperturbers(integrator,perturberargs):
    assert len(perturberargs) == 6
    integrator.initperturbers(*perturberargs)
    return integrator


def load_targetted_gcs(GCname,MWpotential,montecarlokey):
    path_suspects = gcs.path_handler.PerturberSuspects(GCname,MWpotential,montecarlokey)
    suspects=[]
    with open(path_suspects, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            suspects.append(row['suspects']) 
    return suspects



# MASS GRID AND RADIUS GRID
SIZE_GRID = 5
N_MASS_SAMPLING = SIZE_GRID
N_RADIUS_SAMPLING = SIZE_GRID # square grid
MASS_GRID = np.logspace(4,5.2,N_MASS_SAMPLING) # in Msun
RADIUS_GRID = np.linspace(5,30,N_RADIUS_SAMPLING)/1000 # in kpc
DONTCOMPUTE=False
if __name__ == "__main__" : 
    # MASS_INDEX = int(sys.argv[2])
    # RADIUS_INDEX = int(sys.argv[3])
    linear_index = int(sys.argv[1])
    montecarloindex = int(9)
    # iterate over the masses and then the radii
    MASS_INDEX = linear_index//SIZE_GRID
    RADIUS_INDEX = linear_index%SIZE_GRID

    assert len(sys.argv) == 2, "Usage: python3 plummer_pal5_mass_radius_grid.py linear_index"
    assert linear_index < N_MASS_SAMPLING*N_RADIUS_SAMPLING, "linear_index must be less than {:d}".format(N_MASS_SAMPLING*N_RADIUS_SAMPLING)
    assert MASS_INDEX < N_MASS_SAMPLING, "MASS_INDEX must be less than {:d}".format(N_MASS_SAMPLING)
    assert RADIUS_INDEX < N_RADIUS_SAMPLING, "RADIUS_INDEX must be less than {:d}".format(N_RADIUS_SAMPLING)



    MASS                =   MASS_GRID[MASS_INDEX]*u.Msun
    HALF_MASS_RADIUS    =   RADIUS_GRID[RADIUS_INDEX]*u.kpc
    GCname              =   "Pal5"
    montecarlokey       =   "monte-carlo-"+str(montecarloindex).zfill(3)
    internal_dynamics   =   "isotropic-plummer_mass_radius_grid"
    GCorbits_potential  =   "pouliasis2017pii-GCNBody"
    MWpotential         =   "pouliasis2017pii"
    NP                  =   int(5)
    T0                  =   -5e9*u.yr
    integrationtime     =   5e9*u.yr
    dt                  =   1e4*u.yr
    NSKIP               =   int(10000)
    
    text_append_filename = "_mass_{:2d}_radius_{:2d}".format(MASS_INDEX,RADIUS_INDEX)

    attributes = {
        "README":README,
        "NP":NP,
        "T":T0.value,
        "dt":dt.value,
        "NSKIP":NSKIP,
        "GCname":GCname,
        "montecarlokey":montecarlokey,
        "internal_dynamics":internal_dynamics,
        "GCorbits_potential":GCorbits_potential,
        "MWpotential":MWpotential,
        "MASS":MASS.value,
        "HALF_MASS_RADIUS":HALF_MASS_RADIUS.value
    }
    
    ##### i/o files ####

    STRINGMASS = "{:d}".format(int(MASS.value)).zfill(3)
    STRINGRADIUS = "{:d}".format(int(HALF_MASS_RADIUS.value)).zfill(3)
    outfilename=ph.StreamMassRadius(GCname,NP,GCorbits_potential,internal_dynamics,montecarlokey,MASS_INDEX,RADIUS_INDEX)
    snapshotfilename = ph.StreamShapShotsMassRadius(GCname,NP,GCorbits_potential,internal_dynamics,montecarlokey,MASS_INDEX,RADIUS_INDEX)
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

    if DONTCOMPUTE == True:
        print("DONTCOMPUTE is set to True. Exiting")
        print("outfilename",outfilename)
        print("snapshotfilename",snapshotfilename)

    else:

        #### SET TIME STEPS ####
        unitV = u.km/u.s
        unitL = u.kpc
        unitT = unitL/unitV
        NSTEP = int((integrationtime.to(u.yr).value)/dt.to(u.yr).value)
        T0=T0.to(unitT)
        integrationtime = integrationtime.to(unitT)
        dt = dt.to(unitT)
        tsampling = np.arange(T0.value,T0.value+integrationtime.value+dt.value,dt.value)


        ###############################################
        ########### LOAD THE ARGUMENTS ###############
        ###############################################

        # get the potential of the MW
        MWparams = tstrippy.Parsers.potential_parameters.pouliasis2017pii()
        staticgalaxy = (MWpotential,MWparams)

        # set the integration parameters
        integrationparameters = (T0.value,dt.value,NSTEP)
        tsampling = np.linspace(T0.value,T0.value+NSTEP*dt.value,NSTEP+1)
        # sample the new plummer distribution
        G = MWparams[0]
        ## make a new sampling of the plummer sphere for the new mass

        xp,yp,zp,vxp,vyp,vzp = tstrippy.ergodic.isotropicplummer(G,MASS.value,HALF_MASS_RADIUS.value,NP)
        rplummer= gcs.misc.half_mass_to_plummer(HALF_MASS_RADIUS.value)

        # Extract the orbit  
        orbit_file_name                             =   ph.GC_orbits(GCorbits_potential,GCname)
        tH,xH,yH,zH,vxH,vyH,vzH                     =   gcs.extractors.GCOrbits.extract_whole_orbit(orbit_file_name,montecarlokey)    
        xHost,yHost,zHost,vxHost,vyHost,vzHost      =   gcs.misc.interpolate_finer_grid(tsampling,tH,xH,yH,zH,vxH,vyH,vzH)
        inithostperturber = (tsampling,xHost,yHost,zHost,vxHost,vyHost,vzHost,MASS.value,rplummer)
        # place the particle positions 
        initialkinematics = (xp+xHost[0],yp+yHost[0],zp+zHost[0],vxp+vxHost[0],vyp+vyHost[0],vzp+vzHost[0])
        
        # get the perturbers, only those that have close encounters
        # GCnames     = load_GCnames_except_for_the_target(GCname)
        GCnames    = load_targetted_gcs(GCname,MWpotential,montecarlokey)
        perturbers  = load_perturbers(GCnames,GCorbits_potential,montecarloindex)
        
        ###############################################
        ########### INITIALIZE THE INTEGRATOR #########
        ###############################################
        integrator = tstrippy.integrator
        integrator.setstaticgalaxy(*staticgalaxy)
        integrator.setintegrationparameters(*integrationparameters)
        integrator.setinitialkinematics(*initialkinematics)
        integrator.inithostperturber(*inithostperturber)
        integrator.initperturbers(*perturbers)
        integrator.initwritestream(NSKIP,GCname,tempdir)
        print("Writing binary files to", tempdir)
        ###############################################
        ########### PERFORM THE INTEGRATION ###########
        ###############################################
        starttime=datetime.datetime.now()
        integrator.leapfrogtofinalpositions()
        ########### extract the data ##################
        xf,yf,zf            = integrator.xf.copy(),integrator.yf.copy(),integrator.zf.copy()
        vxf,vyf,vzf         = integrator.vxf.copy(),integrator.vyf.copy(),integrator.vzf.copy()
        tesc                = integrator.tesc.copy()
        stream_final        = np.array([xf,yf,zf,vxf,vyf,vzf])
        endtime             = datetime.datetime.now()
        computation_time    = endtime-starttime
        print("Integration took",endtime-starttime)
        ################################################
        ############ THE SAVIOR OF THE DATA ############
        ################################################
        attributes["GCnames"]           =   GCnames
        attributes["computation_time"]  =   (computation_time.seconds)
        gcs.writers.Stream.stream(outfilename,stream_final,tesc,attributes)
        print(outfilename, "saved")
        ################################################
        ############## SAVE THE SNAP SHOTS #############
        ################################################
        snapshottimesampling            =   tsampling[::NSKIP]
        gcs.writers.Stream.StreamSnapShots(snapshotfilename,snapshottimesampling,tesc,attributes,tempdir)
        print(snapshotfilename, "saved")
    
        