import gcs
from gcs import path_handler as ph
from astropy import units as u
import datetime
import sys 
import os
import tstrippy
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



if __name__ == "__main__" : 
    # montecarloindex = int(sys.argv[1])
    montecarloindex = 1

    MASS                =   1e5*u.Msun
    GCname              =   "Pal5"
    montecarlokey       =   "monte-carlo-"+str(montecarloindex).zfill(3)
    internal_dynamics   =   "isotropic-plummer_mass_increase"
    GCorbits_potential  =   "pouliasis2017pii-GCNBody"
    MWpotential         =   "pouliasis2017pii"
    NP                  =   int(1e0)
    T0                  =   -5e9*u.yr
    integrationtime     =   5e6*u.yr
    dt                  =   1e4*u.yr
    NSKIP               =   int(100)
    
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
        "MASS":MASS.value
    }
    
    ##### i/o files ####
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
    _,rh_m,_,_,_,_,_,_=gcs.extractors.MonteCarloObservables.extract_all_GC_observables([GCname],montecarloindex)
    xp,yp,zp,vxp,vyp,vzp = tstrippy.ergodic.isotropicplummer(G,MASS.value,rh_m[0],NP)
    rplummer= gcs.misc.half_mass_to_plummer(rh_m[0])

    # Extract the orbit  
    orbit_file_name                             =   ph.GC_orbits(GCorbits_potential,GCname)
    tH,xH,yH,zH,vxH,vyH,vzH                     =   gcs.extractors.GCOrbits.extract_whole_orbit(orbit_file_name,montecarlokey)    
    xHost,yHost,zHost,vxHost,vyHost,vzHost      =   gcs.misc.interpolate_finer_grid(tsampling,tH,xH,yH,zH,vxH,vyH,vzH)
    inithostperturber = (tsampling,xHost,yHost,zHost,vxHost,vyHost,vzHost,MASS.value,rplummer)
    # place the particle positions 
    initialkinematics = (xp+xHost[0],yp+yHost[0],zp+zHost[0],vxp+vxHost[0],vyp+vyHost[0],vzp+vzHost[0])
    
    # get the perturbers 
    GCnames     = load_GCnames_except_for_the_target(GCname)
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
    # extract the data
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
    gcs.writers.Stream.stream(outfilename,stream_final,tesc,attributes)
    print(outfilename, "saved")
    ################################################
    ############## SAVE THE SNAP SHOTS #############
    ################################################
    attributes["computation_time"]  =   float(computation_time.seconds)
    snapshottimesampling            =   tsampling[::NSKIP]
    gcs.writers.Stream.StreamSnapShots(snapshotfilename,snapshottimesampling,attributes,tempdir)
    print(snapshotfilename, "saved")
    