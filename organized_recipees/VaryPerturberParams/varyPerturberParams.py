"""
This script is in response to the reviewer comments for the paper Ferrone et al. 2025.

We want to investigate if the gap would be more present by changing the mass and size of the perturbing cluster.

"""

author = "Salvatore Ferrone"
author_affiliation = "Sapienza University of Rome"
author_email = "salvatore.ferrone@uniroma1.it"
description = "Modeling a cluster by varying the mass and size of the perturbing cluster"
import gcs
from gcs import path_handler as ph
from astropy import units as u
import datetime
import sys 
import os
import tstrippy
import csv
import numpy as np 

README = "We want to know what happens if we change the Mass and Size of the Perturbing cluster"


# MASS GRID AND RADIUS GRID
SIZE_GRID = 5
N_MASS_SAMPLING = SIZE_GRID
N_RADIUS_SAMPLING = SIZE_GRID # square grid
MASS_GRID = np.logspace(4,5.2, N_MASS_SAMPLING) # in Msun
RADIUS_GRID = np.logspace(np.log10(2),np.log10(30),N_RADIUS_SAMPLING)/1000 # in kpc
DONTCOMPUTE=False


# data params 
GCname              =   "Pal5"
internal_dynamics   =   "isotropic-plummer_mass_radius_grid"
GCorbits_potential  =   "pouliasis2017pii-GCNBody"
stream_potential    =   "pouliasis2017pii-NGC7078"
MWpotential         =   "pouliasis2017pii"
T0                  =   -5e9*u.yr
integrationtime     =   5e9*u.yr
dt                  =   1e4*u.yr
NSKIP               =   int(100)
VARIABLE_PERTURBER  =   ["NGC7078"]
perturber_names     =   ["NGC5272"]
writestreamsnapshots=   False



def main(NP,HostMass,HostHalfMassRadius,PerturberMass,PerturberRadius,montecarloindex,
         GCname              =   "Pal5",
         internal_dynamics   =   "isotropic-plummer_mass_radius_grid",
         GCorbits_potential  =   "pouliasis2017pii-GCNBody",
         stream_potential    =   "pouliasis2017pii-NGC7078",
         MWpotential         =   "pouliasis2017pii",
         T0                  =   -5e9*u.yr,
         integrationtime     =   5e9*u.yr,
         dt                  =   1e4*u.yr,
         NSKIP               =   int(1000),
         VARIABLE_PERTURBER  =   ["NGC7078"],
         perturber_names     =   ["NGC5272"],
         writestreamsnapshots=   False         
         ):

    assert isinstance(NP,int), "NP must be an integer but was {:}".format(type(NP))
    assert isinstance(HostMass,float), "HostMass must be a float but was {:}".format(type(HostMass))    
    assert isinstance(PerturberMass,float), "PerturberMass must be an integer but was {:}".format(type(PerturberMass))
    assert isinstance(HostHalfMassRadius,float), "HostHalfMassRadius must be a float but was {:}".format(type(HostHalfMassRadius))
    assert isinstance(PerturberMass,float), "MASS_RADIUS must be a float but was {:}".format(type(PerturberMass))
    assert isinstance(PerturberRadius,float), "PerturberRadius must be a float but was {:}".format(type(PerturberRadius))
    assert isinstance(montecarloindex,int), "montecarloindex must be an integer but was {:}".format(type(montecarloindex))
    
    montecarlokey       =   "monte-carlo-"+str(montecarloindex).zfill(3)
    
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
        "stream_potential":stream_potential,
        "MWpotential":MWpotential,
        "PerturberMass":PerturberMass,
        "PerturberRadius":PerturberRadius,
        "HostMass":HostMass,
        "HostHalfMassRadius":HostHalfMassRadius,
        "mass_radius_grid_mass":MASS_GRID,
        "mass_radius_grid_radius":RADIUS_GRID,
        "script_name":__file__,
    }
    
    for key in attributes:
        print(key,attributes[key])
    ##### i/o files ####
    outfilename=ph.StreamMassRadiusVaryPerturber(GCname=GCname,
                                     NP=NP,
                                     potential_env=stream_potential,
                                     internal_dynamics=internal_dynamics,
                                     montecarlokey=montecarlokey,
                                     HostMass=int(HostMass),
                                     HostRadius=int(1000*HostHalfMassRadius),
                                     PerturberName=VARIABLE_PERTURBER[0],
                                     PerturberMass=int(PerturberMass),
                                     PerturberRadius=int(1000*PerturberRadius))
    # snapshotfilename = ph.StreamSnapShotsMassRadius(GCname,NP,stream_potential,internal_dynamics,montecarlokey,int(HostMass),int(1000*HostHalfMassRadius))
    snapshotfilename = ph.StreamSnapShotsMassRadiusVaryPerturber(GCname=GCname,
                                     NP=NP,
                                     potential_env=stream_potential,
                                     internal_dynamics=internal_dynamics,
                                     montecarlokey=montecarlokey,
                                     HostMass=int(HostMass),
                                     HostRadius=int(1000*HostHalfMassRadius),
                                     PerturberName=VARIABLE_PERTURBER[0],
                                     PerturberMass=int(PerturberMass),
                                     PerturberRadius=int(1000*PerturberRadius))
    cond = False
    if os.path.exists(snapshotfilename):
        print(snapshotfilename, "Already exists. \n Skipping!")
        cond = True
    if os.path.exists(outfilename):
        print(outfilename, "Already exists. \n Skipping!")
        cond = True
    if cond:
        sys.exit(0)
    tempdir=ph._TemporaryStreamSnapShotsMassRadiusGrid(stream_potential,GCname,NP,internal_dynamics,montecarlokey,MASS_INDEX,RADIUS_INDEX)

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

        xp,yp,zp,vxp,vyp,vzp = tstrippy.ergodic.isotropicplummer(G,HostMass,HostHalfMassRadius,NP)
        rplummer= gcs.misc.half_mass_to_plummer(HostHalfMassRadius)

        # Extract the orbit  
        orbit_file_name                             =   ph.GC_orbits(GCorbits_potential,GCname)
        tH,xH,yH,zH,vxH,vyH,vzH                     =   gcs.extractors.GCOrbits.extract_whole_orbit(orbit_file_name,montecarlokey)    
        xHost,yHost,zHost,vxHost,vyHost,vzHost      =   gcs.misc.interpolate_finer_grid(tsampling,tH,xH,yH,zH,vxH,vyH,vzH)
        inithostperturber = (tsampling,xHost,yHost,zHost,vxHost,vyHost,vzHost,HostMass,rplummer)
        # place the particle positions 
        initialkinematics = (xp+xHost[0],yp+yHost[0],zp+zHost[0],vxp+vxHost[0],vyp+vyHost[0],vzp+vzHost[0])
        
        # get the perturbers, only those that have close encounters
        perturber_names= VARIABLE_PERTURBER + perturber_names # concatenate lists 
        perturbers  = load_perturbers(perturber_names,GCorbits_potential,montecarloindex)
        # replace the perturber with the one we are interested in
        ts,xs,ys,zs,Masses,r_plums = perturbers
        Masses[0]= PerturberMass
        r_plums[0] = PerturberRadius
        perturbers = ts,xs,ys,zs,Masses,r_plums

        ###############################################
        ########### INITIALIZE THE INTEGRATOR #########
        ###############################################
        integrator = tstrippy.integrator
        integrator.setstaticgalaxy(*staticgalaxy)
        integrator.setintegrationparameters(*integrationparameters)
        integrator.setinitialkinematics(*initialkinematics)
        integrator.inithostperturber(*inithostperturber)
        integrator.initperturbers(*perturbers)
        if writestreamsnapshots:
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
        attributes["perturber_names"]            =  perturber_names
        attributes["computation_time"]          =   computation_time.seconds
        attributes['perturber_masses']          =   perturbers[4]
        attributes['peturber_plummer_radii']    =   perturbers[5]
        gcs.writers.Stream.stream(outfilename,stream_final,tesc,attributes)
        print(outfilename, "saved")
        ################################################
        ############## SAVE THE SNAP SHOTS #############
        ################################################
        if writestreamsnapshots:
            snapshottimesampling            =   tsampling[::NSKIP]
            gcs.writers.Stream.StreamSnapShots(snapshotfilename,snapshottimesampling,tesc,attributes,tempdir)
            print(snapshotfilename, "saved")



def load_GCnames_except_for_the_target(GCname):
    GCnames =list(tstrippy.Parsers.baumgardtMWGCs().data['Cluster'][:])
    GCnames.remove(GCname)
    return GCnames


def load_perturbers(GCnames,GCorbits_potential,montecarloindex):
    assert isinstance(GCnames,(list,np.ndarray)), "GCnames must be a list or numpy array but was {:}".format(type(GCnames))
    assert isinstance(GCorbits_potential,str), "GCorbits_potential must be a string but was {:}".format(type(GCorbits_potential))
    assert isinstance(montecarloindex,int), "montecarloindex must be an integer but was {:}".format(type(montecarloindex))
    montecarlokey = "monte-carlo-"+str(montecarloindex).zfill(3)
    ts,xs,ys,zs,_,_,_=gcs.extractors.GCOrbits.extract_orbits_from_all_GCS(GCnames,GCorbits_potential,montecarlokey)
    _,_,_,_,_,_,Masses,rh_mes=gcs.extractors.MonteCarloObservables.extract_all_GC_observablesPre2025(GCnames,montecarloindex)
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




if __name__ == "__main__" : 

    # keep the same, since for the moment I'm just launching one case 
    montecarloindex=9
    
    assert len(sys.argv) == 4, "Usage: python3 varyPerturberParams.py NP INTERNAL_DYNAMICS_INDEX PERTURBER_INDEX"

    NP = int(sys.argv[1])
    INTERNAL_DYNAMICS_INDEX = int(sys.argv[2])
    PERTURBER_INDEX = int(sys.argv[3])

    ## THE HOST MASS AND RADIUS
    nmass = 5
    nradius = 5
    assert INTERNAL_DYNAMICS_INDEX < nmass*nradius, "input must be less than 25"
    MASS_INDEX = INTERNAL_DYNAMICS_INDEX // nradius
    RADIUS_INDEX = INTERNAL_DYNAMICS_INDEX % nradius
    HostMass = MASS_GRID[MASS_INDEX]
    HostHalfMassRadius = RADIUS_GRID[RADIUS_INDEX]
    
    ### THE PERTURBER MASS AND RADIUS
    sub_halo_mass_radius=np.loadtxt("sub_halo_mass_radius.txt",dtype=float)
    assert PERTURBER_INDEX < sub_halo_mass_radius.shape[0], "input must be less than the number of masses {:d}".format(len(masses))

    PerturberMass = sub_halo_mass_radius[PERTURBER_INDEX,0]
    PerturberRadius = sub_halo_mass_radius[PERTURBER_INDEX,1]

    # CONVERT TO PC
    PerturberRadius = PerturberRadius/1000
    
    print("starting main")
    main(NP,HostMass,HostHalfMassRadius,PerturberMass,PerturberRadius,montecarloindex)
    print("Done with",MASS_INDEX,RADIUS_INDEX,NP)

