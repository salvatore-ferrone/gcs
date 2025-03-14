import tstrippy
import gcs
from gcs import path_handler as ph
from astropy import units as u
from astropy import constants as const
import datetime
import sys 
import os
import numpy as np
import h5py


def load_perturbers(GCnames,GCorbits_potential,montecarloindex):
    " get the orbits of the pertubers, masses, sizes, and also extract the file names"
    montecarlokey="monte-carlo-"+str(montecarloindex).zfill(3)
    fnamesPerturbers = []
    masses = np.zeros(len(GCnames))
    rplums = np.zeros(len(GCnames))
    for i in range(len(GCnames)):
        fname = gcs.path_handler.GC_orbits(GCorbits_potential,GCnames[i])
        fnamesPerturbers.append(fname)
        with h5py.File(fname,'r') as myfile:
            Mass=myfile[montecarlokey]['initialConditions']['Mass'][0]
            rh_m=myfile[montecarlokey]['initialConditions']['rh_m'][0]
            rplummer=gcs.misc.half_mass_to_plummer(rh_m)
            if i==0:
                ts = myfile[montecarlokey]['t'][:]
                nPertuberOrbitTimeStamps = myfile[montecarlokey]['xt'].shape[0]
                xs = np.zeros((len(GCnames),nPertuberOrbitTimeStamps))
                ys = np.zeros((len(GCnames),nPertuberOrbitTimeStamps))
                zs = np.zeros((len(GCnames),nPertuberOrbitTimeStamps))
            xs[i] = myfile[montecarlokey]['xt'][:]
            ys[i] = myfile[montecarlokey]['yt'][:]
            zs[i] = myfile[montecarlokey]['zt'][:]
            masses[i] = Mass    
            rplums[i] = rplummer
    perturbers = (ts, xs,ys,zs,masses,rplums)
    return perturbers,fnamesPerturbers
                    



if __name__ == "__main__" : 

    
    NP                  =   int(sys.argv[1])
    montecarloindex     =   0
    GCname              =   "Pal5"
    montecarlokey       =   "monte-carlo-"+str(montecarloindex).zfill(3)
    internal_dynamics   =   "isotropic-plummer"
    GCorbits_potential  =   "pouliasis2017pii-GCNBody"
    MWpotential         =   "pouliasis2017pii"
    integrationTime     =   5e9*u.yr
    T0                  =   -5e9*u.yr
    dt                  =   1e4*u.yr
    WRITESTREAM         =   True
    NSKIP               =   int(100)

    # GCnames             =   ["NGC104","NGC7078", "NGC2808", "NGC5139", "NGC5272"] # if you want to pick specific perturbers
    # LOAD ALL GLOBUALR CLUSTER NAMES 
    GCnames=tstrippy.Parsers.baumgardtMWGCs().data['Cluster']
    targedIndex = np.where(GCnames==GCname)[0][0]
    GCnames = np.delete(GCnames, targedIndex)    
    # establish the integration units 

    unitV = u.km/u.s
    unitL = u.kpc
    unitM = u.Msun
    unitT = unitL/unitV
    unitG = (unitV**2) * (unitL) / unitM
    G=const.G.to(unitG).value
    
    attributes = {
        "GCname":GCname,
        "montecarlokey":montecarlokey,
        "internal_dynamics":internal_dynamics,
        "GCorbits_potential":GCorbits_potential,
        "MWpotential":MWpotential,
        "NP":NP,
        "integrationTime":integrationTime.value,
        "T0":T0.value,
        "dt":dt.value,
        "NSKIP":NSKIP,
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
    ########### LOAD THE ARGUMENTS ################
    ###############################################


    # LOAD THE STATIC GALAXY
    staticgalaxy_params = tstrippy.Parsers.potential_parameters.pouliasis2017pii()
    staticgalaxy = (MWpotential,staticgalaxy_params)
    
    # SET THE INTEGRATION PARAMETERS
    NSTEP = int(integrationTime.value/dt.value)
    T0_ = T0.to(unitT).value # to integration units 
    dt_ = dt.to(unitT).value
    integrationTime_ = integrationTime.to(unitT).value
    integrationparameters = (T0_,dt_,NSTEP)
    tsampling = np.linspace(T0_,T0_+integrationTime_,NSTEP)
    
    # LOAD THE HOST'S ORBIT 
    orbit_file_name             =   ph.GC_orbits(GCorbits_potential,GCname)
    tH,xH,yH,zH,vxH,vyH,vzH     =   gcs.extractors.GCOrbits.extract_whole_orbit(orbit_file_name,montecarlokey)
    xHost,yHost,zHost,vxHost,vyHost,vzHost        =   gcs.misc.interpolate_finer_grid(tsampling,tH,xH,yH,zH,vxH,vyH,vzH)    

    # SAMPLE THE HOST'S PLUMMER SPHERE
    _,_,_,_,_,_,massHost,radiusHost,=gcs.extractors.MonteCarloObservables.extract_GC_observablesPre2025(
                gcs.path_handler.MonteCarloObservablesPre2025(GCname),montecarloindex,)
    radiusPlummerHost=gcs.misc.half_mass_to_plummer(radiusHost)
    xp,yp,zp,vxp,vyp,vzp=tstrippy.ergodic.isotropicplummer(G,massHost,radiusPlummerHost,NP)

    # place the plummer sphere at the host's position
    initialkinematics   = (xp+xHost[0],yp+yHost[0],zp+zHost[0],vxp+vxHost[0],vyp+vyHost[0],vzp+vzHost[0])
    # package the host's info
    inithostperturber   = (tsampling,xHost,yHost,zHost,vxHost,vyHost,vzHost,massHost,radiusPlummerHost)
    perturbers,fnamesPerturbers  = load_perturbers(GCnames,GCorbits_potential,montecarloindex)

    # update the attributes
    for i in range(len(GCnames)):
        attributes["perturber_"+GCnames[i]] = fnamesPerturbers[i]
    attributes["perturber_masses"] = perturbers[4]
    attributes["perturber_sizes"] = perturbers[5]
    attributes['T0_']=T0_
    attributes['dt_']=dt_
    attributes['NSTEP']=NSTEP


    ###############################################
    ########### INITIALIZE THE INTEGRATOR #########
    ###############################################
    tstrippy.integrator.setstaticgalaxy(*staticgalaxy)
    tstrippy.integrator.setintegrationparameters(*integrationparameters)
    tstrippy.integrator.setinitialkinematics(*initialkinematics)
    tstrippy.integrator.inithostperturber(*inithostperturber)
    tstrippy.integrator.initperturbers(*perturbers)
    # do write stream 
    if WRITESTREAM:
        tstrippy.integrator.initwritestream(NSKIP,GCname,tempdir,1000)
        print("Writing binary files to", tempdir)
    ###############################################
    ########### PERFORM THE INTEGRATION ###########
    ###############################################
    starttime=datetime.datetime.now()
    tstrippy.integrator.leapfrogtofinalpositions()
    endtime=datetime.datetime.now()
    computation_time = endtime-starttime
    print("Integration took",computation_time)
    attributes['computation_time'] = computation_time.seconds
    # DEALLOCATE THE INTEGRATOR 
    xf  = tstrippy.integrator.xf
    yf  = tstrippy.integrator.yf
    zf  = tstrippy.integrator.zf
    vxf = tstrippy.integrator.vxf
    vyf = tstrippy.integrator.vyf
    vzf = tstrippy.integrator.vzf
    tesc= tstrippy.integrator.tesc
    phase_space = np.array([xf,yf,zf,vxf,vyf,vzf])
    tstrippy.integrator.deallocate()    
    

    ################################################
    ############ THE SAVIOR OF THE DATA ############
    ################################################
    attributes["GCnames"]=GCnames
    gcs.writers.Stream.stream(outfilename,phase_space,tesc,attributes)
    print(outfilename, "saved")
    ################################################
    ############## SAVE THE SNAP SHOTS #############
    ################################################
    if WRITESTREAM:
        snapshottimesampling = tsampling[::NSKIP]
        gcs.writers.Stream.StreamSnapShots(snapshotfilename,snapshottimesampling,tesc,attributes,tempdir)
        print(snapshotfilename, "saved")
        