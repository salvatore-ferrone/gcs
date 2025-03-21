"""
This script will simulate the entire creation of the stream from a globular cluster in a barred potential.

It will make the stream in time, so we can see the evolution of the stream in the bar.

The output files will have the host orbit, initial conditions, and the evolution of the stream at some timestamps 

"""

import tstrippy
from gcs import path_handler as ph
import gcs
import numpy as np 
import astropy.units as u
from astropy import coordinates as coord
import sys 
import os
import datetime
import h5py
import copy

author = "Salvatore Ferrone"
author_affiliation = "Sapienza University of Rome"
author_email = "salvatore.ferrone@uniroma1.it"
description = "Integrating star-particles with in a globular cluster"


### THERE WILL BE THREE MAIN EXPERIMENT FOLDERS
### 1. VARY BAR MASS
### 2. VARY BAR LENGTH
### 3. VARY BAR AXIS RATIO 
### 4. VARY PATTERN SPEED
### 5. VARY INITIAL ANGLE 

# I will keep the bar axis ratio fixed and used a prolate ellipsoid instead of a triaxial shape model the bar

outfname = "{:s}_barMass_{:d}_barLength_{:d}_barAxisRatio_{:d}_barAngle_{:d}_barPatternSpeed_{:d}.hdf5"
valid_variable_folder_names = ["vary_bar_mass","vary_bar_length","vary_bar_axis_ratio","vary_initial_angle","vary_pattern_speed"]


def make_out_dir(MWpotential,barname,GCname,NP,montecarlokey,variable_folder_name):
    return  ph.paths['simulations'] + "StreamEvolutionInBarredPotential/" + MWpotential + "/" + barname + "/" + GCname + "/" + str(NP) + "/" + montecarlokey + "/" + variable_folder_name + "/"


def main(
        cluster_initial_conditions,
        particle_initial_conditions,
        montecarlokey,  
        variable_folder_name,
        GCname              =   "NGC3201",
        MWpotential         =   "pouliasis2017pii",
        internal_dynamics   =   "isotropic-plummer",
        barname             =   "longmuralibar",
        barparams           =   [22968000000, 4, 1, 1],
        barpoly             =   [0.4363323129985824, 38],
        integrationtime     =   5e9*u.yr,
        T0                  =   -5e9*u.yr,
        dt                  =   1e5*u.yr,
        NSKIP               =   int(10),
        temp_base_name      =   "constant_cluster_initial_conditions",
        description         =   "Integrating star-particles with in a globular cluster in a galaxy with a bar. The cluster initial conditions were passed as an argument, and the particle initial conditions were passed as an argument. Therefore, these results can be compared to others where the bar properties are different.",
        writestream         =   False,
          
):
    if variable_folder_name not in valid_variable_folder_names:
        raise ValueError("variable_folder_name must be one of the following: "+str(valid_variable_folder_names))



    NP = len(particle_initial_conditions[0]) # rederive since it is an input argument 
    # CURRENT OUTPUT FILE
    outdir = make_out_dir(MWpotential,barname,GCname,NP,montecarlokey,variable_folder_name)
    
    barMass = int(barparams[0])
    barLength = int(np.floor(1000*barparams[1]))
    barAxisRatio = int(np.floor(1000*barparams[1]/barparams[2]))
    barAngle = int(np.floor(1000*barpoly[0] * (180/np.pi)))
    barPatternSpeed = int(np.floor(1000*barpoly[1]))
    outname = outfname.format(GCname,barMass,barLength,barAxisRatio,barAngle,barPatternSpeed)
    os.makedirs(outdir,exist_ok=True)
    G=tstrippy.Parsers.potential_parameters.G
    barparams_with_G = [G]+barparams
    ################################################
    #### 1. ORBIT of the HOST GLOBULAR CLUSTER #####
    ################################################
    # LOAD THE STATIC GALAXY
    staticgalaxy_params = tstrippy.Parsers.potential_parameters.pouliasis2017pii()
    staticgalaxy = (MWpotential,staticgalaxy_params)
    # adjust for the bar 
    static_galaxy_adjusted = adjust_disc_mass(staticgalaxy,barparams)    
    # set the integration time.
    unitT = u.s*(u.kpc/u.km)
    NSTEP = int(integrationtime/dt)
    timestep=dt.to(unitT).value
    init_time = T0.to(unitT).value
    integrationparameters = (init_time,timestep,NSTEP)    
    # load the position of the host cluster
    RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=cluster_initial_conditions
    x,y,z,vx,vy,vz = sky_to_galactocentric(RA,DEC,Rsun,RV,mualpha,mu_delta)
    initialkinematics = (x,y,z,vx,vy,vz) # the velocities are positive because setbackwardorbit() makes them negative intenral to tstrippy

    # now integrate the host backward 
    integrator = tstrippy.integrator
    integrator.setstaticgalaxy(*static_galaxy_adjusted)
    integrator.setintegrationparameters(*integrationparameters)
    integrator.setinitialkinematics(*initialkinematics)
    print("barparams_with_G",barparams_with_G)
    print("barpoly",barpoly)
    integrator.initgalacticbar(barname,barparams_with_G,barpoly)
    integrator.setbackwardorbit()
    xt,yt,zt,vxt,vyt,vzt = integrator.leapfrogintime(NSTEP,1)
    timesteps = integrator.timestamps.copy()
    # flip and change the sign of the velocities
    xt=xt[0,::-1]
    yt=yt[0,::-1]
    zt=zt[0,::-1]
    vxt=-vxt[0,::-1]
    vyt=-vyt[0,::-1]
    vzt=-vzt[0,::-1]
    timesteps=timesteps[::-1]
    integrator.deallocate()



    ################################################
    #### 2. ORBIT OF THE STREAM PARTICLES ##########
    ################################################

    # make a plummer sphere
    xp,yp,zp,vxp,vyp,vzp = particle_initial_conditions
    # get plummer radius
    rplummer = gcs.misc.half_mass_to_plummer(rh_m)
    # move about the host
    xp,yp,zp=xp+xt[0],yp+yt[0],zp+zt[0]
    vxp,vyp,vzp=vxp+vxt[0],vyp+vyt[0],vzp+vzt[0]

    # adjust the integration time to set the initial time to the last time of the host
    integrationparameters = (timesteps[0],timestep,NSTEP)
    
    # perpare the integrator
    integrator = tstrippy.integrator
    integrator.setstaticgalaxy(*static_galaxy_adjusted)
    integrator.setintegrationparameters(*integrationparameters)
    integrator.setinitialkinematics(xp,yp,zp,vxp,vyp,vzp)
    integrator.initgalacticbar(barname,barparams_with_G,barpoly)
    integrator.inithostperturber(timesteps,xt,yt,zt,vxt,vyt,vzt,Mass,rplummer)
    # set the writestream
    temporary_dir_binary_files = "No temporary files were used"
    if writestream:
        temporary_dir_binary_files = ph.paths['temporary'] + MWpotential + '/'+ barname + '/'+ GCname + '/' + temp_base_name + '/'
        os.makedirs(temporary_dir_binary_files,exist_ok=True)
        print("Writing binary files to", temporary_dir_binary_files)
        integrator.initwritestream(NSKIP,GCname,temporary_dir_binary_files)
    # perform the integration
    start_time = datetime.datetime.now()
    integrator.leapfrogtofinalpositions()
    end_time = datetime.datetime.now()
    ellapsed_time_in_seconds = (end_time-start_time).total_seconds()
    # copy the final positions
    xf,yf,zf = integrator.xf.copy(),integrator.yf.copy(),integrator.zf.copy()
    vxf,vyf,vzf=integrator.vxf.copy(),integrator.vyf.copy(),integrator.vzf.copy()
    tesc=integrator.tesc.copy()
    stream_final = np.array([xf,yf,zf,vxf,vyf,vzf])
    integrator.deallocate()

    ################################################
    #### 3. SAVE THE DATA ##########################
    ################################################
    # store the host orbit
    hostorbit = np.array([xt,yt,zt,vxt,vyt,vzt])
    # set the attributes of the file
    initialconditions = [RA,DEC,Rsun,RV,mualpha,mu_delta]
    attributes={
        "GCname":GCname,
        "MWpotential": MWpotential,
        "staticgalaxy": static_galaxy_adjusted[1],
        "barname": barname,
        "barparams": barparams,
        "barpoly": barpoly,
        "internal-dynamics": internal_dynamics,
        "NP": NP,
        "integrationtime": integrationtime,
        "dt": dt,
        "NSKIP": NSKIP,
        "temp_base_name": temp_base_name,
        "temporary_dir_binary_files": temporary_dir_binary_files,
        "initialconditions": initialconditions,
        "ellapsed_time_in_seconds": ellapsed_time_in_seconds
    }
    attributes["author"] = author
    attributes["author_affiliation"] = author_affiliation
    attributes["author_email"] = author_email
    attributes["description"] = description
    attributes["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),

    # load in the snap shots 
    

    try: 
        with h5py.File(outdir+outname, 'w') as outfile:
            outfile.create_group("HostOrbit")
            outfile['HostOrbit'].create_dataset("timestamps",data=timesteps)
            outfile['HostOrbit'].create_dataset("orbit",data=hostorbit)
            outfile['HostOrbit'].create_dataset("Mass",data=Mass)
            outfile['HostOrbit'].create_dataset("rh_m",data=rh_m)

            # save the stream snapshots
            outfile.create_group("Stream")
            if writestream:
                snapshottimesampling = timesteps[::NSKIP]
            else:
                snapshottimesampling = np.array([timesteps[-1]])
            NSTAMPS = len(snapshottimesampling)
            outfile['Stream'].create_dataset("timestamps",data=snapshottimesampling)
            outfile["Stream"].create_dataset("NSTAMPS",data=NSTAMPS)
            outfile["Stream"].create_dataset("tesc",data=tesc)
            outfile["Stream"].create_group("StreamSnapShots")
            if writestream:
                filenames,indexes=get_temp_snapshot_filenames(temporary_dir_binary_files)
                for i in range(len(filenames)):
                    phase_space = read_fortran_stream_binary_file(temporary_dir_binary_files+filenames[i])
                    outfile["Stream"]["StreamSnapShots"].create_dataset(str(indexes[i]),data=phase_space)

            outfile["Stream"]["StreamSnapShots"].create_dataset("final",data=stream_final)
            for attr in attributes:
                outfile.attrs[attr] = attributes[attr]
        print(outdir+outname, "saved")
    except Exception as e:
        print(e)
        print("failed to write",outdir+outname)


def adjust_disc_mass(staticgalaxy_,barparameters_):
    """
    Remove some mass from the discs to compensate for the bar mass
    """

    staticgalaxy_bar_adjusted = copy.deepcopy(staticgalaxy_)
    M_disc1=staticgalaxy_[1][5]
    M_disc2=staticgalaxy_[1][8]
    Mbar = barparameters_[0]

    M_disc1_new = M_disc1 - Mbar/2
    M_disc2_new = M_disc2 - Mbar/2

    staticgalaxy_bar_adjusted[1][5] = M_disc1_new
    staticgalaxy_bar_adjusted[1][8] = M_disc2_new
    return staticgalaxy_bar_adjusted


def sky_to_galactocentric(RA,DEC,Rsun,RV,mualpha,mu_delta):
    skycoords = coord.SkyCoord(ra=RA*u.deg, dec=DEC*u.deg, distance=Rsun*u.kpc, pm_ra_cosdec=mualpha*u.mas/u.yr, pm_dec=mu_delta*u.mas/u.yr, radial_velocity=RV*u.km/u.s)
    refframe=tstrippy.Parsers.potential_parameters.MWreferenceframe()
    galcentric=skycoords.transform_to(refframe)
    x,y,z,vx,vy,vz=galcentric.cartesian.x.to(u.kpc).value,galcentric.cartesian.y.to(u.kpc).value,galcentric.cartesian.z.to(u.kpc).value,galcentric.velocity.d_x.to(u.km/u.s).value,galcentric.velocity.d_y.to(u.km/u.s).value,galcentric.velocity.d_z.to(u.km/u.s).value
    return x,y,z,vx,vy,vz



def read_fortran_stream_binary_file(filename):
    """The format of the files is first a 4-byte integer that gives the length of the size record.
    Next, is an array of single precision floats. This is the data

    Parameters
    ----------
    filename : str
        path to temporary file that stores the stream snapsho

    Returns
    -------
    phase_space: np.ndarray
        the phase space coordinates of the star-particles within the file
    """
    # Open the binary file
    with open(filename, 'rb') as file:
        # Read the length of the size record
        record_length = np.frombuffer(file.read(4), dtype=np.int32)[0]
        
        # Read the size information
        mysize = np.frombuffer(file.read(record_length), dtype=np.int32)
        
        # Skip the trailing record length of the size information
        file.read(4)
        
        # Read the length of the real numbers record
        record_length = np.frombuffer(file.read(4), dtype=np.int32)[0]
        
        # Read the real numbers as single precision floats
        temp = np.frombuffer(file.read(record_length), dtype=np.single)
        
        # Skip the trailing record length of the real numbers
        file.read(4)
        
        # Reshape the array according to the size information
        phase_space = temp.reshape(mysize)
        
        # Print the information
    return phase_space







def get_temp_snapshot_filenames(tempdir):
    # sort the by chronological order, not alphanumerical order 
    myfiles = np.array(os.listdir(tempdir))
    indexes = np.array([int(file.split("-")[1].split(".bin")[0]) for file in myfiles])
    sortdexes=np.argsort(indexes)
    orderedindexes = indexes[sortdexes]
    myfiles = myfiles[sortdexes]    
    return myfiles,orderedindexes




