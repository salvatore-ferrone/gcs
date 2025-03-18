"""
execute the constant cluster initial conditions
from stream_evolution_in_barred_potential.py

"""

import stream_evolution_in_barred_potential as sebp
from astropy import units as u
import numpy as np 
import tstrippy
import sys 
import gcs
import multiprocessing as mp 
import datetime

#########################################################
##### THE RANGE OF PARAMETERS FOR THE WE WILL PROBE #####
#########################################################
# the bar length 
BARLENGTHS = np.array([3,3.4,3.5,3.6,3.9,4,4.5,5])
n_bar_lengths = len(BARLENGTHS)
# the bar mass
bar_mass_peak = 1.6e10  # Peak value from your original selection
bar_mass_width = 8e9    # Width parameter - adjust as needed
num_points = 11         # Number of sampling points - adjust as needed
BAR_MASSES = np.linspace(bar_mass_peak-bar_mass_width, bar_mass_peak+bar_mass_width, num_points)
n_bar_masses = len(BAR_MASSES)
# the bar axis ratio 
AXIS_RATIOS = np.array([1/4,1/8,1/10])
n_axis_ratios = len(AXIS_RATIOS)
# the bar angle
BAR_ANGLES = np.arange(20,35,1)*np.pi/180
n_bar_angles = len(BAR_ANGLES)
# the bar pattern speed 
omega_min, omega_max = 25,66
omega_step = 0.25
PATTERN_SPEEDS = np.arange(omega_min,omega_max+omega_step,omega_step)
npatternspeeds  = len(PATTERN_SPEEDS)





# THE DEFAULT PARAMETERS
GCname              =   "NGC6397"
MWpotential         =   "pouliasis2017pii"
internal_dynamics   =   "isotropic-plummer"
barname             =   "longmuralibar"
integrationtime     =   5e9*u.yr
dt                  =   1e5*u.yr
NSKIP               =   int(10)  
description         =   "Integrating star-particles with in a globular cluster in a galaxy with a bar. The cluster initial conditions were passed as an argument, and the particle initial conditions were passed as an argument. Therefore, these results can be compared to others where the bar properties are different."
writestream         =   False  
DOMULTIPROCESSING   =   False

def wrapper(variable_folder_name, GCname, montecarloindex, NP, bar_angle_index, bar_pattern_speed_index, bar_mass_index, bar_length_index, bar_axis_ratio_index):
    """
        Intended through parallelization with a slrum and the slurm job array
    """

    assert isinstance(bar_pattern_speed_index, int), "bar_pattern_speed_index must be an integer"
    assert isinstance(montecarloindex, int), "monte_carlo_index must be an integer"

    assert bar_pattern_speed_index < npatternspeeds, "bar_pattern_speed_index must be less than the number of pattern speeds,\n {:d} was given and needs to be less than {:d}".format(bar_pattern_speed_index,npatternspeeds)
    assert montecarloindex < 1000, "monte_carlo_index must be less than 1000, {:d} was given".format(montecarloindex)
    assert bar_angle_index < n_bar_angles, "bar_angle_index must be less than the number of bar angles,\n {:d} was given and needs to be less than {:d}".format(bar_angle_index,n_bar_angles)
    assert bar_mass_index < n_bar_masses, "bar_mass_index must be less than the number of bar masses,\n {:d} was given and needs to be less than {:d}".format(bar_mass_index,n_bar_masses)
    assert bar_length_index < n_bar_lengths, "bar_length_index must be less than the number of bar lengths,\n {:d} was given and needs to be less than {:d}".format(bar_length_index,n_bar_lengths)
    assert bar_axis_ratio_index < n_axis_ratios, "bar_axis_ratio_index must be less than the number of bar axis ratios,\n {:d} was given and needs to be less than {:d}".format(bar_axis_ratio_index,n_axis_ratios)
    montecarlokey = "monte-carlo-"+str(montecarloindex).zfill(3)
    
    # set the bar parameters
    barmass = BAR_MASSES[bar_mass_index]
    barlength = BARLENGTHS[bar_length_index]
    baraxisratio = AXIS_RATIOS[bar_axis_ratio_index]
    barangle = BAR_ANGLES[bar_angle_index]
    barpatternspeed = PATTERN_SPEEDS[bar_pattern_speed_index]

    # shape
    barparams = [barmass,barlength,barlength*baraxisratio,barlength*baraxisratio]
    # rotation 
    barpoly = [barangle, barpatternspeed]

    temp_base_name = "{:s}-{:s}-NP-{:d}-mass-{:d}-length-{:d}-axisRatio-{:d}-angle-{:d}-patternSpeed-{:d}".format(
        GCname,montecarlokey,NP,int(barmass),int(1000*barlength),int(1000*(1/baraxisratio)),int(1000*180*barangle/np.pi),int(1000*barpatternspeed))
    ###### GET THE INITIAL CONDITIONS OF THE CLUSTER IN PHASE SPACE ######
    # get the name of the master file of initial conditions
    initial_conditions_file_name = gcs.path_handler.MonteCarloObservables(GCname)
    RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=gcs.extractors.MonteCarloObservables.extract_GC_observables(initial_conditions_file_name,montecarloindex)
    cluster_initial_conditions=RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m

    ###### MAKE THE PLUMMER SPHERE ######
    G=tstrippy.Parsers.potential_parameters.G
    xp,yp,zp,vxp,vyp,vzp = tstrippy.ergodic.isotropicplummer(G,Mass,rh_m,NP)
    particle_initial_conditions = np.array([xp,yp,zp,vxp,vyp,vzp])

    print("Plummer sphere with {:d} particles created".format(NP))


    # temp_base_name = "{:d}_m_kpc_s".format(int(1000*bar_pattern_speeds[bar_pattern_speed_index]))
    
    starttime = datetime.datetime.now()
    sebp.main(
        cluster_initial_conditions=cluster_initial_conditions,
        particle_initial_conditions=particle_initial_conditions,
        montecarlokey=montecarlokey,
        GCname=GCname,
        MWpotential=MWpotential,
        internal_dynamics=internal_dynamics,
        barname=barname,
        barparams=barparams,
        barpoly=barpoly,
        integrationtime=integrationtime,
        dt=dt,
        NSKIP=NSKIP,
        temp_base_name=temp_base_name,
        description=description,
        writestream=writestream,
        variable_folder_name=variable_folder_name)


    endtime = datetime.datetime.now()
    print("Elapsed time: ", endtime-starttime)
    print("Done with monte carlo index {:d}".format(montecarloindex))
    
    return None




# ### THESE TWO FUNCTIONS WERE USED FOR DEVELOPMENT ON HAMPI
# def multiprocessingloop(monte_carlo_index):
#     """
#         This is for multiprocessing, intended mainly for parallelization on nhampi. 

#     """

#     assert isinstance(monte_carlo_index, int), "monte_carlo_index must be an integer"
#     ###### GET THE INITIAL CONDITIONS OF THE CLUSTER IN PHASE SPACE ######
#     # get the name of the master file of initial conditions
#     initial_conditions_file_name = gcs.path_handler.MonteCarloObservables(GCname)
#     RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=gcs.extractors.MonteCarloObservables.extract_GC_observables(initial_conditions_file_name,monte_carlo_index)
#     cluster_initial_conditions=RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m

#     ###### MAKE THE PLUMMER SPHERE ######
#     G=tstrippy.Parsers.potential_parameters.G
#     xp,yp,zp,vxp,vyp,vzp = tstrippy.ergodic.isotropicplummer(G,Mass,rh_m,NP)
#     particle_initial_conditions = np.array([xp,yp,zp,vxp,vyp,vzp])


#     # iterate over the bar pattern speeds
#     starttime = datetime.datetime.now()

#     ncpu = mp.cpu_count()
#     Pool=mp.Pool(ncpu)
#     for index in range(len(bar_pattern_speeds)):
#         barpoly = [barpoly_ferrone_2023[0], bar_pattern_speeds[index]]
#         temp_base_name = "constant_cluster_initial_conditions_{:d}_bar_pattern_speed_{:d}_m_kpc_s".format(monte_carlo_index,bar_pattern_speeds_m_kpc_s[index])
#         Pool.apply_async(sebp.main, args=(cluster_initial_conditions,particle_initial_conditions,GCname,MWpotential,internal_dynamics,barname,barparams,barpoly,integrationtime,dt,NSKIP,temp_base_name,description,writestream))
#     Pool.close()
#     Pool.join()
#     Pool.terminate()

#     endtime = datetime.datetime.now()
#     print("Done with monte carlo index {:d}".format(monte_carlo_index))
#     print("Elapsed time: ", endtime-starttime)
#     return None


# def simpleloop(monte_carlo_index):
#     """
#         If you want to run this in a simple loop, without multiprocessing
#         this is looping over the different pattern speeds 

#     """

#     assert isinstance(monte_carlo_index, int), "monte_carlo_index must be an integer"
#     ###### GET THE INITIAL CONDITIONS OF THE CLUSTER IN PHASE SPACE ######
#     # get the name of the master file of initial conditions
#     initial_conditions_file_name = gcs.path_handler.MonteCarloObservables(GCname)
#     RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=gcs.extractors.MonteCarloObservables.extract_GC_observables(initial_conditions_file_name,monte_carlo_index)
#     cluster_initial_conditions=RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m

#     ###### MAKE THE PLUMMER SPHERE ######
#     G=tstrippy.Parsers.potential_parameters.G
#     xp,yp,zp,vxp,vyp,vzp = tstrippy.ergodic.isotropicplummer(G,Mass,rh_m,NP)
#     particle_initial_conditions = np.array([xp,yp,zp,vxp,vyp,vzp])


#     # iterate over the bar pattern speeds
#     starttime = datetime.datetime.now()

#     for index in range(len(bar_pattern_speeds)):
#         barpoly = [barpoly_ferrone_2023[0], bar_pattern_speeds[index]]
#         temp_base_name = "constant_cluster_initial_conditions_{:d}_bar_pattern_speed_{:d}_m_kpc_s".format(monte_carlo_index,bar_pattern_speeds_m_kpc_s[index])
#         sebp.main(
#             cluster_initial_conditions=cluster_initial_conditions,
#             particle_initial_conditions=particle_initial_conditions,
#             GCname=GCname,
#             MWpotential=MWpotential,
#             internal_dynamics=internal_dynamics,
#             barname=barname,
#             barparams=barparams,
#             barpoly=barpoly,
#             integrationtime=integrationtime,
#             dt=dt,
#             NSKIP=NSKIP,
#             temp_base_name=temp_base_name,
#             description=description,
#             writestream=writestream)

#     endtime = datetime.datetime.now()
#     print("Done with monte carlo index {:d}".format(monte_carlo_index))
#     print("Elapsed time: ", endtime-starttime)
#     return None



if __name__=="__main__":
    print(len(sys.argv), sys.argv, "arguments given")

    variable_folder_name = sys.argv[1]
    GCname = sys.argv[2]
    montecarloindex = int(sys.argv[3])
    NP = int(sys.argv[4])
    bar_angle_index = int(sys.argv[5])
    bar_pattern_speed_index = int(sys.argv[6])
    bar_mass_index = int(sys.argv[7])
    bar_length_index = int(sys.argv[8])
    bar_axis_ratio_index = int(sys.argv[9])


 
    if len (sys.argv)==10:
        scriptname = sys.argv[0]



        # print("Running {:s} with the following parameters".format(scriptname))
        # print("GCname: {:s}".format(GCname))
        # print("montecarloindex: {:d}".format(montecarloindex))
        # print("NP: {:d}".format(NP))
        # print("bar_angle_index: {:d}".format(bar_angle_index))
        # print("bar_pattern_speed_index: {:d}".format(bar_pattern_speed_index))
        # print("bar_mass_index: {:d}".format(bar_mass_index))
        # print("bar_length_index: {:d}".format(bar_length_index))
        # print("bar_axis_ratio_index: {:d}".format(bar_axis_ratio_index))

        wrapper(variable_folder_name,GCname,montecarloindex, NP, bar_angle_index, bar_pattern_speed_index, bar_mass_index, bar_length_index, bar_axis_ratio_index)
    else:
        print ("NOT ENOUGH ARGUMENTS GIVEN, ")
        for arg in sys.argv:
            print(arg)

        # if DOMULTIPROCESSING:
        #     multiprocessingloop(montecarloindex)
        # else:
        #     simpleloop(montecarloindex)