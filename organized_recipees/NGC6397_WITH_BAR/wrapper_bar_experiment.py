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


# the bar mass

# the bar axis ratio 

# the bar angle

# the bar pattern speed 
# the bar pattern speeds
omega_min, omega_max = 25,66
omega_step = 0.25
bar_pattern_speeds = np.arange(omega_min,omega_max+omega_step,omega_step)


# The parameters of of the code 
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

def wrapper(bar_pattern_speed_index, bar_angle_index, bar_mass_index, bar_length_index, montecarloindex, NP):
    """
        Intended through parallelization with a slrum and the slurm job array
    """

    assert isinstance(bar_pattern_speed_index, int), "bar_pattern_speed_index must be an integer"
    assert isinstance(monte_carlo_index, int), "monte_carlo_index must be an integer"
    npatternspeeds=len(bar_pattern_speeds)
    assert bar_pattern_speed_index < npatternspeeds, "bar_pattern_speed_index must be less than the number of pattern speeds,\n {:d} was given and needs to be less than {:d}".format(bar_pattern_speed_index,npatternspeeds)
    
    montecarlokey = "monte-carlo-"+str(monte_carlo_index).zfill(3)
    ###### GET THE INITIAL CONDITIONS OF THE CLUSTER IN PHASE SPACE ######
    # get the name of the master file of initial conditions
    initial_conditions_file_name = gcs.path_handler.MonteCarloObservables(GCname)
    RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=gcs.extractors.MonteCarloObservables.extract_GC_observables(initial_conditions_file_name,monte_carlo_index)
    cluster_initial_conditions=RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m

    ###### MAKE THE PLUMMER SPHERE ######
    G=tstrippy.Parsers.potential_parameters.G
    xp,yp,zp,vxp,vyp,vzp = tstrippy.ergodic.isotropicplummer(G,Mass,rh_m,NP)
    particle_initial_conditions = np.array([xp,yp,zp,vxp,vyp,vzp])

    print("Plummer sphere with {:d} particles created".format(NP))



    barpoly = [barpoly_ferrone_2023[0], bar_pattern_speeds[bar_pattern_speed_index]]
    temp_base_name = "constant_cluster_initial_conditions_bar_pattern_speed_{:d}_m_kpc_s".format(int(1000*bar_pattern_speeds[bar_pattern_speed_index]))
    
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
        writestream=writestream)


    endtime = datetime.datetime.now()
    print("Elapsed time: ", endtime-starttime)
    print("Done with monte carlo index {:d}".format(monte_carlo_index))
    
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

    barparams           =   [22968000000, 4, 1, 1]
    barpoly_ferrone_2023=   [0.4363323129985824, 38]    
    if len (sys.argv)==4:
        scriptname = sys.argv[0]
        montecarloindex=int(sys.argv[1])
        patternspeedindex=int(sys.argv[2])
        NP = int(sys.argv[3])
        print("Running with script name {:s}".format(scriptname))
        print("Running with monte carlo index {:d}".format(montecarloindex))
        print("Running with pattern speed index {:d}".format(patternspeedindex))
        print("Running with NP {:d}".format(NP))
        single_pattern_speed(montecarloindex,patternspeedindex,NP)
    else:
        print ("NOT ENOUGH ARGUMENTS GIVEN, ")
        for arg in sys.argv:
            print(arg)

        # if DOMULTIPROCESSING:
        #     multiprocessingloop(montecarloindex)
        # else:
        #     simpleloop(montecarloindex)