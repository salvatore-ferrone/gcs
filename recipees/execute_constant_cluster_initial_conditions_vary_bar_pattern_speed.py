"""
execute the constant cluster initial conditions
from stream_evolution_in_barred_potential.py

"""




import stream_evolution_in_barred_potential as sebp
from astropy import units as u
import numpy as np 
import tstrippy
import sys 
import multiprocessing as mp 
import datetime

omega_min, omega_max = 25,66
omega_step = 0.25
bar_pattern_speeds = np.arange(omega_min,omega_max+omega_step,omega_step)
bar_pattern_speeds_m_kpc_s=np.array(1000*bar_pattern_speeds, dtype=int)


DOMULTIPROCESSING = False
def main(monte_carlo_index):

    # do this a handful of times 
    GCname              =   "Pal5"
    MWpotential         =   "pouliasis2017pii"
    internal_dynamics   =   "isotropic-plummer"
    barname             =   "longmuralibar"
    barparams           =   [22968000000, 4, 1, 0.5]
    barpoly_ferrone_2023=   [0.4363323129985824, 38]
    NP                  =   int(5e1)
    integrationtime     =   5e9*u.yr
    dt                  =   1e5*u.yr
    NSKIP               =   int(10)  
    description         =   "Integrating star-particles with in a globular cluster in a galaxy with a bar. The cluster initial conditions were passed as an argument, and the particle initial conditions were passed as an argument. Therefore, these results can be compared to others where the bar properties are different."
    writestream         =   False  

    GCdata = tstrippy.Parsers.baumgardtMWGCs()
    ###### GET THE INITIAL CONDITIONS OF THE CLUSTER IN PHASE SPACE ######
    means,cov=GCdata.getGCCovarianceMatrix(GCname)
    if monte_carlo_index==0:
        RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m = means
    else:   
        temp=np.random.multivariate_normal(means,cov,1)
        RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m=temp[0]
    
    cluster_initial_conditions=RA,DEC,Rsun,RV,mualpha,mu_delta,Mass,rh_m

    ###### MAKE THE PLUMMER SPHERE ######
    G=tstrippy.Parsers.potential_parameters.G
    xp,yp,zp,vxp,vyp,vzp = tstrippy.ergodic.isotropicplummer(G,Mass,rh_m,NP)
    particle_initial_conditions = np.array([xp,yp,zp,vxp,vyp,vzp])


    # iterate over the bar pattern speeds
    starttime = datetime.datetime.now()

    if DOMULTIPROCESSING:
        ncpu = mp.cpu_count()
        Pool=mp.Pool(ncpu)

        for index in range(len(bar_pattern_speeds)):
            barpoly = [barpoly_ferrone_2023[0], bar_pattern_speeds[index]]
            temp_base_name = "constant_cluster_initial_conditions_{:d}_bar_pattern_speed_{:d}_m_kpc_s".format(monte_carlo_index,bar_pattern_speeds_m_kpc_s[index])
            Pool.apply_async(sebp.constant_cluster_initial_conditions, args=(cluster_initial_conditions,particle_initial_conditions,GCname,MWpotential,internal_dynamics,barname,barparams,barpoly,integrationtime,dt,NSKIP,temp_base_name,description,writestream))
        Pool.close()
        Pool.join()
        Pool.terminate()
    else:

        for index in range(len(bar_pattern_speeds)):

            barpoly = [barpoly_ferrone_2023[0], bar_pattern_speeds[index]]
            temp_base_name = "constant_cluster_initial_conditions_{:d}_bar_pattern_speed_{:d}_m_kpc_s".format(monte_carlo_index,bar_pattern_speeds_m_kpc_s[index])

            sebp.constant_cluster_initial_conditions(
                cluster_initial_conditions=cluster_initial_conditions,
                particle_initial_conditions=particle_initial_conditions,
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
    print("Done with monte carlo index {:d}".format(monte_carlo_index))
    print("Elapsed time: ", endtime-starttime)
    
    return None
if __name__=="__main__":
    montecarloindex=sys.argv[1]
    main(int(montecarloindex))