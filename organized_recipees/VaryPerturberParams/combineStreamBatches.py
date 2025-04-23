"""
Since the code does batches of particles,
we can stack them together into one data file per simulation
This script does that 
"""

import gcs
import os
import sys
import numpy as np
import varyPerturberParams as vpp
import arguments as arguments
import h5py
import datetime
import multiprocessing as mp 

GLOBAL_START_TIME = datetime.datetime.now()

def extract_individual_data(fnames,NPs):
    """ 
        assumes all fnames are valid files of the same format,
        these are the data fields that vary accross the files
    """
    # set the indicies
    cummulative_NPs = np.cumsum(NPs)
    cummulative_NPs = np.insert(cummulative_NPs, 0, 0)
    total_particles = np.sum(NPs)
    # the data arrays
    phase_space = np.zeros((6,total_particles))
    tesc=np.zeros(total_particles)
    # the attributes
    computation_time = [ ]
    dates = [ ]
    for i in range(len(fnames)):
        with h5py.File(fnames[i],"r") as f:
            phase_space[:,cummulative_NPs[i]:cummulative_NPs[i+1]] = f["phase_space"][:]
            tesc[cummulative_NPs[i]:cummulative_NPs[i+1]] = f['tesc'][:]
            # extract the attributes
            computation_time.append(f.attrs['computation_time'])
            dates.append(f.attrs['date'][0])
    return phase_space, tesc, computation_time, dates

def get_stream_file_names_for_each_NP(GCname, NPs, 
                          stream_potential, internal_dynamics, montecarlokey, 
                          HostMass, HostHalfMassRadius,
                          PerturberName,PerturberMass, PerturberRadius):
    """
    Get the stream file names for each NP.
    Only grab the ones who computed successfully.
    """
    valid_nps = []
    fnames = []
    for ii in range(len(NPs)):
        filename=gcs.path_handler.StreamMassRadiusVaryPerturber(GCname=GCname,
                                        NP=int(NPs[ii]),
                                        potential_env=stream_potential,
                                        internal_dynamics=internal_dynamics,
                                        montecarlokey=montecarlokey,
                                        HostMass=int(HostMass),
                                        HostRadius=int(1000*HostHalfMassRadius),
                                        PerturberName=PerturberName,
                                        PerturberMass=int(PerturberMass),
                                        PerturberRadius=int(PerturberRadius))
        if os.path.exists(filename):
            fnames.append(filename)
            valid_nps.append(int(NPs[ii]))
        else:
            print("file does not exist",filename)
    return fnames, valid_nps



def main(montecarloindex,hostRadiusIndex,hostMassIndex,perturberIndex):


    # the data parameters that we are iterating over for stacking 
    montecarlokey = "monte-carlo-"+str(montecarloindex).zfill(3)
    PerturberMass,PerturberRadius=int(arguments.mass_radius[perturberIndex,0]),int(arguments.mass_radius[perturberIndex,1])
    HostMass=vpp.MASS_GRID[hostMassIndex]
    HostHalfMassRadius=vpp.RADIUS_GRID[hostRadiusIndex]
    NPs = [int(NP) for NP in arguments.NPs]
    # obtain the file names 
    fnames, valid_nps = get_stream_file_names_for_each_NP(
        GCname=vpp.GCname,
        NPs=NPs,
        stream_potential=vpp.stream_potential,
        internal_dynamics=vpp.internal_dynamics,
        montecarlokey=montecarlokey,
        HostMass=HostMass,
        HostHalfMassRadius=HostHalfMassRadius,
        PerturberName=vpp.VARIABLE_PERTURBER[0],
        PerturberMass=PerturberMass,
        PerturberRadius=PerturberRadius
    ) 

    # prepare the attributes for the outfile by copying the first file
    with h5py.File(fnames[0],"r") as f:
        # extract the attributes
        attributes=dict(f.attrs)
    
    # stack the phase space 
    phase_space, tesc, computation_time, date = extract_individual_data(fnames,NPs=valid_nps)

    # update the attributes
    attributes['NP'] = np.sum(valid_nps)

    # create the output file name
    outfilename=gcs.path_handler.StreamMassRadiusVaryPerturber(GCname=vpp.GCname,
                                        NP=int(np.sum(valid_nps)),
                                        potential_env=vpp.stream_potential,
                                        internal_dynamics=vpp.internal_dynamics,
                                        montecarlokey=montecarlokey,
                                        HostMass=int(HostMass),
                                        HostRadius=int(HostHalfMassRadius),
                                        PerturberName=vpp.VARIABLE_PERTURBER[0],
                                        PerturberMass=int(PerturberMass),
                                        PerturberRadius=int(PerturberRadius))
    
    # create the output file
    with h5py.File(outfilename,"w") as f:
        # create the datasets
        f.create_dataset("phase_space",data=phase_space)
        f.create_dataset("tesc",data=tesc)
        # add the attributes
        for key in attributes.keys():
            f.attrs[key] = attributes[key]

        # add the stacked information 
        f.create_group("stacked")
        f["stacked"].create_dataset("computation_time",data=computation_time)
        f["stacked"].create_dataset("date",data=date)
        f["stacked"].create_dataset("NPs",data=valid_nps)
        f['stacked'].create_dataset("fnames",data=fnames)

        # update the date and  computation time 
        # this is the time it took to run the script
        # this is not the itme it took to run the simulation
        dt=datetime.datetime.now()-GLOBAL_START_TIME
        attributes['computation_time'] = dt.total_seconds()
        attributes["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print("DONE: ",outfilename)

if __name__ == "__main__":
    # Get the command line arguments
    # if len(sys.argv) != 5:
        # print("Usage: python combineStreamBatches.py <montecarloindex> <hostRadiusIndex> <hostMassIndex> <perturberIndex>")
        # sys.exit(1)

    montecarloindex = 9
    # hostRadiusIndex = int(sys.argv[2])
    # hostMassIndex = int(sys.argv[3])
    # perturberIndex = int(sys.argv[4])
    

    N_perturber_configs = arguments.mass_radius.shape[0]
    n_host_radius=vpp.RADIUS_GRID.shape[0]
    n_host_mass=vpp.MASS_GRID.shape[0]

    args = []
    for i in range(n_host_radius):
        for j in range(n_host_mass):
            for k in range(N_perturber_configs):
                args.append((montecarloindex, i, j, k))

    pool = mp.Pool(processes=mp.cpu_count())
    pool.starmap(main, args)
    pool.close()
    pool.join()
    # Uncomment the line below to run the main function directly
    # main(montecarloindex=montecarloindex, hostRadiusIndex=hostRadiusIndex, hostMassIndex=hostMassIndex, perturberIndex=perturberIndex)