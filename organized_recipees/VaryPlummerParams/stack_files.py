"""
Since the code does batches of particles,
we can stack them together into one data file per simulation
This script does that 
"""

import gcs
import os
import sys
import numpy as np
import varyPlummerParams as vpp
import arguments as arguments
import h5py
import datetime


GLOBAL_START_TIME = datetime.datetime.now()

def get_streamSnapShotsFileName(GCname, NP, potential_env, internal_dynamics, montecarlokey, MASS, RADIUS, mytype='physical'):
    """
    Get the filename for the StreamSnapShots file based on the GC name, Monte Carlo key, mass index, and radius index.
    
    Parameters:
    ----------
    GCname : str
        Name of the globular cluster.
    montecarlokey : str
        Key for the Monte Carlo simulation.
    MASS_INDEX : int
        Index for the mass.
    RADIUS : int
        Index for the radius.
    Returns:
    -------
    str
        The formatted filename.
    """ 
    valid_types = ['index', "physical"]
    if mytype not in valid_types:
        raise ValueError("type must be one of the following: {:s}".format(str(valid_types)))
    if mytype == 'index':
        assert isinstance(MASS, int), "MASS must be an integer but was {:}".format(type(MASS))
        assert isinstance(RADIUS, int), "RADIUS must be an integer but was {:}".format(type(RADIUS))
        path=gcs.path_handler._StreamSnapShots(GCname=GCname,NP=NP,potential_env=potential_env,internal_dynamics=internal_dynamics)
        fname = "{:s}-StreamSnapShots-{:s}_mass_{:s}_radius_{:s}.hdf5".format(GCname, montecarlokey, str(MASS).zfill(3), str(RADIUS).zfill(3))
        fname = path + fname
    elif mytype == 'physical':
        assert isinstance(MASS, int), "MASS must be a int IN SOLAR MASSES but was {:}".format(type(MASS))
        assert isinstance(RADIUS, int), "RADIUS must be a int in PARSECS but was  {:}".format(type(RADIUS))
        fname = gcs.path_handler.StreamSnapShotsMassRadius(GCname, NP, potential_env, internal_dynamics, montecarlokey, MASS, RADIUS)
    return fname


def grab_valid_fnames(GCname, NPs, potential_env, internal_dynamics, montecarlokey, MASS, RADIUS, mytype='physical'):

    """
    iterate over the NPs and check if the file exists
    Parameters
    ----------
    GCname : str
        Name of the globular cluster.
    NPs : list
        List of number of particles.
    potential_env : str
        Name of the potential environment.          
    internal_dynamics : str
        Name of the internal dynamics.
    montecarlokey : str
        Key for the Monte Carlo simulation.
    MASS : int
        Mass of the globular cluster in solar masses.
    RADIUS : int
        Half Mass Radius of the cluster in parsecs.
    mytype : str
        Type of the file name, either 'index' or 'physical'.
    Returns
    -------
    fnames : list
        List of valid file names.
    valid_NPs : list
        List of valid number of particles corresponding to the file names.
    """

    fnames = []
    valid_NPs=[]
    
    for i in range(len(NPs)):
        # get the file name
        fpath=get_streamSnapShotsFileName(GCname, NPs[i], potential_env, internal_dynamics, montecarlokey, MASS, RADIUS, mytype=mytype)
        # check if the file exists
        if os.path.exists(fpath):
            fnames.append(fpath)
            valid_NPs.append(NPs[i])
        else:
            print("file does not exist",fpath)
    valid_NPs=np.array(valid_NPs)
    return fnames, valid_NPs


def stack_phase_space(fnames,NPs,time_of_interest=0):
    """ assumes all fnames are valid files of the same format 
    """
    # set the indicies
    cummulative_NPs = np.cumsum(NPs)
    cummulative_NPs = np.insert(cummulative_NPs, 0, 0)
    # initiate the output arrays 
    phase_space = np.zeros((6,NPs.sum()))
    snapshottime=np.zeros(len(fnames))
    for i in range(len(fnames)):
        with h5py.File(fnames[i],"r") as f:
            target_index = np.argmin(np.abs(f['time_stamps'][:]-time_of_interest))
            phase_space[:,cummulative_NPs[i]:cummulative_NPs[i+1]] = f["StreamSnapShots"][str(target_index)][:]
            snapshottime[i]=f['time_stamps'][target_index]
    return phase_space, snapshottime

def stack_tesc(fnames,NPs,time_of_interest=0):
    """ assumes all fnames are valid files of the same format 
    """
    # set the indicies
    cummulative_NPs = np.cumsum(NPs)
    cummulative_NPs = np.insert(cummulative_NPs, 0, 0)
    # initiate the output arrays 
    tesc = np.zeros((NPs.sum()))
    snapshottime=np.zeros(len(fnames))
    for i in range(len(fnames)):
        with h5py.File(fnames[i],"r") as f:
            target_index = np.argmin(np.abs(f['time_stamps'][:]-time_of_interest))
            tesc[cummulative_NPs[i]:cummulative_NPs[i+1]] = f["tesc"][:]
            snapshottime[i]=f['time_stamps'][target_index]
    return tesc, snapshottime

def extract_and_stack_all_phase_space(outfile,time_stamps,valid_fnames,valid_NPs):
    """
    The bottle neck for sure of the code 
    """
    ntimestamps=len(time_stamps)
    for i in range(ntimestamps):
        if i%50==0:
            print("Extracting phase space for time stamp ",i," of ",len(time_stamps))
            dt=datetime.datetime.now()-GLOBAL_START_TIME
            dt=dt.total_seconds()
            print(dt," s")
        # convert the time_stamps to days
        time_of_interest= time_stamps[i]
        phase_space,_=stack_phase_space(valid_fnames,valid_NPs,time_of_interest=time_of_interest)
        outfile.create_dataset("StreamSnapShots/{:d}".format(i),data=phase_space)    


def extract_and_stack_all_phase_space_efficient(outfile, time_stamps, valid_fnames, valid_NPs):
    """
    Efficient version that loads all data into memory at once to avoid repeated file I/O.
    """
    print("Loading all data from files into memory...")
    start_time = datetime.datetime.now()
    
    # Calculate indices for placing data
    cummulative_NPs = np.cumsum(valid_NPs)
    cummulative_NPs = np.insert(cummulative_NPs, 0, 0)
    total_particles = valid_NPs.sum()
    ntimestamps = len(time_stamps)
    
    # Dictionary to store all file data: {file_index: file_data}
    all_data = {}
    file_time_indices = {}
    
    # Open all files once and load data into memory
    for i, fname in enumerate(valid_fnames):
        with h5py.File(fname, "r") as f:
            # Get the file's timestamps
            file_timestamps = f['time_stamps'][:]
            
            # For each timestamp in our target list, find the closest match in this file
            file_indices = {}
            for t_idx, t in enumerate(time_stamps):
                target_index = np.argmin(np.abs(file_timestamps - t))
                file_indices[t_idx] = target_index
            
            file_time_indices[i] = file_indices
            
            # Load all stream snapshots for this file
            file_data = {}
            for t_idx, file_t_idx in file_indices.items():
                file_data[t_idx] = f["StreamSnapShots"][str(file_t_idx)][:]
            
            all_data[i] = file_data
            dt = datetime.datetime.now() - GLOBAL_START_TIME
            dt = dt.total_seconds()
            print(dt, " s")
            print(f"Loaded data from {fname} in {(datetime.datetime.now() - start_time).total_seconds():.2f} seconds")

    
    print(f"All data loaded in {(datetime.datetime.now() - start_time).total_seconds():.2f} seconds")
    print(f"Processing {ntimestamps} timestamps for {total_particles} particles...")
    
    # Now process each timestamp and write to the output file
    for t_idx in range(ntimestamps):
        if t_idx % 50 == 0:
            dt = datetime.datetime.now() - start_time
            print(f"Processing timestamp {t_idx} of {ntimestamps} ({dt.total_seconds():.2f} s)")
        
        # Allocate the combined array for this timestamp
        combined_data = np.zeros((6, total_particles))
        
        # Fill in data from each file
        for file_idx, file_data in all_data.items():
            start_idx = cummulative_NPs[file_idx]
            end_idx = cummulative_NPs[file_idx + 1]
            combined_data[:, start_idx:end_idx] = file_data[t_idx]
        
        # Write this timestamp's data to the output file
        outfile.create_dataset(f"StreamSnapShots/{t_idx}", data=combined_data)
    
    print(f"Total processing time: {(datetime.datetime.now() - start_time).total_seconds():.2f} seconds")

def main(montecarlovalue,radiusindex,massindex):
    RADIUS=int(np.floor(1000*vpp.RADIUS_GRID[radiusindex]))
    MASS=int(vpp.MASS_GRID[massindex])
    montecarlovalue=9
    montecarlokey="monte-carlo-{:s}".format(str(montecarlovalue).zfill(3))    
    
    # extract number of particles per batch
    NPs=arguments.mysequence

    # get the valid fnames 
    valid_fnames,valid_NPs=grab_valid_fnames(
        GCname=vpp.GCname,
        NPs=NPs,
        potential_env=vpp.GCorbits_potential,
        internal_dynamics=vpp.internal_dynamics,
        montecarlokey=montecarlokey,
        MASS=MASS,
        RADIUS=RADIUS,
        mytype='physical'
    )    
    


    # extract the time_stamps and attributes from the first file
    main_keys=[]
    with h5py.File(valid_fnames[0],"r") as f:
        time_stamps=f['time_stamps'][:]
        # extract the attributes
        attributes=dict(f.attrs)
        for key in f.keys():
            print(key, f[key])
            main_keys.append(key)
    # delete the date and date from this, since it is file specific 
    del attributes['computation_time']
    del attributes['date']
    
    # update to the total number of particles
    attributes["NP"] = valid_NPs.sum()    

    # extract the dates and computation times
    n_files = len(valid_fnames)
    computation_time = np.zeros(n_files)
    dates = []
    for i in range(n_files):
        with h5py.File(valid_fnames[i],"r") as f:
            computation_time[i] = f.attrs['computation_time']
            dates.append(f.attrs['date'][0])
    attributes["computation_time"] = computation_time.sum()
    attributes["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Extract the escape time since it's constant across all the files
    tesc,_=stack_tesc(valid_fnames,valid_NPs)

    # create the new file name 
    outfilename=gcs.path_handler.StreamSnapShotsMassRadius(\
        GCname=vpp.GCname,
        NP=valid_NPs.sum(),
        potential_env=vpp.GCorbits_potential,
        internal_dynamics=vpp.internal_dynamics,
        montecarlokey=montecarlokey,
        Mass=int(MASS),
        radius=RADIUS)    
    
    # check if the file already exists
    if os.path.exists(outfilename):
        print(outfilename, "Already exists. \n Skipping!")
        return
    # create the file 
    outfile=h5py.File(outfilename,"w")
    # extract the phase space
    starttime=datetime.datetime.now()
    dt=starttime-GLOBAL_START_TIME
    dt=dt.total_seconds()
    print(dt,"Start stacking phase space ")
    # extract_and_stack_all_phase_space(outfile,time_stamps,valid_fnames,valid_NPs)
    extract_and_stack_all_phase_space_efficient(outfile,time_stamps,valid_fnames,valid_NPs)

    endtime=datetime.datetime.now()


    print("Time taken to stack phase space: ", endtime-starttime)
    # add the following attributes
    attributes["extraction_time"] = (endtime-starttime).total_seconds()
    # add the time_stamps 
    outfile.create_dataset("time_stamps",data=time_stamps)
    # add the tesc
    outfile.create_dataset("tesc",data=tesc)
    # add the stacked
    outfile.create_group("stacked")
    # add the computation time
    outfile.create_dataset("stacked/computation_time",data=computation_time)
    # add the date
    outfile.create_dataset("stacked/date",data=dates)
    # add the file names
    outfile.create_dataset("stacked/file_names",data=valid_fnames)
    # add the NPs
    outfile.create_dataset("stacked/NPs",data=valid_NPs)

    # add the attributes
    for key in attributes.keys():
        outfile.attrs[key] = attributes[key]

    # close the file
    outfile.close()
    dt=datetime.datetime.now()-GLOBAL_START_TIME
    dt=dt.total_seconds()
    print(dt," s")
    print("Done stacking phase space")
    print("Saved to ", outfilename)


    
if __name__ == "__main__":
    # get the mass and radius index from the command line
    montecarlovalue=9
    myinput = int(sys.argv[1])
    nmass = 5
    nradius = 5
    assert myinput < nmass*nradius, "input must be less than 25"
    MASS_INDEX = myinput // nradius
    RADIUS_INDEX = myinput % nradius

    # run the main function
    main(montecarlovalue,RADIUS_INDEX,MASS_INDEX)