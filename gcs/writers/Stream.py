import os
import h5py 
import datetime
import numpy as np 

author = "Salvatore Ferrone"
author_affiliation = "Sapienza University of Rome"
author_email = "salvatore.ferrone@uniroma1.it"
description = "Integrating star-particles with in a globular cluster"




def stream(filename,phase_space,tesc,attributes):
    """
    This function writes the stream to a file
    """
    assert stream.shape[1] == int(attributes["NP"])
    

    
    phase_space_columns = ["x","y","z","vx","vy","vz"]
    
    attributes["author"] = author
    attributes["author_affiliation"] = author_affiliation
    attributes["author_email"] = author_email
    attributes["description"] = description
    attributes["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    with h5py.File(filename, 'w') as f:
        dset = f.create_dataset("phase_space", data=phase_space)
        dset.attrs['phase_space_columns'] = phase_space_columns
        f.create_dataset("tesc", data=tesc)
        for attr in attributes:
            f.attrs[attr] = attributes[attr]
    return None


def get_temp_snapshot_filenames(tempdir):
    # sort the files 
    myfiles = np.array(os.listdir(tempdir))
    indexes = np.array([int(file.split("-")[1].split(".bin")[0]) for file in myfiles])
    sortdexes=np.argsort(indexes)
    orderedindexes = indexes[sortdexes]
    myfiles = myfiles[sortdexes]
    return myfiles,orderedindexes


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


def assemble_StreamSnapShots(outfilename,t_sampling,attributes,inputdir):
    """
    This function assembles the stream snapshots into a single file
    """
    filenames,indexes=get_temp_snapshot_filenames(inputdir)
    with h5py.File(outfilename, 'w') as outfile:
        outfile.create_dataset("time_stamps",data=t_sampling)
        outfile.create_group("StreamSnapShots")
        for i in range(len(filenames)):
            phase_space = read_fortran_stream_binary_file(inputdir+filenames[i])
            outfile["StreamSnapShots"].create_dataset(str(indexes[i]),data=phase_space)
        for attr in attributes:
            outfile.attrs[attr] = attributes[attr]
            
    return None
    