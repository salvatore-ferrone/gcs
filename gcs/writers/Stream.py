import h5py 
import datetime


author = "Salvatore Ferrone"
author_affiliation = "Sapienza University of Rome"
author_email = "salvatore.ferrone@uniroma1.it"
description = "Integrating star-particles with in a globular cluster"




def stream(filename,stream,tesc,attributes):
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
        dset = f.create_dataset("phase_space", data=stream)
        dset.attrs['phase_space_columns'] = phase_space_columns
        f.create_dataset("tesc", data=tesc)
        for attr in attributes:
            f.attrs[attr] = attributes[attr]
    return None

