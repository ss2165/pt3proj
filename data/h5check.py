import h5py
import numpy as np

with h5py.File('delphes_images.hdf5', 'r') as f:
    for key in f.keys():
        print(key, f[key][:].shape)