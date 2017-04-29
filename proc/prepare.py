"""
combine wprime and qcd datasets in to friendly form for training
"""
import h5py
import numpy as np

def main():
    with h5py.File('../data/comb_qcd.hdf', 'r') as f:
        qims = f['images'][:]
    l1 = qims.shape[0]
    images = np.zeros((l1*2,25, 25))

    with h5py.File('../data/comb_wprime.hdf', 'r') as f:
        wims = f['images'][:l1]

    images[:l1, :, :] = qims
    images[l1:, :, :] = wims
    signal = np.zeros(l1*2)
    signal[l1:] += 1

    with h5py.File('../data/prepared.hdf', 'w') as f:
        dset = f.create_dataset('image', data=images)
        sigs = f.create_dataset('signal', data=signal)
main()