"""Usage:
    datacheck.py <data_file>  <dataset> [-p]
    datacheck.py -h | --help

Arguments:
    <data_file>  HDF5 file to extract from
    <dataset>  set to use

Options:
    -h --help        Show this screen
    -p               Plot average

"""
from docopt import docopt
import h5py
import numpy as np
from jetimage.analysis import  plot_jet, average_image

def main(fname, plot, dset):
    arrays = []
    with h5py.File(fname, 'r') as f:
        print("(Key, shape)")
        for key in f.keys():
            # print(key, f[key][:].shape)
            arrays.append(f[key][:])

    images = np.vstack(arrays)
    with h5py.File(dset, 'w') as f:
        dset = f.create_dataset('images', data=images, chunks=True)
    if plot:
        av = average_image(ar)
        plot_jet(av)



if __name__ == '__main__':
    arguments = docopt(__doc__, help=True)
    print(arguments)

    main(arguments['<data_file>'], arguments['-p'], arguments['<dataset>'])