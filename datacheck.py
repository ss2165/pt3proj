"""Usage:
    datacheck.py <data_file>  <dataset>
    datacheck.py -h | --help

Arguments:
    <data_file>  HDF5 file to extract from
    <dataset>  set to use

Options:
    -h --help        Show this screen

"""
from docopt import docopt
import h5py
import numpy as np
from jetimage.analysis import  plot_jet, average_image

def main(fname, keys, dset):
    with h5py.File(fname, 'r') as f:
        print("(Key, shape)")
        for key in f.keys():
            print(key, f[key][:].shape)
        if dset is not None:
            ar = f[dset][:]
    if not keys:
        av = average_image(ar)
        plot_jet(av)



if __name__ == '__main__':
    arguments = docopt(__doc__, help=True)
    print(arguments)

    main(arguments['<data_file>'], False, arguments['<dataset>'])