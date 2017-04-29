"""Usage:
    combine.py <data_file>  <dataset>
    combine.py -h | --help

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

def main(fname, output):
    arrays = []
    with h5py.File(fname, 'r') as f:
        print("(Key, shape)")
        for key in f.keys():
            # print(key, f[key][:].shape)
            arrays.append(f[key][:])

    images = np.vstack(arrays)
    with h5py.File(output, 'w') as f:
        d_set = f.create_dataset('images', data=images, chunks=True)



if __name__ == '__main__':
    arguments = docopt(__doc__, help=True)
    print(arguments)

    main(arguments['<data_file>'], arguments['<dataset>'])