"""Usage:
    gen.py <in_file> [-w] [-o <out_file>] [--ptmin=<ptmin>] [--ptmax=<ptmax>]
    gen.py -h | --help

Arguments:
    <in_file> Root file to extract from

Options:
    -h --help        Show this screen
    -o <out_file>    file to save to
    -w               Write to new file rather than append
    --ptmin=<ptmin>  Minimum pT of jets in GeV [default: 250]
    --ptmax=<ptmax>  Maximum pT of jets in GeV [default: 300]

"""

from docopt import docopt
import numpy as np
import os
import h5py

from tqdm import tqdm
from jetimage.readjet import RootEvents, JetImage
from jetimage.analysis import plot_jet, average_image, image_set


def main(fname, output, ptmin, ptmax):
    print("Reading ROOT file.")
    jets = RootEvents(fname, ptmin=ptmin, ptmax=ptmax)
    images = []

    print("Pre-processing Jet images.")
    for jet in tqdm(jets):
        im0 = JetImage(jet)
        im0.sj2_rotate()
        im0.normalise()
        im0.flip()
        # plot_jet(im0)
        images.append(im0)

    # print("Images processed: {}".format(len(images)))

    if output is not None:
        print("Saving.")
        im_ar=image_set(images)
        # savefile=os.path.abspath(os.path.join('~','pt3proj','data', output))
        if len(im_ar) >0:
            savefile=os.path.abspath(output)
            with h5py.File(savefile, "a") as f:
                sname = set_name(f.keys())
                dset = f.create_dataset(sname, data=im_ar, chunks=True)
                                        #chunks=(100, im_ar.shape[1], im_ar.shape[2]))
        else:
            print("No images read.")
    else:
        #if no output is specified just plot average image
        av = average_image(images)
        plot_jet(av)

def set_name(keys):
    if len(keys) == 0:
        return "1"
    ints = [int(key) for key in keys]
    m = max(ints)
    return str(m+1)

if __name__ == '__main__':
    arguments = docopt(__doc__, help=True)
    main(arguments['<in_file>'], arguments['-o'],
         float(arguments['--ptmin']), float(arguments['--ptmax']))