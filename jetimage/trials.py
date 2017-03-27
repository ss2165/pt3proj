"""Usage:
    trials.py <in_file> [-o=<out_file>] [--ptmin=<ptmin>] [--ptmax=<ptmax>]
    trials.py -h | --help

Arguments:
    <in_file> Root file to extract from

Options:
    -h --help        Show this screen
    -o=<out_file>    file to save to
    --ptmin=<ptmin>  Minimum pT of jets in GeV [default: 250]
    --ptmax=<ptmax>  Maximum pT of jets in GeV [default: 300]

"""

from docopt import docopt
import numpy as np
import os

from readjet import RootEvents, JetImage, average_image, plot_jet, image_set

def main(fname, output, ptmin, ptmax):
    jets = RootEvents(fname, ptmin=ptmin, ptmax=ptmax)
    images = []

    for jet in jets:
        im0 = JetImage(jet)
        im0.sj2_rotate()
        im0.normalise()
        im0.flip()
        # plot_jet(im0)
        images.append(im0)

    print("Images processed: {}".format(len(images)))
    av = average_image(images)
    plot_jet(av)
    if output is not None:
        im_ar=image_set(images)
        savefile=os.path.abspath(os.path.join('..','data', output))
        np.save(savefile, im_ar)



if __name__ == '__main__':
    arguments = docopt(__doc__, help=True)
    main(arguments['<in_file>'], arguments['-o'],
         float(arguments['--ptmin']), float(arguments['--ptmax']))