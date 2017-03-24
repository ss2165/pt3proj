"""Usage: trials.py <file> [--out=<out_file>]

Arguments:
    <file> Root file to extract from

Options:
    --out=<out_file> file to save to
"""

from docopt import docopt
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns

from readjet import RootEvents, JetImage, average_image


def draw_jet(jet, etaran=6.0, phiran=4.0, step=0.1):
    """
    jet --- Jet object to be drawn
    etaran, phiran -- range of plot in eta, phi
    step -- pizel step size in coordinates
    """
    jet = np.array(jet) #convert to array of cell values
    eta, phi, E = jet.T

    # number of pixels per axis
    Ne = int(2 * etaran / step)
    Np = int(2 * phiran / step)

    im = np.zeros((Ne, Np))  # array containing image
    phi_y = np.arange(-phiran, phiran, step)
    eta_x = np.arange(-etaran, etaran, step)

    for p, e, E in zip(phi, eta, E):
        # find pixel closest to coordinate values
        i = int(np.floor((e + etaran) / step))
        j = int(np.floor((p + phiran) / step))

        im[i][j] += E

    # im = im/np.cosh(eta_x[:, None]) #convert from E to pT using pt = E/cosh(eta)

    plt.imshow(im, norm=LogNorm(vmin=0.00001, vmax=1), interpolation='nearest') #from jettools.py
    # plt.plot(eta, phi, 'x')
    # plt.plot(J[0][:2], J[1][:2], 'xr')
    plt.show()




def main(fname, output):

    ents = RootEvents(fname, 100, 300)
    images = []
    for ev in ents:
        im0 = JetImage(ev, dim=(25,25), phiran=(-1.25, 1.25), etaran=(-1.25, 1.25))
        im0.sj2_rotate()
        im0.normalise()
        images.append(im0)
    print(len(images))
    av = average_image(images)
    av.draw()



def plot_jet(rec, title='Jet Image', log=True):
    fig = plt.figure(figsize=(8, 8), dpi=100)
    ax = fig.add_subplot(111)
    if log:
        im = ax.imshow(rec, norm=LogNorm(vmin=0.00001, vmax=1), interpolation='nearest')
    else:
        im = ax.imshow(rec, interpolation='nearest')
    plt.title(r'' + title)
    return fig


if __name__ == '__main__':
    arguments = docopt(__doc__)
    main(arguments['<file>'], arguments['--out'])