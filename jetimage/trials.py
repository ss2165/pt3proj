import numpy as np

from rootpy.io import root_open
from rootpy.tree import Tree, TreeModel, FloatCol, FloatArrayCol

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
fname = 'ww2.root'

from readjet import RootEvents


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




def main():
    # with root_open(fname) as f:
    #     d = f.Delphes.to_array(branches = ['Tower.Eta', 'Tower.Phi', 'Tower.E'])
    #     J = f.Delphes.to_array(branches=['Jet.Eta', 'Jet.Phi', 'Jet.PT'])
    #     nevents = len(d)
    #     d = d[2] #first event only, for now
    #     print(d)
    #     J = J[0]
    #     eta, phi, E = d
    #     #range of coordinates
    #     etran = 6.0
    #     phran = 4.0
    #     step = 0.1
    #     #number of pixels per axis
    #     Ne = int(2*etran/step)
    #     Np = int(2*phran/step)
    #
    #     im = np.zeros((Ne, Np)) #array containing image
    #     phi_y = np.arange(-phran, phran, step)
    #     eta_x = np.arange(-etran, etran, step)
    #
    #     for p, e, E in zip(phi, eta, E):
    #         #find pixel closest to coordinate values
    #         i = int(np.floor((e + etran)/step))
    #         j = int(np.floor((p + phran)/step))
    #
    #         im[i][j] += E
    #
    #     # im = im/np.cosh(eta_x[:, None]) #convert from E to pT using pt = E/cosh(eta)
    #
    #     plt.imshow(im, interpolation='none')
    #     # plt.plot(eta, phi, 'x')
    #     # plt.plot(J[0][:2], J[1][:2], 'xr')
    #     plt.show()

    ents = RootEvents('ww2.root')
    ev2 = ents[2]
    j1_ev2 = ev2[1]

    draw_jet(j1_ev2)

if __name__ == '__main__':
    main()