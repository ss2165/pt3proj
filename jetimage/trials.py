import numpy as np

from rootpy.io import root_open
from rootpy.tree import Tree, TreeModel, FloatCol, FloatArrayCol

import matplotlib.pyplot as plt
fname = 'ww.root'


def main():
    with root_open(fname) as f:
        d = f.Delphes.to_array(branches = ['Tower.Eta', 'Tower.Phi', 'Tower.E'])[0]
        J = f.Delphes.to_array(branches=['Jet.Eta', 'Jet.Phi', 'Jet.PT'])[0]
        eta, phi, E = d
        #range of coordinates
        etran = 6.0
        phran = 4.0
        step = 0.1
        #number of pixels per axis
        Ne = int(2*etran/step)
        Np = int(2*phran/step)

        im = np.zeros((Ne, Np)) #array containing image
        phi_y = np.arange(-phran, phran, step)
        eta_x = np.arange(-etran, etran, step)

        for p, e, E in zip(phi, eta, E):
            #find pixel closest to coordinate values
            i = int(np.floor((e + etran)/step))
            j = int(np.floor((p + phran)/step))

            im[i][j] += E

        # im = im/np.cosh(eta_x[:, None]) #convert from E to pT using pt = E/cosh(eta)

        plt.imshow(im, interpolation='none')
        # plt.plot(eta, phi, 'x')
        # plt.plot(J[0][:2], J[1][:2], 'xr')
        plt.show()

main()