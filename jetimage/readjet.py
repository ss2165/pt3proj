import ROOT
ROOT.gSystem.Load("~/MG5_aMC_v2_2_3/Delphes/libDelphes") #need to set location of libDelphes manually
import ROOT.Tower
import ROOT.Jet
from ROOT import TLorentzVector

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib import cm
import skimage.transform as sk
#TODO catch CompBase exception

class Cell:

    """Cell object. Captures energy, azimuthal angle and energy from tower cell constituent of a jet.
    Can be used as array of said variables"""
    def __init__(self, constituent_tower):
        self.E = constituent_tower.E
        self.eta = constituent_tower.Eta
        self.phi = constituent_tower.Phi
        self.pT = self.E/np.cosh(self.eta)

    # def __float__(self):
    #     return self.__array__()

    def __array__(self):
        return np.array([self.eta, self.phi, self.E], dtype='float_')

    def __len__(self):
        return len(self.__array__())

    def notempty(self):
        return len(self)>1

class Jet:
    def __init__(self, jet_obj, subjet_array):
        self.pT = jet_obj.PT
        self.eta = jet_obj.Eta
        self.phi = jet_obj.Phi
        self.mass = jet_obj.Mass
        self.momentum = TLorentzVector()
        self.momentum.SetPtEtaPhiM(self.pT, self.eta, self.phi, self.mass)
        self.full_obj = jet_obj
        momenta = [array2lorentz(a) for a in subjet_array]
        self.trimmed_momentum = momenta[0]
        self.sj_momentum = momenta[1:]
        self.cells = self.read_constituents(self.full_obj.Constituents)

    def read_constituents(self, constituents):
        cells = []
        for const in constituents:
            if isinstance(const, ROOT.Tower):
                c = Cell(const)
                if c.notempty():
                    cells.append(c)
        return cells



    def __array__(self):
        ar = np.zeros((len(self.cells), 3))
        for i, cell in enumerate(self.cells):
            ar[i,:] = np.array(cell)
        return ar
    def __len__(self):
        return len(self.cells)
    def __getitem__(self, item):
        return self.cells[item]



class RootEvents:
    def __init__(self, fname, ptmin, ptmax):
        trimp4 = self.extract_subjets(fname)
        chain = ROOT.TChain("Delphes")

        chain.Add(fname)
        self.tree_read = ROOT.ExRootTreeReader(chain)
        nev = self.tree_read.GetEntries()
        bjet = self.tree_read.UseBranch("Jet")

        btow = self.tree_read.UseBranch("Tower")

        self.events = []

        # print('hit')
        for entry, subjets in trimp4:
            self.tree_read.ReadEntry(entry)
            # check jet cutoff
            # jets = []
            # for j in range(bjet.GetEntriesFast()):
            #     jet = Jet(bjet.At(j))
            #     if len(jet)>0:
            #         jets.append(jet)
            # if len(jets)>0:
            #     self.events.append(jets)

            #take leading jet and its subjets
            j = bjet.At(0)
            if ptmin <= j.PT <= ptmax:
                self.events.append(Jet(j, subjets))

    def extract_subjets(self, fname):
        """Run root macro to extract 4 momenta of subjets"""
        script_name = 'trim_macro.cxx'
        root_command = ['root', '-q', '-b', "{}(\"{}\")".format(script_name, fname)]
        proc = subprocess.Popen(root_command, stdout=subprocess.PIPE)
        tmp = proc.stdout.read()
        seek_text = 'Processing {}("{}")...\n'.format(script_name, fname)
        idx = tmp.find(seek_text)
        ar = np.fromstring(tmp[idx + len(seek_text):-1], dtype=float, sep=",")
        trimp4 = [(int(ar[x]), ar[x + 1: x + 21].reshape((5, 4))) for x in range(0,ar.size,21)]
        return trimp4

    def __len__(self):
        return len(self.events)
    def __getitem__(self, item):
        return self.events[item]

class JetImage:
    def __init__(self, jet, dim=(25,25), etaran=(-1.25,1.25), phiran=(-1.25,1.25)):
        self.momentum = jet.momentum
        self.subjets = jet.sj_momentum
        jet = np.array(jet)  # convert to array of cell values
        # eta, phi, E = jet.T

        # number of pixels per axis
        Ne, Np = dim

        im = np.zeros((Np, Ne))  # array containing image
        # phi_y = np.arange(-phiran, phiran, step)
        # eta_x = np.arange(-etaran, etaran, step)

        etas = np.linspace(etaran[0], etaran[1], Ne)

        self.ec, self.pc = (self.subjets[1].Eta(), self.subjets[1].Phi())
        for e, p, E in jet:
            # find pixel closest to coordinate values
            # j = self.find_index(e, etaran, Ne)
            # i = Np - 1 - self.find_index(p, phiran, Np) #flip so that eta axis goes vertically
            translated = self.centre(e, p)
            j, i = self.find_index(translated, (etaran, phiran), dim)
            if (0 <= i < Np) and (0 <= j < Ne):
                im[i][j] += E/np.cosh(e) #pT, invariant under translations
        self.image_array = im

        self.dim = dim
        self.phiran = phiran
        self.etaran = etaran
        #store discrete sets of coordinates
        self.etavals = self.generate_coords(etaran, Ne, Np)
        self.phivals = np.flipud(self.generate_coords(phiran, Np, Ne).T)

    @staticmethod
    def generate_coords(rang, N, N_other, flip=False):
        #generate array of coordinate values in shape of image, over range=rang,
        #N steps, N_other is size of other axis
        step = (rang[1]-rang[0])/N
        vals= np.linspace(rang[0]+step/2, rang[1]-step/2, N)

        return np.tile(vals, (N_other,1))

    def centre(self, e, p):
        #translate coordinates to central value
        return e-self.ec, delta_phi(p,self.pc)

    def sj2_rotate(self):
        #rotate second subjet to -pi.2
        e, p = (self.subjets[2].Eta(), self.subjets[2].Phi())
        if e < -10 or p < -10:
            e, p = self.pca_dir()
        else:
            e, p = self.centre(e, p)
        angle = np.arctan(p / e) + np.pi/2

        if (-np.sin(angle) * e + np.cos(angle) * p) > 0:
            angle += -np.pi
        self.image_array = sk.rotate(self.image_array, np.rad2deg(-angle), order = 3)
        return self.image_array

    def pca_dir(self):
        print("pca used")
        I = self.image_array
        norm=np.sum(I)
        mux = np.sum(self.etavals*I)/norm
        muy = np.sum(self.phivals*I)/norm

        x = self.etavals-mux
        y = self.phivals-muy
        xbar = np.sum(x*I)
        ybar = np.sum(y*I)
        x2bar = np.sum(x*x*I)
        y2bar = np.sum(y*y*I)
        xybar = np.sum(x*y*I)

        sigmax2 = x2bar/norm - mux**2
        sigmay2 = y2bar/norm - muy**2
        sigmaxy = xybar/norm - mux*muy

        lamb_min = 0.5*( sigmax2 + sigmay2 - np.sqrt( (sigmax2-sigmay2)*(sigmax2-sigmay2) + 4*sigmaxy*sigmaxy) )
        lamb_max = 0.5 * (sigmax2 + sigmay2 + np.sqrt((sigmax2 - sigmay2) * (sigmax2 - sigmay2) + 4 * sigmaxy * sigmaxy))

        dir_x = sigmax2 + sigmaxy - lamb_min
        dir_y = sigmay2 + sigmaxy - lamb_min

        #The first PC is only defined up to a sign.  Let's have it point toward the side of the jet with the most energy.
        dotprod = dir_x*x + dir_y*y
        if np.sum(I[dotprod>0]) > np.sum(I[dotprod<0]):
            dir_x *= -1
            dir_y *= -1

        return dir_x, dir_y

    def flip(self, side='r'):
        """
            Flips image so that most energetic half is on 'side'
            """
        jet = self.image_array
        weight = jet.sum(axis=0)

        halfway = jet.shape[0] / 2.
        l, r = np.int(np.floor(halfway)), np.int(np.ceil(halfway))
        l_weight, r_weight = np.sum(weight[:l]), np.sum(weight[r:])

        if 'r' in side:
            if r_weight > l_weight:
                return jet
            return np.fliplr(jet)
        else:
            if l_weight > r_weight:
                return jet
            return np.fliplr(jet)

    def draw(self):
        f = plt.imshow(self.image_array, norm=LogNorm(vmin=0.00001, vmax=1), interpolation='nearest')
        f.set_cmap('nipy_spectral')
        plt.colorbar()
        plt.show()

    def __array__(self):
        return self.image_array

    @staticmethod
    def find_index(var, ran, dim):
        j = int(np.floor(dim[0] * (var[0] - ran[0][0]) / (ran[0][1] - ran[0][0])))
        i = dim[1] -1 - int(np.floor(dim[1] * (var[1] - ran[1][0]) / (ran[1][1] - ran[1][0])))
        return j, i

    def normalise(self):
        self.image_array /= np.linalg.norm(self.image_array)
        return self.image_array


def plot_diff(phiran, etaran, image_array, output_name=None, extr=None, title='', cmap=cm.PRGn_r):
    """
    Adapted from Luke de Oliveira, & Michela Paganini. (2017). lukedeo/adversarial-jets: Initial Release [Data set]. Zenodo. http://doi.org/10.5281/zenodo.400708

        Function to help you visualize the difference between two sets of jet images on a linear scale
        Args:
        -----
           content : numpy array of dimensions 25x25, first arg to imshow, content of the image
                     e.g.: generated_images.mean(axis=0) - real_images.mean(axis=0) --> difference between avg generated and avg Pythia image
                           etc...
           output_name : string, name of the output file where the plot will be saved. Note: it will be located in ../plots/
           extr : (default = None) float, magnitude of the upper and lower bounds of the pixel intensity scale before saturation (symmetric around 0)
           title : (default = '') string, title of the plot, to be displayed on top of the image
           cmap : (default = matplotlib.cm.PRGn_r) matplotlib colormap, ideally white in the middle
        Outputs:
        --------
           no function returns
           saves file in ../plots/output_name
        """
    fig, ax = plt.subplots(figsize=(6, 6))
    extent = phiran + etaran
    if extr is None:
        extr = max(abs(image_array.min()), abs(image_array.max()))
    im = ax.imshow(image_array,
                   interpolation='nearest', norm=Normalize(vmin=-extr, vmax=+extr), extent=extent,
                   cmap=cmap)
    plt.colorbar(im, fraction=0.05, pad=0.05)
    plt.xlabel(r'[Transformed] Pseudorapidity $(\eta)$')
    plt.ylabel(r'[Transformed] Azimuthal Angle $(\phi)$')
    plt.title(title)
    if output_name is None:
        plt.show()
    else:
        plt.savefig(output_name)


def discrete_pT(image_array, phivals):
    cos_term = image_array * np.cos(phivals)
    sin_term = image_array * np.sin(phivals)
    return np.sqrt(np.sum(cos_term)**2 + np.sum(sin_term)**2)

def discrete_m(image_array, etavals, phivals):
    sinh_term = image_array * np.sinh(etavals)
    return np.sqrt(np.sum(image_array) ** 2 - discrete_pT(image_array, phivals) - np.sum(sinh_term) ** 2)


def tau1(image_array, etavals, phivals, dim):
    """"
    Calculates the normalized tau1 from a pixelated jet image
    Args:
    -----
        jet_image: numpy ndarray of dim (1, 25, 25)
    Returns:
    --------
        float, normalized jet tau1
    """
    # find coordinate of most energetic pixel, then use formula to compute tau1
    eta = etavals
    phi = phivals
    tau1_axis_eta = eta.ravel()[np.argmax(image_array)]
    tau1_axis_phi = phi.ravel()[np.argmax(image_array)]
    t1 = np.sum(image_array*
                  np.sqrt(np.square(tau1_axis_eta - eta) + np.square(
                      [delta_phi(tau1_axis_phi, p) for p in phi.ravel()]).reshape(dim))
                  )
    return t1 / np.sum(image_array)  # normalize by the total intensity

def tau2(image_array, etavals, phivals):
    """
    Calculates the normalized tau2 from a pixelated jet image
    Args:
    -----
        jet_image: numpy ndarray of dim (1, 25, 25)
    Returns:
    --------
        float, normalized jet tau2
    Notes:
    ------
        slow implementation
    """
    jet_image = image_array
    eta= etavals
    phi= phivals
    proto = np.array(zip(jet_image[jet_image != 0],
                         eta[jet_image != 0],
                         phi[jet_image != 0]))

    while len(proto) > 2:
        candidates = [
            (
                (i, j),
                (min(pt1, pt2) ** 2) * ((eta1 - eta2) ** 2 + (phi1 - phi2) ** 2)
            )
            for i, (pt1, eta1, phi1) in enumerate(proto)
            for j, (pt2, eta2, phi2) in enumerate(proto)
            if j > i
            ]

        index, value = zip(*candidates)
        pix1, pix2 = index[np.argmin(value)]
        if pix1 > pix2:
            # swap
            pix1, pix2 = pix2, pix1

        (pt1, eta1, phi1) = proto[pix1]
        (pt2, eta2, phi2) = proto[pix2]

        e1 = pt1 / np.cosh(eta1)
        e2 = pt2 / np.cosh(eta2)
        choice = e1 > e2

        eta_add = (eta1 if choice else eta2)
        phi_add = (phi1 if choice else phi2)
        pt_add = (e1 + e2) * np.cosh(eta_add)

        proto[pix1] = (pt_add, eta_add, phi_add)

        proto = np.delete(proto, pix2, axis=0).tolist()

    (_, eta1, phi1), (_, eta2, phi2) = proto
    np.sqrt(np.square(eta - eta1) + np.square(phi - phi1))

    grid = np.array([
        np.sqrt(np.square(eta - eta1) + np.square(phi - phi1)),
        np.sqrt(np.square(eta - eta2) + np.square(phi - phi2))
    ]).min(axis=0)

    return np.sum(jet_image * grid) / np.sum(jet_image)  # normalize by the total intensity


def tau21(image_array, etavals, phivals, dim):
    """
    Calculates the tau21 from a pixelated jet image using the functions above
    Args:
    -----
        jet_image: numpy ndarray of dim (1, 25, 25)
    Returns:
    --------
        float, jet tau21
    Notes:
    ------
        slow implementation
    """
    t1 = tau1(image_array, etavals, phivals, dim)
    if tau1 <= 0:
        return 0
    else:
        t2 = tau2(image_array, etavals, phivals)
        return t2 / t1



def plot_jet(image_array, etaran=(-1.25,1.25), phiran=(-1.25,1.25), output_name=None, vmin=1e-6, vmax=300, title=''):
    """
        Adapted from Luke de Oliveira, & Michela Paganini. (2017). lukedeo/adversarial-jets: Initial Release [Data set]. Zenodo. http://doi.org/10.5281/zenodo.400708
        Function to help you visualize a jet image on a log scale
        Args:
        -----
           content : numpy array of dimensions 25x25, first arg to imshow, content of the image
                     e.g.: generated_images.mean(axis=0) --> the average generated image
                           real_images.mean(axis=0) --> the average Pythia image
                           generated_images[aux_out == 1].mean(axis=0) --> the average generated image labeled as real by the discriminator
                           etc...
           output_name : string, name of the output file where the plot will be saved. Note: it will be located in ../plots/
           vmin : (default = 1e-6) float, lower bound of the pixel intensity scale before saturation
           vmax : (default = 300) float, upper bound of the pixel intensity scale before saturation
           title : (default = '') string, title of the plot, to be displayed on top of the image
        Outputs:
        --------
           no function returns
           saves file in ../plots/output_name
        """
    fig, ax = plt.subplots(figsize=(7, 6))
    extent = phiran + etaran
    im = ax.imshow(image_array, interpolation='nearest', norm=LogNorm(vmin=vmin, vmax=vmax), extent=extent)
    cbar = plt.colorbar(im, fraction=0.05, pad=0.05)
    cbar.set_label(r'Pixel $p_T$ (GeV)', y=0.85)
    plt.xlabel(r'[Transformed] Pseudorapidity $(\eta)$')
    plt.ylabel(r'[Transformed] Azimuthal Angle $(\phi)$')
    plt.title(title)
    if output_name is None:
        plt.show()
    else:
        plt.savefig(output_name)


def average_image(images):
    #avereage image from list of JetImage
    im_ar = image_set(images)
    av_im = np.mean(im_ar, axis=0)
    #TODO average momenta
    final_image = images[0]
    final_image.image_array = av_im
    return final_image

def maxim(images):
    #return the JetImage object with highest discrete pT in list
    return max(images, key=lambda item: discrete_pT(item.image_array, item.phivals))

def image_set(images):
    #convert list of JetImages to array of image arrays shape (N_images, Nphi, Neta)
    return np.array([image.image_array for image in images])

def array2lorentz(arr):
    #convert (1,4) array lorentz vector to TLorentzVector object
    px, py, pz, E = arr
    return TLorentzVector(px, py, pz, E)

def theta2eta(theta):
    return -np.log(np.tan(theta/2))
def eta2theta(eta):
    return 2*np.arctan(np.exp(-eta))
def delta_phi(p1, p2):
    #calculate phi separation accounting for discontinuity
    return np.arccos(np.cos(abs(p1 - p2)))


