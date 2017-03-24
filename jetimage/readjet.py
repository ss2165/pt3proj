import numpy as np
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
    def __init__(self, jet, dim=(25,25), phiran=(-1.25,1.25), etaran=(-1.25,1.25)):
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
        self.etavals = self.generate_coords(etaran, Ne, Np).T
        self.phivals = self.generate_coords(phiran, Np, Ne)


    @staticmethod
    def generate_coords(rang, N, N_other):
        step = (rang[1]-rang[0])/N
        vals= np.linspace(rang[0]+step/2, rang[1]-step/2, N)
        return np.tile(vals, (N_other,1))

    def centre(self, e, p):
        return e-self.ec, delta_phi(p,self.pc)

    def sj2_rotate(self):
        e, p = (self.subjets[2].Eta(), self.subjets[2].Phi())
        e, p = self.centre(e, p)
        angle = np.arctan(p / e) + np.pi/2

        if (-np.sin(angle) * e + np.cos(angle) * p) > 0:
            angle += -np.pi
        self.image_array = sk.rotate(self.image_array, np.rad2deg(-angle), order = 3)
        return self.image_array

    def draw(self):
        f = plt.imshow(self.image_array, norm=LogNorm(vmin=0.00001, vmax=1), interpolation='nearest')
        f.set_cmap('nipy_spectral')
        plt.colorbar()
        plt.show()

    def plot(self, output_name=None, vmin=1e-6, vmax=300, title=''):
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
        extent = self.phiran+self.etaran
        im = ax.imshow(self.image_array, interpolation='nearest', norm=LogNorm(vmin=vmin, vmax=vmax), extent=extent)
        cbar = plt.colorbar(im, fraction=0.05, pad=0.05)
        cbar.set_label(r'Pixel $p_T$ (GeV)', y=0.85)
        plt.xlabel(r'[Transformed] Pseudorapidity $(\eta)$')
        plt.ylabel(r'[Transformed] Azimuthal Angle $(\phi)$')
        plt.title(title)
        if output_name is None:
            plt.show()
        else:
            plt.savefig(output_name)

    def plot_diff(self, output_name=None, extr=None, title='', cmap=cm.PRGn_r):
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
        extent = self.phiran+self.etaran
        if extr is None:
            extr = max(abs(self.image_array.min()), abs(self.image_array.max()))
        im = ax.imshow(self.image_array,
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

    def discrete_pT(self):
        cos_term = self.image_array*np.cos(self.phivals)
        sin_term = self.image_array*np.sin(self.phivals)
        return np.sqrt(np.sum(cos_term)**2 + np.sum(sin_term)**2)

    def discrete_m(self):
        sinh_term = self.image_array*np.sinh(self.etaval)
        return np.sqrt(np.sum(self.image_array)**2 - self.pT2() - np.sum(sinh_term)**2)

    def tau1(self):
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
        jet_image = self.image_array
        eta = self.etavals
        phi = self.phivals
        tau1_axis_eta = eta.ravel()[np.argmax(jet_image)]
        tau1_axis_phi = phi.ravel()[np.argmax(jet_image)]
        tau1 = np.sum(jet_image*
                      np.sqrt(np.square(tau1_axis_eta - eta) + np.square(
                          [delta_phi(tau1_axis_phi, p) for p in phi.ravel()]).reshape(self.dim))
                      )
        return tau1 / np.sum(jet_image)  # normalize by the total intensity

    def tau2(self):
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
        jet_image = self.image_array
        eta=self.etavals
        phi=self.phivals
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

    def tau21(self):
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
        tau1 = self.tau1()
        if tau1 <= 0:
            return 0
        else:
            tau2 = self.tau2()
            return tau2 / tau1

def average_image(images):
    im_ar = np.array([im.image_array for im in images])
    av_im = np.mean(im_ar, axis=0)
    #TODO average momenta
    final_image = images[0]
    final_image.image_array = av_im
    return final_image

def image_set(images):
    return np.array([image.image_array for image in images])

def array2lorentz(arr):
    px, py, pz, E = arr
    return TLorentzVector(px, py, pz, E)

def theta2eta(theta):
    return -np.log(np.tan(theta/2))
def eta2theta(eta):
    return 2*np.arctan(np.exp(-eta))
def delta_phi(p1, p2):
    #calculate phi separation accounting for discontinuity
    return np.arccos(np.cos(abs(p1 - p2)))


