import numpy as np
import ROOT

ROOT.gSystem.Load("~/MG5_aMC_v2_2_3/Delphes/libDelphes") #need to set location of libDelphes manually
import ROOT.Tower
import ROOT.Jet
from ROOT import TLorentzVector
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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
            if (i>=0 and i<Np) and (j>=0 and j<Ne):
                im[i][j] += E/np.cosh(e)
        self.image_array = im
        # self.image_array = im/np.cosh(etas[None,:]) #replace energy with pT, invariant under translation
        self.dim = dim
        self.phiran = phiran
        self.etaran = etaran

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

    def __array__(self):
        return self.image_array

    def find_index(self, var, ran, dim):
        j = int(np.floor(dim[0] * (var[0] - ran[0][0]) / (ran[0][1] - ran[0][0])))
        i = dim[1] -1 - int(np.floor(dim[1] * (var[1] - ran[1][0]) / (ran[1][1] - ran[1][0])))
        return (j, i)
    def normalise(self):
        self.image_array /= np.linalg.norm(self.image_array)
        return self.image_array

def average_image(images):
    dim = images[0].dim
    im_ar = np.array([im.image_array for im in images])
    av_im = np.mean(im_ar, axis=0)
    #TODO average momenta
    final_image = images[0]
    final_image.image_array = av_im
    return final_image

def array2lorentz(arr):
    px, py, pz, E = arr
    return TLorentzVector(px, py, pz, E)

def theta2eta(theta):
    return -np.log(np.tan(theta/2))
def eta2theta(eta):
    return 2*np.arctan(np.exp(-eta))
def delta_phi(p1, p2):
    #calculate phi separation accounting for discontinuity
    #adapted from fastjet PseudoJet method
    dphi = p1-p2
    if dphi>np.pi:
        dphi -= 2*np.pi
    if dphi < -np.pi:
        dphi += 2*np.pi
    return dphi


