import numpy as np
import ROOT

ROOT.gSystem.Load("~/MG5_aMC_v2_2_3/Delphes/libDelphes") #need to set location of libDelphes manually
import ROOT.Tower
import ROOT.Jet
from ROOT import TLorentzVector
import subprocess
import numpy as np
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
        self.momentum.SetPtEtaPhiE(self.pT, self.eta, self.phi, jet_obj.E)
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
    def __init__(self, fname):
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
            self.events.append(Jet(bjet.At(0), subjets))

    def extract_subjets(self, fname):
        """Run root macro to extract 4 momenta of subjets"""
        script_name = 'trim_macro.cxx'
        root_command = ['root', '-q', '-b', "{}(\"{}\")".format(script_name, fname)]
        proc = subprocess.Popen(root_command, stdout=subprocess.PIPE)
        tmp = proc.stdout.read()
        seek_text = 'Processing {}("{}")...\n'.format(script_name, fname)
        idx = tmp.find(seek_text)
        ar = np.fromstring(tmp[idx + len(seek_text):-1], dtype=float, sep=",")
        trimp4 = [(int(ar[21 * x]), ar[21 * x + 1: 21 * (x + 1)].reshape((5, 4))) for x in range(ar.size / 21)]
        return trimp4

    def __len__(self):
        return len(self.events)
    def __getitem__(self, item):
        return self.events[item]

def array2lorentz(arr):
    px, py, pz, E = arr
    return TLorentzVector(px, py, pz, E)


