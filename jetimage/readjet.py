import numpy as np
import ROOT

ROOT.gSystem.Load("~/MG5_aMC_v2_2_3/Delphes/libDelphes") #need to set location of libDelphes manually
import ROOT.Tower
import ROOT.Jet

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
    def __init__(self, jet_obj):
        self.pT = jet_obj.PT
        self.eta = jet_obj.Eta
        self.phi = jet_obj.Phi
        self.mass = jet_obj.Mass
        self.full_obj = jet_obj

        self.cells = self.read_constituents(self.full_obj)

    def read_constituents(self, obj):
        cells = []
        for k in range(obj.Constituents.GetEntriesFast()):
            const = obj.Constituents.At(k)
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
        chain = ROOT.TChain("Delphes")

        chain.Add(fname)
        self.tree_read = ROOT.ExRootTreeReader(chain)
        nev = self.tree_read.GetEntries()
        bjet = self.tree_read.UseBranch("Jet")
        btow = self.tree_read.UseBranch("Tower")

        self.events = []

        # print('hit')
        for entry in range(nev):
            self.tree_read.ReadEntry(entry)
            # check jet cutoff
            jets = []
            for j in range(bjet.GetEntriesFast()):
                jet = Jet(bjet.At(j))
                if len(jet)>0:
                    jets.append(jet)
            if len(jets)>0:
                self.events.append(jets)

    def __len__(self):
        return len(self.events)
    def __getitem__(self, item):
        return self.events[item]

# def root2jetary(fname):
#     chain = ROOT.TChain("Delphes")
#
#     chain.Add(fname)
#     treer = ROOT.ExRootTreeReader(chain)
#     nev = treer.GetEntries()
#     bjet = treer.UseBranch("Jet")
#     btow = treer.UseBranch("Tower")
#
#
#     events = []
#
#     for entry in range(nev):
#         treer.ReadEntry(entry)
#         #check jet cutoff
#         jets = []
#         for j in range(bjet.GetEntriesFast()):
#             jet = bjet.At(j)
#
#             const_lst = []
#             for k in range(jet.Constituents.GetEntriesFast()):
#                 const = jet.Constituents.At(k)
#                 if isinstance(const, Tower):
#                     const_lst.append([const.Eta, const.Phi, const.E])
#
#             jets.append(const_lst)
#
#         events.append(jets)
#
#     ev_ary = np.array(events)
#     return(ev_ary)

# root2jetary('ww2.root')

def dank():
    return RootEvents('ww2.root')