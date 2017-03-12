# import ROOT
#
# ROOT.gSystem.Load("~/MG5_aMC_v2_2_3/Delphes/libDelphes") #need to set location of libDelphes manually
# import ROOT.Tower, ROOT.Tower.Phi
# import ROOT.Jet
#
# chain = ROOT.TChain("Delphes")
# fname = 'ww2.root'
#
# chain.Add(fname)
# tree_read = ROOT.ExRootTreeReader(chain)
# nev = tree_read.GetEntries()
# # bjet = self.tree_read.UseBranch("Jet")
# btow = tree_read.UseBranch("Tower")
#
# print(btow.Phi)

# f = ROOT.Tfile('mytree.root', 'RECREATE')
# tree = ROOT.TTree('T', 'Just A Tree')
# tphi =  ROOT.Tower.Phi()
# tree.Branch('Tower.Phi', tphi)
#
# tphi = btow.Phi



# for entry in range(nev):
#     self.tree_read.ReadEntry(entry)
#     for j in range(btow.GetEntriesFast()):
#         jet = Jet(bjet.At(j))
#         if len(jet) > 0:
#             jets.append(jet)
#     if len(jets) > 0:
#         self.events.append(jets)


from readjet import RootEvents
import array
rev= RootEvents('ww11.root')
# rev[14][0].full_obj.TrimmedP4.Print()
# for ev in rev.events:
#     for j in ev:
#         j.full_obj.TrimmedP4.Print()