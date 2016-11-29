print "Importing root ..."
import ROOT
print " ... done"

print "Loading libDelphes ..."

ROOT.gSystem.cd("~/MG5_aMC_v2_2_3/Delphes")
ROOT.gSystem.Load("libDelphes")
ROOT.gSystem.cd("~/pt3proj/ext/cl")
print " ... done"

print "Opening events.root ..."
f = ROOT.TFile("events.root")
print " ... done"

tree = f.Get("Delphes")

"""
Electron
Photon
GenJet
Jet
Muon
  fUniqueID
  fBits
  PT
  Eta
  Phi
  T
  Charge
  Particle
  IsolationVar
  IsolationVarRhoCorr
  SumPtCharged
  SumPtNeutral
  SumPtChargedPU
  SumPt
  @size
GenMissingET
ScalarHT
  HT
"""

cross_sec_fb = 24.3 # For the lepton(e/nu) n1 top process, according to MadGraph 
desired_int_lumi_eventPerfb = 20.3
generated_events = 500000

lumiSF = desired_int_lumi_eventPerfb*cross_sec_fb/generated_events

def scaleCut(scale,cut):
        return "("+str(scale)+")*("+cut+")"

print "Scanning tree ..."
#tree.Scan()

canvas = ROOT.TCanvas("x","Plots for "+str(desired_int_lumi_eventPerfb)+"/fb",1000,800)
canvas.Divide(3,3)

#browser = ROOT.TBrowser()

canvas.cd(1);
tree.Draw("Muon_size",scaleCut(lumiSF,"1"),"hist")

canvas.cd(2);
tree.Draw("Electron_size",scaleCut(lumiSF,"1"),"hist")

canvas.cd(3);
tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]<0"),"hist")
tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]>0"),"histsame")

canvas.cd(4);

canvas.cd(5);
tree.Draw("Muon.Charge[0]",scaleCut(lumiSF,"Muon_size==1"),"hist")

canvas.cd(6);
tree.Draw("Electron.Charge[0]",scaleCut(lumiSF,"Electron_size==1"),"hist")

canvas.cd(7);
tree.Draw("Electron.PT[0]+Muon.PT[0]",scaleCut(lumiSF,"Electron_size==1&&Muon_size==1&&Electron.Charge[0]>0&&Muon.Charge[0]<0"),"hist")
tree.Draw("Electron.PT[0]+Muon.PT[0]",scaleCut(lumiSF,"Electron_size==1&&Muon_size==1&&Electron.Charge[0]<0&&Muon.Charge[0]>0"),"histsame")

canvas.cd(8);
tree.Draw("Electron.Eta[0]-Muon.Eta[0]",scaleCut(lumiSF,"Electron_size==1&&Muon_size==1&&Electron.Charge[0]>0&&Muon.Charge[0]<0"),"hist")
tree.Draw("Electron.Eta[0]-Muon.Eta[0]",scaleCut(lumiSF,"Electron_size==1&&Muon_size==1&&Electron.Charge[0]<0&&Muon.Charge[0]>0"),"histsame")

canvas.cd(8);
thing="Electron.Phi[0]-Muon.Phi[0]"
tree.Draw(thing,scaleCut(lumiSF,"Electron_size==1&&Muon_size==1&&Electron.Charge[0]>0&&Muon.Charge[0]<0"),"hist")
tree.Draw(thing,scaleCut(lumiSF,"Electron_size==1&&Muon_size==1&&Electron.Charge[0]<0&&Muon.Charge[0]>0"),"histsame")

#browser = ROOT.TBrowser()

while True:
        pass
