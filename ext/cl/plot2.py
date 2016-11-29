randflat="((rand()+1.0-1.0)/2147483647)"

print "Importing root ..."
import ROOT
print " ... done"

print "Loading libDelphes ..."

ROOT.gROOT.ForceStyle()
ROOT.gSystem.cd("Delphes")
ROOT.gSystem.Load("libDelphes")
ROOT.gSystem.cd("..")
print " ... done"

print "Opening root files..."
class Dataset:
   def __init__(self, shortName, filename, crossSecPb, colour, legendText=None):
      if not legendText:
         legendText = shortName
      self.shortName=shortName
      self.legendText = legendText
      self.filename = filename
      self.colour = colour
      self.crossSecPb = crossSecPb
      self.crossSecFb = crossSecPb*1000.0
      self.rootfile = ROOT.TFile(filename)
      self.tree = self.rootfile.Get("Delphes")
      self.numberOfEvents = self.tree.GetEntries()
      self.hists = {}
      print "I believe that there are "+str(self.numberOfEvents)+" events in "+filename

   def project(self, shortHistName, texName, nbin, xmin, xmax, var, selection, desiredIntLumi_eventPerFb, legend):
      rootHistogramName = shortHistName+"_from_"+self.shortName
      self.hists[ shortHistName ] = ROOT.TH1F(rootHistogramName, texName, nbin, xmin, xmax)
      #hist.SetDirectory(0)  NEVER DO THIS
      #hist.SetStats(0)
      #self.hists[ shortHistName ].SetFillStyle(3001) # Light hatching
      self.hists[ shortHistName ].SetLineWidth(0)
      self.hists[ shortHistName ].SetLineColor(self.colour)
      self.hists[ shortHistName ].SetFillColor(self.colour)
      if legend:
          legend.AddEntry(self.hists[ shortHistName ], self.legendText)
      scaling = desiredIntLumi_eventPerFb*self.crossSecFb/self.numberOfEvents
      print "RESCALING for "+self.shortName+" is ",scaling," since desiredIntLumi_eventPerFb*self.crossSecFb/self.numberOfEvents = ",desiredIntLumi_eventPerFb, " * ",self.crossSecFb," / ",self.numberOfEvents
      if scaling<1:
         scaledSelection = "("+selection+")*(("+randflat+"<"+str(scaling)+")?1.0:0.0)"
      else: 
         scaledSelection = "("+selection+")*("+str(scaling)+")"
      self.tree.Project(rootHistogramName, var, scaledSelection, "hist")
      return self.hists[ shortHistName ]


def lesterProject(datasets, shortHistName, nbin, xmin, xmax, var, selection, desiredIntLumi_eventPerFb, texName=None):
    if not texName:
       texName = selection+' ; '+var
    leg = ROOT.TLegend(0.5,0.7,0.95,0.90)
    stack = ROOT.THStack(
                      shortHistName+"_stack",
                      texName,
                        )
    for dataset in datasets:
       hist = dataset.project(shortHistName, texName, nbin, xmin, xmax, var, selection, desiredIntLumi_eventPerFb,leg)
       #hist.Draw("hist")
       stack.Add(hist)
    #hist.Draw("hist")
    print "Printing the stack"
    stack.Print()
    print "Have printed the stack"
    print hist.GetFillColor()
    #stack.Draw("histf")
    stack.SetMinimum(0.5)
    return stack, leg

datasets = [
 Dataset("ttbar","ttbar_leptonic/Events/run_01/tag_1_delphes_events.root",crossSecPb=23.13*(1./3.), colour=ROOT.kRed),
 Dataset("ttbar","ttbar_leptonic/Events/run_02/tag_1_delphes_events.root",crossSecPb=23.12*(2./3.), colour=ROOT.kRed),
 #Dataset("tw","lester_test_tw_leptonicw/Events/run_01/tag_1_delphes_events.root",crossSecPb=12.08, colour=ROOT.kYellow), # Low stats run
 Dataset("tw","lester_test_tw_leptonicw/Events/run_02/tag_1_delphes_events.root",crossSecPb=12.08, colour=ROOT.kYellow), # High stats run 
 #Dataset("RPV","lester_test_leptonemu_n1_top/Events/run_01/tag_1_delphes_events.root", crossSecPb=0.02433*100, colour=ROOT.kGreen), # Note the *100 elevates the lambda'231 coupling to unity from 0.1.  This was a mixture of left and right smuons with masses around 320 GeV, and a neutralino around 90 GeV.  By the looks of later simulations (see below) this was 99.99% smuon right production.
 Dataset("RPV_1_150_250_X","lester_test_leptonemu_n1_top_1_150_250_X/Events/run_02/tag_3_delphes_events.root", crossSecPb=2.238e-41, colour=ROOT.kMagenta+2), # numbers mean, respectively: lambda'231, neutralino mass (GeV), smuon right mass (GeV), smuon left mass (GeV).  X indicates infinity
 Dataset("RPV_1_150_X_250","lester_test_leptonemu_n1_top_1_150_X_250/Events/run_02/tag_3_delphes_events.root", crossSecPb=1.674, colour=ROOT.kMagenta), 
 Dataset("RPV_1_90_300_X","lester_test_leptonemu_n1_top_1_90_300_X/Events/run_02/tag_3_delphes_events.root", crossSecPb=2.43e-41, colour=ROOT.kGreen+2), 
 Dataset("RPV_1_90_X_300","lester_test_leptonemu_n1_top_1_90_X_300/Events/run_02/tag_3_delphes_events.root", crossSecPb=2.843, colour=ROOT.kGreen), # Smuon left
]
print " ... done"




"""
Electron
Photon
GenJet
Jet
  BTag (0,1)
  Mass
  Flavor
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
MissingET
  MET
  Eta
  Phi
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

canvas.cd(2).SetLogy()
stack2, leg2 = lesterProject(datasets, "testHist002", 10, -2.5, 2.5, "Muon.Eta[0]","Muon_size>0&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
#tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]>0"),"histsame")
stack2.Draw("histf")
leg2.Draw()

canvas.cd(3).SetLogy()
stack3, leg3 = lesterProject(datasets, "testHist003", 10, -2.5, 2.5, "Muon.Eta[0]","Muon_size>0&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
#tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]>0"),"histsame")
stack3.Draw("histf")
leg3.Draw()


canvas.cd(1).SetLogy()
stack1, leg1 = lesterProject(datasets, "testHist001", 10, 0, 10, "Muon_size", "1==1", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack1.Draw("histf")
leg1.Draw()

osdflep="(Muon_size==1&&Electron_size==1&&Muon.Charge[0]*Electron.Charge[0]<0)"

canvas.cd(5).SetLogy()
stack5, leg5 = lesterProject(datasets, "testHist005", 10, -2.5, 2.5, "Muon.Eta[0]",osdflep+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
#tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]>0"),"histsame")
stack5.Draw("histf")
leg5.Draw()

canvas.cd(6).SetLogy()
stack6, leg6 = lesterProject(datasets, "testHist006", 10, -2.5, 2.5, "Muon.Eta[0]",osdflep+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
#tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]>0"),"histsame")
stack6.Draw("histf")
leg6.Draw()

canvas.cd(4).SetLogy()
stack4, leg4 = lesterProject(datasets, "testHist004", 10, 0, 10, "Jet_size", "1==1", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack4.Draw("histf")
leg4.Draw()

canvas.cd(7).SetLogy()
stack7, leg7 = lesterProject(datasets, "testHist007", 40, 0, 400, "MissingET.MET[0]", "MissingET_size==1", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack7.Draw("histf")
leg7.Draw()

osdflepandhighptmiss = "("+osdflep+"&&MissingET_size==1&&MissingET.MET[0]>150)"

canvas.cd(8).SetLogy()
stack8, leg8 = lesterProject(datasets, "testHist008", 10, -2.5, 2.5, "Muon.Eta[0]",osdflepandhighptmiss+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
#tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]>0"),"histsame")
stack8.Draw("histf")
leg8.Draw()

canvas.cd(9).SetLogy()
stack9, leg9 = lesterProject(datasets, "testHist009", 10, -2.5, 2.5, "Muon.Eta[0]",osdflepandhighptmiss+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
#tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]>0"),"histsame")
stack9.Draw("histf")
leg9.Draw()


canvas.Print("Moo1.pdf")

canvas.Clear()
canvas.Divide(3,3)


canvas.cd(1).SetLogy()
stack11, leg11 = lesterProject(datasets, "testHist011", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack11.Draw("histf")
leg11.Draw()

canvas.cd(2).SetLogy()
stack12, leg12 = lesterProject(datasets, "testHist012", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack12.Draw("histf")
leg12.Draw()

canvas.cd(3).SetLogy()
bit  = randflat+"<-1||"
stack13, leg13 = lesterProject(datasets, "testHist013", 24, -0.1, 1.1, randflat, bit+bit+"1==1", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack13.Draw("histf")
leg13.Draw()

canvas.cd(4) #.SetLogy()
stack14, leg14 = lesterProject(datasets, "testHist014", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack14.Draw("histf")
leg14.Draw()

canvas.cd(5) #.SetLogy()
stack15, leg15 = lesterProject(datasets, "testHist015", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack15.Draw("histf")
leg15.Draw()

canvas.cd(7).SetLogy()
stack17, leg17 = lesterProject(datasets, "testHist017", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack17.Draw("histf")
leg17.Draw()

canvas.cd(8).SetLogy()
stack18, leg18 = lesterProject(datasets, "testHist018", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb)
stack18.Draw("histf")
leg18.Draw()

canvas.Print("Moo2.pdf")

import sys
#while True:
#  import time
#  time.sleep(1)

sys.exit(1)

canvas.cd(1);
tree.Draw("Muon_size",scaleCut(lumiSF,"1"),"hist")

canvas.cd(2);
tree.Draw("Electron_size",scaleCut(lumiSF,"1"),"hist")

canvas.cd(3);
tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]<0"),"hist")
tree.Draw("Muon.Eta[0]",scaleCut(lumiSF,"Muon_size>0&&Muon.Charge[0]>0"),"histsame")


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
