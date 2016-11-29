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

def subtractHistograms(one, two):
   ans = one.Clone()
   ans.SetName(one.GetName()+"_minus_"+two.GetName())
   ans.Add(two,-1.0)
   return ans

def addHistograms(one, two):
   ans = one.Clone()
   ans.SetName(one.GetName()+"_plus_"+two.GetName())
   ans.Add(two)
   return ans

def processHistograms(one, two, func, err):
   ans = one.Clone()
   ans.SetName(one.GetName()+"_combinedwith_"+two.GetName())
   for bin in range(ans.GetNbinsX()+2): #+2 accounts for over- and under-flow bins
      val1 = one.GetBinContent(bin)
      val2 = two.GetBinContent(bin)
      ans.SetBinContent(bin,func(val1,val2))
      ans.SetBinError(bin,err(val1,val2))
   return ans

class Bundle:
   def draw(self):
      self.stack.Draw("histf")
      self.stacktot.Draw("hist,p,same")
      for i in self.nonStackable:
         i.Draw("hist,same")
      self.leg.Draw()
      
   def drawDiffWith(self,other):
      self.diffTot = subtractHistograms(self.stacktot, other.stacktot)
      self.diffNonStackable = [ subtractHistograms(a,b) for a,b in zip(self.nonStackable, other.nonStackable) ]
      self.diffTot.Draw("hist,e,p")
      for i in self.diffNonStackable:
         i.Draw("hist,same")
      self.leg.Draw()
      
   def drawComparisonWith(self,other,ymax=+5.0,ymin=None):
      if ymin is None:
        ymin=-ymax
      import math
      comparison = lambda a,b : ( (a-b)/math.sqrt(a+b) if a+b != 0 else 0 )
      error = lambda a,b : 1 
     
      self.compTot = processHistograms(self.stacktot, other.stacktot, comparison, error)
      self.compTot.SetMaximum(ymax)
      self.compTot.SetMinimum(ymin)
      self.compTot.SetMarkerStyle(0)
      self.compTot.SetFillStyle(0)
      self.compTot.SetLineColor(ROOT.kBlack)
      self.compTot.SetLineWidth(1)
      self.compTot.Draw("hist,e")
      self.compNonStackable = []
      for i,j in zip(self.nonStackable, other.nonStackable):
         currentComparison = processHistograms(addHistograms(i,self.stacktot), addHistograms(j,other.stacktot), comparison, error)
         self.compNonStackable.append(currentComparison)
         currentComparison.Draw("hist,e,same")
      self.leg.Draw()
      

print "Opening root files..."
class Dataset:
   def __init__(self, shortName, filenameList, crossSecPb, colour, legendText=None, stackable=False):
      if not legendText:
         legendText = shortName
      self.stackable = stackable
      self.shortName=shortName
      self.legendText = legendText
      self.filenameList = filenameList
      self.rootfileList = []
      self.treeList = []
      self.colour = colour
      self.crossSecPb = crossSecPb
      self.crossSecFb = crossSecPb*1000.0
      self.numberOfEvents = 0
      for filename in filenameList:
          rootfile = ROOT.TFile(filename)
          self.rootfileList.append(rootfile)
          tree = rootfile.Get("Delphes")
          self.treeList.append(tree)
          self.numberOfEvents += tree.GetEntries()
      self.hists = {}
      print "I believe that there are "+str(self.numberOfEvents)+" events in "+str(filenameList)

   def project(self, shortHistName, texName, nbin, xmin, xmax, var, selection, desiredIntLumi_eventPerFb, legend):
      rootHistogramName = shortHistName+"_from_"+self.shortName
      self.hists[ shortHistName ] = ROOT.TH1F(rootHistogramName, texName, nbin, xmin, xmax)
      temporaryClone = self.hists[ shortHistName ].Clone()
      tempName = "temporaryHistogram"
      temporaryClone.SetName(tempName)
      #hist.SetDirectory(0)  NEVER DO THIS
      #hist.SetStats(0)
      #self.hists[ shortHistName ].SetFillStyle(3001) # Light hatching
      #self.hists[ shortHistName ].SetLineWidth(0)
      self.hists[ shortHistName ].SetLineColor(self.colour)
      if self.stackable:
         self.hists[ shortHistName ].SetFillColor(self.colour)
         self.hists[ shortHistName ].SetLineWidth(0)
      else:
         self.hists[ shortHistName ].SetLineWidth(1)
      if legend:
          legend.AddEntry(self.hists[ shortHistName ], self.legendText)
      scaling = desiredIntLumi_eventPerFb*self.crossSecFb/self.numberOfEvents
      print "RESCALING for "+self.shortName+" is ",scaling," since desiredIntLumi_eventPerFb*self.crossSecFb/self.numberOfEvents = ",desiredIntLumi_eventPerFb, " * ",self.crossSecFb," / ",self.numberOfEvents
      if scaling<1:
         scaledSelection = "("+selection+")*(("+randflat+"<"+str(scaling)+")?1.0:0.0)"
      else: 
         scaledSelection = "("+selection+")*("+str(scaling)+")"
      for tree in self.treeList:
         tree.Project(tempName, var, scaledSelection, "hist")
         self.hists[ shortHistName ].Add(temporaryClone)
          
      return self.hists[ shortHistName ]


def lesterProject(datasets, shortHistName, nbin, xmin, xmax, var, selection, desiredIntLumi_eventPerFb, texName=None):
    bundle = Bundle()

    if not texName:
       texName = selection+' ; '+var
    bundle.leg = ROOT.TLegend(0.5,0.7,0.95,0.90)
    bundle.stackable = []
    bundle.nonStackable = []
    bundle.stack = ROOT.THStack(
                      shortHistName+"_stack",
                      texName,
                        )
    for dataset in datasets:
       hist = dataset.project(shortHistName, texName, nbin, xmin, xmax, var, selection, desiredIntLumi_eventPerFb,bundle.leg)
       if dataset.stackable:
          bundle.stack.Add(hist)
          bundle.stackable.append(hist)
       else:
          bundle.nonStackable.append(hist)
          
    bundle.stack.SetMinimum(0.5)
    bundle.stacktot = bundle.stack.GetStack().Last().Clone()
    bundle.stacktot.SetMarkerStyle(20)
    bundle.stacktot.SetMarkerSize(0.7)
    bundle.stacktot.SetMarkerColor(ROOT.kBlack)
    return bundle

datasets = [
#      Dataset("singleW",[
#"lester_test_singlewleptonic/Events/run_01/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 100
#"lester_test_singlewleptonic/Events/run_02/tag_1_delphes_events.root", # High stats run - 160000 events. Seed 200
#],crossSecPb=1.933e4, colour=ROOT.kCyan,stackable=True), 

      Dataset("W_4h_up",[
"lester_test_singlewleptonic/Events/run_06/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 0.  metmin 400 metmax none
],crossSecPb=0.02849, colour=ROOT.kCyan,stackable=True), 
      Dataset("W_3h_4h",[
"lester_test_singlewleptonic/Events/run_03/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 0.  metmin 300 metmax 400
],crossSecPb=0.06487, colour=ROOT.kCyan+1,stackable=True), 
      Dataset("W_2h_3h",[
"lester_test_singlewleptonic/Events/run_04/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 0.  metmin 200 metmax 300
],crossSecPb=0.3459, colour=ROOT.kCyan-3,stackable=True), 
      Dataset("W_1h_2h",[
"lester_test_singlewleptonic/Events/run_05/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 0.  metmin 100 metmax 200
],crossSecPb=4.669, colour=ROOT.kCyan-6,stackable=True), 
      Dataset("W_70_1h",[
"lester_test_singlewleptonic/Events/run_08/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 0.  metmin 70 metmax 100
],crossSecPb=13.6, colour=ROOT.kCyan-8,stackable=True), 
      Dataset("W_50_70",[
"lester_test_singlewleptonic/Events/run_09/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 0.  metmin 50 metmax 70
],crossSecPb=67.15, colour=ROOT.kCyan-1,stackable=True), 
      Dataset("W_30_50",[
"lester_test_singlewleptonic/Events/run_10/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 0.  metmin 30 metmax 50
],crossSecPb=11290, colour=ROOT.kCyan+3,stackable=True), 
      Dataset("W_0_30",[
"lester_test_singlewleptonic/Events/run_10/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 0.  metmin 30 metmax 50
],crossSecPb=8033, colour=ROOT.kCyan+4,stackable=True), 
#      Dataset("W_0h_1h",[
#"lester_test_singlewleptonic/Events/run_07/tag_1_delphes_events.root", # Low  stats run -  10000 events. Seed 0.  metmin 000 metmax 100
#],crossSecPb=1.936e4, colour=ROOT.kCyan+4,stackable=True), 

 #Dataset("tw",["lester_test_tw_leptonicw/Events/run_01/tag_1_delphes_events.root"],crossSecPb=12.08, colour=ROOT.kYellow,stackable=True), # Low stats run
 Dataset("tw",["lester_test_tw_leptonicw/Events/run_02/tag_1_delphes_events.root"],crossSecPb=12.08, colour=ROOT.kYellow,stackable=True), # High stats run 
 Dataset("ttbar",[
         "ttbar_leptonic/Events/run_03/tag_1_delphes_events.root", # 1,000,000 events with different seed
         "ttbar_leptonic/Events/run_02/tag_1_delphes_events.root", # 1,000,000 events
         #"ttbar_leptonic/Events/run_01/tag_1_delphes_events.root", #   500,000 events # Don't use -- uses same seed as run_02
         ],crossSecPb=23.12, colour=ROOT.kRed, stackable=True),
 #Dataset("RPV",["lester_test_leptonemu_n1_top/Events/run_01/tag_1_delphes_events.root"], crossSecPb=0.02433*100, colour=ROOT.kGreen), # Note the *100 elevates the lambda'231 coupling to unity from 0.1.  This was a mixture of left and right smuons with masses around 320 GeV, and a neutralino around 90 GeV.  By the looks of later simulations (see below) this was 99.99% smuon right production.
 #Dataset("RPV_1_150_250_X",["lester_test_leptonemu_n1_top_1_150_250_X/Events/run_02/tag_3_delphes_events.root"], crossSecPb=2.238e-41, colour=ROOT.kMagenta), # numbers mean, respectively: lambda'231, neutralino mass (GeV), smuon right mass (GeV), smuon left mass (GeV).  X indicates infinity
 Dataset("RPV_1_150_X_250",[
       "lester_test_leptonemu_n1_top_1_150_X_250/Events/run_02/tag_3_delphes_events.root", # 50001 events -- zero seed
       "lester_test_leptonemu_n1_top_1_150_X_250/Events/run_03/tag_3_delphes_events.root", # 90000 events -- different seed
], crossSecPb=1.674, colour=ROOT.kMagenta+2), 
 #Dataset("RPV_1_90_300_X",["lester_test_leptonemu_n1_top_1_90_300_X/Events/run_02/tag_3_delphes_events.root"], crossSecPb=2.43e-41, colour=ROOT.kGreen), 
 Dataset("RPV_1_90_X_300",[
#"lester_test_leptonemu_n1_top_1_90_X_300/Events/run_02/tag_3_delphes_events.root", # Low stats - but don't use, as has same seed as run_04
"lester_test_leptonemu_n1_top_1_90_X_300/Events/run_04/tag_3_delphes_events.root", # High stats (200,000 events)
], crossSecPb=2.843, colour=ROOT.kGreen+2), # Smuon left
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

bundles = []

osdflep="(Muon_size==1&&Electron_size==1&&Muon.Charge[0]*Electron.Charge[0]<0)"
osdflepandhighptmiss = "("+osdflep+"&&MissingET_size==1&&MissingET.MET[0]>150)"

canvas = ROOT.TCanvas("x","Plots for "+str(desired_int_lumi_eventPerfb)+"/fb",1000,800)



if True:
   canvas.Divide(3,3)
   canvas.cd(1).SetLogy(0)
   bundles.append(lesterProject(datasets, "testHist-001", 25, 0, 800, "MissingET.MET[0]","MissingET_size==1&&Muon_size>0&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()
   
   
   canvas.cd(2).SetLogy(0)
   bundles.append(lesterProject(datasets, "testHist-002", 25, 0, 800, "MissingET.MET[0]","MissingET_size==1&&Muon_size>0&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()
   
   
   canvas.cd(3)
   bundles[-1].drawComparisonWith(bundles[-2],ymin=-5.0,ymax=25.0)
   canvas.Update()


   canvas.Print("Moo0-3.pdf")


canvas.Clear()


canvas.Divide(3,3)


if True:
   canvas.cd(1).SetLogy()
   bundles.append(lesterProject(datasets, "testHist001", 25, 0, 800, "MissingET.MET[0]","MissingET_size==1&&Muon_size>0&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()
   
   
   canvas.cd(2).SetLogy()
   bundles.append(lesterProject(datasets, "testHist002", 25, 0, 800, "MissingET.MET[0]","MissingET_size==1&&Muon_size>0&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()
   
   
   canvas.cd(3)
   bundles[-1].drawComparisonWith(bundles[-2],ymin=-5.0,ymax=25.0)
   canvas.Update()

   canvas.cd(4).SetLogy()
   bundles.append(lesterProject(datasets, "testHist004", 10, 0, 10, "Muon_size", "1==1", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()


   canvas.cd(5).SetLogy()
   bundles.append(lesterProject(datasets, "testHist005", 10, -2.5, 2.5, "Muon.Eta[0]",osdflep+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()
   
   canvas.cd(6).SetLogy()
   bundles.append(lesterProject(datasets, "testHist006", 10, -2.5, 2.5, "Muon.Eta[0]",osdflep+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()
   
   #canvas.cd(4).SetLogy()
   #bundles.append(lesterProject(datasets, "testHist004", 10, 0, 10, "Jet_size", "1==1", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   #bundles[-1].draw()
   #canvas.Update()
   
   canvas.cd(7).SetLogy()
   bundles.append(lesterProject(datasets, "testHist007", 40, 0, 400, "MissingET.MET[0]", "MissingET_size==1", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()
   
   
   canvas.cd(8).SetLogy()
   bundles.append(lesterProject(datasets, "testHist008", 10, -2.5, 2.5, "Muon.Eta[0]",osdflepandhighptmiss+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()
   
   canvas.cd(9).SetLogy()
   bundles.append(lesterProject(datasets, "testHist009", 10, -2.5, 2.5, "Muon.Eta[0]",osdflepandhighptmiss+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].draw()
   canvas.Update()

   canvas.Print("Moo1-3.pdf")


canvas.Clear()
canvas.Divide(3,3)


canvas.cd(1).SetLogy()
bundles.append(lesterProject(datasets, "testHist011", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
bundles[-1].draw()
canvas.Update()

canvas.cd(2).SetLogy()
bundles.append(lesterProject(datasets, "testHist012", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
bundles[-1].draw()
canvas.Update()

bit  = randflat+"<-1||"

canvas.cd(3).SetLogy()
bundles.append(lesterProject(datasets, "testHist013", 24, -0.1, 1.1, randflat, bit+bit+"1==1", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
bundles[-1].draw()
canvas.Update()

canvas.cd(4) #.SetLogy()
bundles.append(lesterProject(datasets, "testHist014", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
bundles[-1].draw()
canvas.Update()

canvas.cd(5) #.SetLogy()
bundles.append(lesterProject(datasets, "testHist015", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
bundles[-1].draw()
canvas.Update()

canvas.cd(6)
#bundles[-1].drawDiffWith(bundles[-2])
bundles[-1].drawComparisonWith(bundles[-2])
canvas.Update()

canvas.cd(7).SetLogy()
bundles.append(lesterProject(datasets, "testHist017", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
bundles[-1].draw()
canvas.Update()

canvas.cd(8).SetLogy()
bundles.append(lesterProject(datasets, "testHist018", 20, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
bundles[-1].draw()
canvas.Update()

canvas.cd(9)
#bundles[-1].drawDiffWith(bundles[-2])
bundles[-1].drawComparisonWith(bundles[-2])
canvas.Update()

canvas.Print("Moo2-3.pdf")

canvas.Clear()
canvas.Divide(3,3)


for i in range(1,10):
   canvas.cd(i)
   bundles.append(lesterProject(datasets, "testHist02"+str(i)+"a", 10+i*5, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]>0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles.append(lesterProject(datasets, "testHist02"+str(i)+"b", 10+i*5, 0, 800, "MissingET.MET[0]", "MissingET_size==1&&"+osdflep+"&&Muon.Charge[0]<0", desiredIntLumi_eventPerFb = desired_int_lumi_eventPerfb))
   bundles[-1].drawComparisonWith(bundles[-2])
   canvas.Update()

canvas.Print("Moo3-3.pdf")



import sys
#while True:
#  import time
#  time.sleep(1)

sys.exit(1)
