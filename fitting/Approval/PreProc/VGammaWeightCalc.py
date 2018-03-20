##

import ROOT
from ROOT import *
def FindAndSetMax(someset):
        maximum = 0.0
        for i in someset:
                i.SetStats(0)
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in someset:
                j.GetYaxis().SetRangeUser(0,maximum*1.35)

def PlotSR(name, F, W, C):
	H = TH1F("H"+name, ";AK_{8} Soft Drop Mass (GeV);events (1fb^{-1})", 80, 0, 200)
	H.SetStats(0)
	H.SetLineColor(C)
	H.SetLineWidth(3)
	chain = ROOT.TChain("tree")
	chain.Add(F)
	chain.Draw("SDM*TheaW>>"+"H"+name, "("+W+")*(SR>0&PhoEta>-2.1&PhoEta<2.1&rho<-2.&rho>-7.5)", "goff")
	return H

HW = PlotSR("W", "WQQ.root", "puW*TrigW*kfNLO*1000.*1.22/1507795.125", kBlue)
HZ = PlotSR("Z", "ZQQ.root", "puW*TrigW*kfNLO*1000.*0.569/398015.4276", kRed)
print str(35.9*HZ.Integral())
HT = PlotSR("T", "ttbar.root", "puW*TrigW*kfNLO*weight", kBlack)
FindAndSetMax([HW, HZ, HT])
L = TLegend(0.65,0.7,0.89,0.85)
L.SetFillColor(0)
L.SetLineColor(0)
L.AddEntry(HT, "ttbar", "L")
L.AddEntry(HW, "W#gamma", "L")
L.AddEntry(HZ, "Z#gamma", "L")
C = TCanvas()
C.cd()
HW.Draw("hist")
HZ.Draw("histsame")
HT.Draw("histsame")
L.Draw()


