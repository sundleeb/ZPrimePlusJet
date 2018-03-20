import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys
import MicroHelperFunctions
from MicroHelperFunctions import *

def GetSRB(s, b, C):
	countS = 0
	countB = 0
	binsCon = 0
	for n in range(b.GetXaxis().GetNbins()+1):
		b.SetBinContent(n, math.sqrt(b.GetBinContent(n)))
		if s.GetBinContent(n) < s.Integral()*C:
			continue
		binsCon += 1
		countB = countB + b.GetBinContent(n)
		countS = countS + s.GetBinContent(n)
	if binsCon == 0:
		return GetSRB(s, b, C*0.75)
	return countS/countB

# Mass Hists:
SR_DATA = TFile("SR15_DDT5_MCF.root")
S10P = SR_DATA.Get("zqq10_pass_msd")
S25P = SR_DATA.Get("zqq25_pass_msd")
S50P = SR_DATA.Get("zqq50_pass_msd")
S75P = SR_DATA.Get("zqq75_pass_msd")
S100P = SR_DATA.Get("zqq100_pass_msd")
S125P = SR_DATA.Get("zqq125_pass_msd")
S10F = SR_DATA.Get("zqq10_fail_msd")
S25F = SR_DATA.Get("zqq25_fail_msd")
S50F = SR_DATA.Get("zqq50_fail_msd")
S75F = SR_DATA.Get("zqq75_fail_msd")
S100F = SR_DATA.Get("zqq100_fail_msd")
S125F = SR_DATA.Get("zqq125_fail_msd")

QCDP5 = SR_DATA.Get("qcd_pass_msd")

SR10_DATA = TFile("SR15_DDT10_MCF.root")
St10P = SR10_DATA.Get("zqq10_pass_msd")
St25P = SR10_DATA.Get("zqq25_pass_msd")
St50P = SR10_DATA.Get("zqq50_pass_msd")
St75P = SR10_DATA.Get("zqq75_pass_msd")
St100P = SR10_DATA.Get("zqq100_pass_msd")
St125P = SR10_DATA.Get("zqq125_pass_msd")


QCDP10 = SR10_DATA.Get("qcd_pass_msd")

QCDP10.Rebin(1)
QCDP5.Rebin(1)

#S10P.Add(S10F)
#S25P.Add(S25F)
#S50P.Add(S50F)
#S75P.Add(S75F)
#S100P.Add(S100F)
#S125P.Add(S125F)

AllPSigs = [S10P, S25P, S50P, S75P, S100P, S125P]
AllFSigs = [S10F, S25P, S50F, S75F, S100F, S125F]
AllPSigsT = [St10P, St25P, St50P, St75P, St100P, St125P]
for S in AllPSigs:
	S.Rebin(1)
	S.SetLineWidth(2)
	S.SetLineColor(kBlue)
	#S.Scale(1./S.Integral())
	S.SetStats(0)
L = []
M = FindAndSetMax(AllPSigs)*1.2
for i in [10., 25., 50., 75., 100., 125.]: 
	l = TLine(i, 0., i, M)
	l.SetLineColor(kRed)
	L.append(l)
	
#T = TCanvas()
#T.cd()
#S10P.Draw("hist")
#for s in AllPSigs:
#	s.Draw("histsame")
#for l in L:
#	l.Draw("same")
#T.Print("CorrectedSignalSDMasses.png")

m = [10.,25.,50.,75.,100.,125.]
xs = [0.6927, 0.5254, 0.3312, 0.2868, 0.2473, 0.2345]
N = [101413, 80878, 103297,99824,96104,96628]

Acc = []
AccT = []
AccPF = []
SRB5 = []
SRB10 = []
for i in range(len(xs)):
	a = AllPSigs[i].Integral()/((35800./15.)*xs[i]) * 100.
	at = AllPSigsT[i].Integral()/((35800./15.)*xs[i]) * 100.
	f = (AllPSigs[i].Integral()+AllFSigs[i].Integral())/((35800./15.)*xs[i]) * 100.
	srb5 = GetSRB(AllPSigs[i], QCDP5, 0.1)
	srb10 = GetSRB(AllPSigsT[i], QCDP10, 0.1)
	Acc.append(a)
	AccT.append(at)
	AccPF.append(f)
	SRB5.append(srb5)
	SRB10.append(srb10)
G = TGraph(6, scipy.array(m), scipy.array(Acc))
G.SetLineColor(kBlue)
G2 = TGraph(6, scipy.array(m), scipy.array(AccPF))
G3 = TGraph(6, scipy.array(m), scipy.array(AccT))
G3.SetLineColor(kRed)
G4 = TGraph(6, scipy.array(m), scipy.array(SRB5))
G4.SetLineColor(kBlue)
G4.SetLineStyle(2)
G5 = TGraph(6, scipy.array(m), scipy.array(SRB10))
G5.SetLineStyle(2)
G5.SetLineColor(kRed)

L = TLegend(0.6,0.6,0.89,0.89)
L.SetLineColor(0)
L.SetFillColor(0)
L.AddEntry(G, "acc. #times eff., 5%^{DDT}", "L")
L.AddEntry(G3, "acc. #times eff., 10%^{DDT}", "L")
#L.AddEntry(G4, "s/#sqrt{b} (a.u.), 5%^{DDT}", "L")
#L.AddEntry(G5, "s/#sqrt{b} (a.u.), 10%^{DDT}", "L")

Ax = TH1F("Ax", ";Z' mass (GeV);acc. #times eff. (%)", 28, 0, 140)
Ax.SetStats(0)
Ax.GetYaxis().SetRangeUser(0.,3.5)
T2 = TCanvas()
T2.cd()
Ax.Draw()
G3.Draw("Csame")
G.Draw("Csame")
#G4.Draw("Csame")
#G5.Draw("Csame")
L.Draw("same")
T2.Print("SignalEfficiencies.png")








