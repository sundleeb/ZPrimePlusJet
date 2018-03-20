import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys
import MicroHelperFunctions
from MicroHelperFunctions import *

def GetPunzi(S, B, N):
	nB = 0
	nS = 0
	for n in range(S.GetXaxis().GetNbins()+1):
		if S.GetBinContent(n) < S.Integral()*0.01:
			continue
		nB = nB + B.GetBinContent(n)
		nS = nS + S.GetBinContent(n)
	del B
	sigEff = nS/N
	print sigEff
	return sigEff/(5/2 + math.sqrt(nB))
	

SR10_DATA = TFile("DATASETS/ZGp10.root")
St10P = SR10_DATA.Get("zqq10_pass_msd")
St25P = SR10_DATA.Get("zqq25_pass_msd")
St50P = SR10_DATA.Get("zqq50_pass_msd")
St75P = SR10_DATA.Get("zqq75_pass_msd")
St100P = SR10_DATA.Get("zqq100_pass_msd")
St125P = SR10_DATA.Get("zqq125_pass_msd")
QCDP10 = SR10_DATA.Get("qcd_pass_msd")
ttP10 = SR10_DATA.Get("tqq_pass_msd")
zP10 = SR10_DATA.Get("wqq_pass_msd")
wP10 = SR10_DATA.Get("zqq_pass_msd")
#QCDP10.Add(ttP10)
#QCDP10.Add(zP10)
#QCDP10.Add(wP10)

S = [St10P, St25P, St50P, St75P, St100P, St125P]
m = [10.,25.,50.,75.,100.,125.]
xs = [0.6927, 0.5254, 0.3312, 0.2868, 0.2473, 0.2345]
N = [101413, 80878, 103297,99824,96104,96628]
MVAScale = [0.9252, 0.90715, 0.9029, 0.8832, 0.87195, 0.8651]
P = []
for i in range(len(S)):
	S[i].Scale(15*N[i]/(35900*xs[i]))
	p = GetPunzi(S[i], QCDP10, N[i])
	print p
	P.append(p)

G = TGraph(6, scipy.array(m), scipy.array(P))
G.SetLineColor(kBlue)

Ax = TH1F("Ax", ";Z' mass (GeV);Punzi Significanc", 28, 0, 140)
Ax.SetStats(0)
Ax.GetYaxis().SetRangeUser(0.,max(P)*1.25)
T = TCanvas()
T.cd()
Ax.Draw()
G.Draw("Csame")



