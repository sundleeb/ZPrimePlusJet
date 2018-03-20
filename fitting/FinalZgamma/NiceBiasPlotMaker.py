import ROOT
from ROOT import *
from optparse import OptionParser
import random
import scipy
import sys
import glob, os

def GetHistFromCanv(C):
	F = TFile(C)
	c = C.split(".")
	Cv = F.Get(c[0])
	H = Cv.GetPrimitive("h")
	return H

Band = TH2F("Band", ";Signal Mass (GeV);#Delta#mu/#sigma", 1, 0., 150., 10, -3., 3.)
Band.SetStats(0)

x = []
y = []
ex = []
ey = []
for file in glob.glob("pull_*"+sys.argv[1]+"*.root"):
	print(file)
	S = file.split("zqq")
	SS = S[1].split(".")[0]
	h = GetHistFromCanv(file)
	h.Fit("gaus")
	F = h.GetFunction("gaus")
	m = F.GetParameter(1)
	s = F.GetParameter(2)
	print m, "   ", s
	x.append(float(SS))
	ex.append(0.)
	y.append(m)
	ey.append(s)

G = TGraphErrors(len(x), scipy.array(x), scipy.array(y), scipy.array(ex), scipy.array(ey))
G.SetMarkerStyle(21)

B = TBox(0,-0.5,150,0.5)
B.SetFillColor(41)
B.SetLineColor(0)
L = TLine(0,0,150,0)
L.SetLineColor(1)
L.SetLineStyle(2)
L.SetLineWidth(2)
L1 = TLine(0,1,150,1)
L1.SetLineColor(1)
L1.SetLineStyle(3)
L1.SetLineWidth(1)
L1m = TLine(0,-1,150,-1)
L1m.SetLineColor(1)
L1m.SetLineStyle(3)
L1m.SetLineWidth(1)
L2 = TLine(0,2,150,2)
L2.SetLineColor(1)
L2.SetLineStyle(3)
L2.SetLineWidth(1)
L2m = TLine(0,-2,150,-2)
L2m.SetLineColor(1)
L2m.SetLineStyle(3)
L2m.SetLineWidth(1)

C = TCanvas()
C.cd()
Band.Draw()
B.Draw("same")
L.Draw("same")
L1.Draw("same")
L1m.Draw("same")
L2.Draw("same")
L2m.Draw("same")
G.Draw("sameP")

#H = []
#for file in glob.glob("pull_*"+sys.argv[1]+"GeV.root"):
#	print(file)
#	H.append(GetHistFromCanv(file))
#for i in range(1,len(H)):
#	H[0].Add(H[i])
#
#C = TCanvas()
#C.cd()
#H[0].Draw("e0")
#C.Print("NicePullPlot"+sys.argv[1]+"GeV.png")