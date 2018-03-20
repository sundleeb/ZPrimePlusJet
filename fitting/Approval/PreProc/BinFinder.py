import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys
import MicroHelperFunctions
from MicroHelperFunctions import *

def make_lines(V, H):
	MAX = H.GetMaximum()*1.14
	L = []
	for i in range(1, len(V)-1):
		tempL = TLine(V[i], 0., V[i], MAX)
		tempL.SetLineColor(kBlue)
		L.append(tempL)
	return L

H_PT = TH1F("PT", ";jet p_{T} (GeV)", 400, 200, 1000)
H_PT.SetStats(0)
H_PT.SetLineColor(2)
H_PT.SetFillColor(2)
H_PT.SetFillStyle(3003)
quickplot("gjets.root", "tree", H_PT, "pT", "rho<-2.&rho>-7.5", "weight*puW*TrigW*kfNLO")

H_PT.Scale(1/H_PT.Integral())
if sys.argv[1] == "2": x = [200, 450, 1000]
if sys.argv[1] == "3": x = [200, 400, 700, 1000]
if sys.argv[1] == "4": x = [200, 300, 400, 600, 1000]
if sys.argv[1] == "5": x = [200, 300, 400, 500, 600, 1000]
if sys.argv[1] == "6": x = [200, 300, 400, 500, 600, 800, 1000]

Optimize(x, H_PT)

L = make_lines(x, H_PT)

H_PT.GetYaxis().SetRangeUser(0,H_PT.GetMaximum()*1.14)

leg = TLegend(0.6,0.7,0.89,0.89)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(H_PT, "non-resonant MCs", "F")
leg.AddEntry(L[0], "#frac{1}{"+sys.argv[1]+"} p_{T} bin edges", "L")

C = TCanvas("C", "", 800, 600)
C.cd()
H_PT.Draw("hist")
for l in L:
	l.Draw("same")
leg.Draw("same")
C.Print("Optimized"+sys.argv[1]+"PtBins.png")
print "----Bin Edges: "
print x
print "------------------------------------------"
print "tot = " + str(H_PT.Integral(H_PT.FindBin(x[0]), H_PT.FindBin(x[-1])))
for i in range(len(x)-1):
	print " -=-"
	print str(x[i]) + ", " + str(x[i+1])
	print str(H_PT.Integral(H_PT.FindBin(x[i]), H_PT.FindBin(x[i+1])-1)) + " ~ " + str(H_PT.Integral(H_PT.FindBin(x[i]), H_PT.FindBin(x[i+1])-1) /H_PT.Integral(H_PT.FindBin(x[0]), H_PT.FindBin(x[-1])))
