import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys


def make_lines(V, H):
	MAX = H.GetMaximum()*1.14
	L = []
	for i in range(1, len(V)-1):
		tempL = TLine(V[i], 0., V[i], MAX)
		tempL.SetLineColor(kBlue)
		L.append(tempL)
	return L	

def Optimize(x, H):
	Tot = H.Integral(H.FindBin(x[0]), H.FindBin(x[-1]))
	tol = 1./(len(x)-1)
	for i in range(len(x)-2):
		print " ==== "
		print i
		ConMet = False
		while not ConMet:
			print " try " + str(x[i+1])
			print "bins:" + str(H.FindBin(x[i])) + ", " + str(H.FindBin(x[i+1])-1)
			BinIntRel = H.Integral(H.FindBin(x[i]), H.FindBin(x[i+1])-1)/Tot
			print BinIntRel
			if BinIntRel < (tol + tol*0.06) and BinIntRel > (tol - tol*0.06):
				print "bingo"
				ConMet = True
			elif BinIntRel > (tol + tol*0.06):
				x[i+1] = x[i+1] - 2
			else:
				x[i+1] = x[i+1] +2

	

Bkgs =["GJetsmvaEVv3.root", "QCDmvaEVwithExtv3.root"]
#Bkgs =["QCDmvaEVwithExtv3.root"]

H_PT = TH1F("PT", ";jet p_{T} (GeV)", 350, 150, 1000)
H_PT.SetStats(0)
H_PT.SetLineColor(2)
H_PT.SetFillColor(2)
H_PT.SetFillStyle(3003)
for b in Bkgs:
	infile=ROOT.TFile(b)	
	print " "
	tree= infile.Get("Events")
	nent = tree.GetEntries();

	for i in range(tree.GetEntries()):
	
		tree.GetEntry(i)
		if(i % (1 * nent/100) == 0):
			sys.stdout.write("\r" + "="*int(20*i/nent) + "> " + str(round(100.*i/nent,0)) + "% done")
			sys.stdout.flush()

		puweight = tree.puWeight
		fbweight = tree.scale1fb
		weight = puweight*fbweight
		jmsd_8 = tree.AK8Puppijet0_msd
		
		PT = tree.AK8Puppijet0_pt
		if not PT > 0.: 
			continue
		if jmsd_8 <= 0: jmsd_8 = 0.01

		RHO = math.log(jmsd_8*jmsd_8/PT/PT)

		if RHO < -6.6 or RHO > -2.0: 
			continue

		if tree.neleLoose == 0 and tree.nmuLoose == 0 and tree.ntau==0 and tree.AK8Puppijet0_isTightVJet ==1:
			if tree.AK8Puppijet0_N2sdb1 > 0. and PT > 200. and PT < 1000.:
				H_PT.Fill(PT, weight)
	print " "

#H_PT.Scale(1/H_PT.Integral())
x = [200, 300, 400, 500, 600, 1000]
Optimize(x, H_PT)

L = make_lines(x, H_PT)

H_PT.GetYaxis().SetRangeUser(0,H_PT.GetMaximum()*1.14)

leg = TLegend(0.6,0.7,0.89,0.89)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(H_PT, "non-resonant MCs", "F")
leg.AddEntry(L[0], "#frac{1}{5} p_{T} bin edges", "L")

C = TCanvas("C", "", 800, 600)
C.cd()
H_PT.Draw("hist")
for l in L:
	l.Draw("same")
leg.Draw("same")
C.Print("OptimizedPtBins.png")
print "----Bin Edges: "
print x
print "------------------------------------------"
print "tot = " + str(H_PT.Integral(H_PT.FindBin(x[0]), H_PT.FindBin(x[-1])))
for i in range(len(x)-1):
	print " -=-"
	print str(x[i]) + ", " + str(x[i+1])
	print str(H_PT.Integral(H_PT.FindBin(x[i]), H_PT.FindBin(x[i+1])-1)) + " ~ " + str(H_PT.Integral(H_PT.FindBin(x[i]), H_PT.FindBin(x[i+1])-1) /H_PT.Integral(H_PT.FindBin(x[0]), H_PT.FindBin(x[-1])))
