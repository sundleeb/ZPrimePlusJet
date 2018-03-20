#

import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys


def MakeSuperPoly(X,Y):

        P = ""
        numbering = 0
        for x in range(X+1):
                ypart = "["+str(numbering)+"]"
                numbering += 1
                for y in range(1,Y+1):
                        ypart += " + ["+str(numbering)+"]"
                        for i in range(y):
                                ypart += "*y"
                        numbering += 1
                if x == 0:
                        string = ypart
                else:
                        string = " + ("+ypart+")"
                for i in range(x):
                        string += "*x"
                P += string


        print P

        G = TF2("ralpha"+str(X)+str(Y), P, 200.,1000.,-6.5-2.)
        for n in range(numbering):
                G.SetParameter(n, 0.)
        return G # Makes a TF2 in 2d with X,Y order

def FindAndSetMax(someset):
        maximum = 0.0
        for i in someset:
                i.SetStats(0)
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in someset:
                j.GetYaxis().SetRangeUser(0,maximum*1.35) # takes a set of histograms, sets their Range to show all of them

def createHists(DDT,tag,filename):
	print "processing " + filename

	mass = TH1F("M"+tag, ";Jet Mass (Gev)", 80, 0, 200)
	mass.SetStats(0)

	HN2 = TH1F("N2"+tag, ";Two Pronged Variable", 50, -0.5, 0.5)
	HN2DDT = TH1F("N2DDT"+tag, ";Two Pronged Variable", 50, -0.5, 0.5)


	infile=ROOT.TFile(filename+".root")	

	tree= infile.Get("Events")
	nent = tree.GetEntries();

	for i in range(tree.GetEntries()):

		tree.GetEntry(i)
		if(i % (1 * nent/100) == 0):
			sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done")
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

		if RHO < -6.5 or RHO > -2.0: 
			continue
		if not (tree.AK8Puppijet0_N2sdb1 > 0. and PT > 200. and PT < 1000.):
			continue
                if isinstance(DDT, TH2F):
                        cur_rho_index = DDT.GetXaxis().FindBin(RHO);
                        cur_pt_index  = DDT.GetYaxis().FindBin(PT);
                        N2DDT = tree.AK8Puppijet0_N2sdb1 - DDT.GetBinContent(cur_rho_index,cur_pt_index);
                elif isinstance(DDT, TF2): 
			N2DDT = tree.AK8Puppijet0_N2sdb1 - DDT.Eval(RHO, PT)
                else: print "THIS ISN'T AN ACCEPTABLE DDT OBJECT"

		HN2.Fill(tree.AK8Puppijet0_N2sdb1 - 0.25, weight)
		HN2DDT.Fill(N2DDT, weight)

		if N2DDT < 0: mass.Fill(jmsd_8, weight)
	
	HN2.Scale(1/HN2.Integral())
	HN2.SetLineColor(kBlack)
	HN2.SetLineWidth(2)
	HN2DDT.Scale(1/HN2DDT.Integral())
	HN2DDT.SetLineColor(kRed)
	HN2DDT.SetLineWidth(2)

	leg = TLegend(0.7,0.7,0.89,0.89)
	leg.SetLineColor(0)
	leg.SetFillColor(0)
	leg.AddEntry(HN2, "N_{2} (adjusted)", "L")
	leg.AddEntry(HN2DDT, "N_{2}^{DDT}", "L")

	FindAndSetMax([HN2, HN2DDT])

	C = TCanvas("Cnvnddt"+tag)
	C.cd()
	HN2.Draw("hist")
	HN2DDT.Draw("histsame")
	leg.Draw("same")
	C.Print("TwoProgComp"+tag+".png")

	C2 = TCanvas("C2"+tag)
	C2.cd()
	mass.Draw("hist")
	C2.Print("TwoProgMass"+tag+".png")


f_h2ddt = TFile("MapDDT.root");
trans_h2ddt = f_h2ddt.Get("N2DDT");

G = MakeSuperPoly(2,1)
G.SetParameter(0, 0.2)
G.SetLineColor(1)
trans_h2ddt.Fit("ralpha21")

#createHists(G,'data_obs','SinglePhotonmvaEV12av3')
createHists(G,'nonres','NonResBkg')
createHists(G,'top','TTmvaEVv3_1000pb_weighted')
createHists(G,'sig','VectorDiJet1Gammamva_M50')
