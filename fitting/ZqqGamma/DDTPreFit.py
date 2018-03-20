# Fitter

import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys

RB = 4
PB = 4

def ComputeDDT(name, point, nPtBins, nRhoBins, H):
	DDT = TH2F(name, "", nRhoBins, -6.0, -2.0, nPtBins, 200, 800)
	DDT.SetStats(0)
	nXb = H.GetXaxis().GetNbins()
	nYb = H.GetYaxis().GetNbins()
	for x in range(nXb):
		for y in range(nYb):
			proj = H.ProjectionZ("H3"+str(x)+str(y),x+1,x+1,y+1,y+1)
			print str(x+1) + "," + str(y+1) + ":    "+ str(proj.Integral())
			p = array('d', [point])
			q = array('d', [0.0]*len(p))
			proj.GetQuantiles(len(p), q, p)
			DDT.SetBinContent( x+1, y+1, q[0] );
	return DDT

def FindAndSetMax(someset):
        maximum = 0.0
        for i in someset:
                i.SetStats(0)
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in someset:
                j.GetYaxis().SetRangeUser(0,maximum*1.35)

def FillEst(DDT, Bkgs, name, skip):

	NonResPass = TH1F("NRP", "", 40, 0, 200)
	NonResFail = TH1F("NRF", "", 40, 0, 200)
	for B in Bkgs:
		F = TFile(B)
		T = F.Get("Events")
		n = T.GetEntries()
		for j in range(0, n): # Here is where we loop over all events.
			if (j) % skip != 0: continue
			T.GetEntry(j)
			weight = T.puWeight*T.scale1fb*T.kfactor*T.kfactorNLO
			PT = T.AK8Puppijet0_pt
			M = T.AK8Puppijet0_msd
			if (T.AK8Puppijet0_N2sdb1 > 0. and M*M/PT/PT) > 0.0 and PT > 200 and PT < 800 and M>0.:
				RHO = 2.*math.log(M/PT)

				N2DDT = T.AK8Puppijet0_N2sdb1 - DDT.Eval(RHO, PT)
				if  RHO < -2.0 and RHO > -6.0:
					if N2DDT < 0: #pass
						NonResPass.Fill(M,weight)
					else: #fail
						NonResFail.Fill(M,weight)

	NonResFail.Scale(5./95.)

	FindAndSetMax([NonResPass, NonResFail])

	leg = TLegend(0.475,0.385,0.85,0.89)
	leg.SetHeader("QCD & #gamma + jets MCs (DDT~5%):")
	leg.AddEntry(NonResPass, "events N_{2}^{DDT} < 0", "PL")
	leg.AddEntry(NonResFail, "#frac{5}{95} #times events N_{2}^{DDT} > 0", "L")
	leg.SetLineColor(0)
	leg.SetFillColor(0)

	NonResPass.Sumw2()
	NonResPass.SetLineColor(1)
	NonResPass.SetFillColor(0)
	NonResPass.SetMarkerColor(1)
	NonResPass.SetMarkerStyle(20)

	NonResFail.SetStats(0)
	NonResFail.SetLineColor(kGreen)
	NonResFail.SetLineWidth(2)
	NonResFail.GetXaxis().SetTitle("AK8 Soft Drop Mass [GeV]")
	NonResFail.GetYaxis().SetTitle("Events")
	NonResFail.GetXaxis().SetTitleSize(0.045)
	NonResFail.GetYaxis().SetTitleSize(0.045)


	C = TCanvas("Cfitcheck"+name, "", 800, 600)
	C.cd()
	CMSLABL = TLatex()
	CMSLABL.SetNDC()
	CMSLABL.SetTextSize(0.045)
	PRELABL = TLatex()
	PRELABL.SetNDC()
	PRELABL.SetTextSize(0.04)
	THILABL = TLatex()
	THILABL.SetNDC()
	THILABL.SetTextSize(0.045)
	NonResFail.Draw("hist")
	NonResPass.Draw("E0 same")
	leg.Draw("same")
	C.Print("DDT_SanityCheck_"+name+"_FIT.png")


H3 = TH3F("H3", "", RB, -6.0, -2.0, PB, 200, 800, 750, 0., 0.75)
H3.SetStats(0)
H3T = TH3F("H3T", "", RB, -6.0, -2.0, PB, 200, 800, 500, 0., 1.)
H3T.SetStats(0)
Bkgs =["GJetsmvaEVv3.root", "QCDmvaEVwithExtv3.root"]
Data = ["SinglePhotonmvaEV12av3.root"]
for B in Bkgs:
	F = TFile(B)
	T = F.Get("Events")
	n = T.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
		T.GetEntry(j)
		if T.neleLoose == 0 and T.nmuLoose == 0 and T.ntau==0  and T.AK8Puppijet0_isTightVJet ==1:
			weight = T.puWeight*T.scale1fb*T.kfactor*T.kfactorNLO
			PT = T.AK8Puppijet0_pt
			if PT > 200.:
			    preRHO = T.AK8Puppijet0_msd*T.AK8Puppijet0_msd/T.AK8Puppijet0_pt/T.AK8Puppijet0_pt
			    if preRHO > 0.:
				RHO = math.log(preRHO)
				if RHO < -2.0 and RHO > -6.0 and T.AK8Puppijet0_msd > 0.:
					if T.AK8Puppijet0_N2sdb1 > 0.: 
						H3.Fill(RHO, PT, T.AK8Puppijet0_N2sdb1, weight)
					if T.AK8Puppijet0_tau21 > 0. :
						H3T.Fill(RHO, PT, T.AK8Puppijet0_tau21, weight)

ddt = ComputeDDT("N2DDT", 0.05, PB, RB, H3)

P = "([0] + [1]*y + [2]*y*y + [3]*y*y*y + ([4] + [5]*y + [6]*y*y + [7]*y*y*y)*x + ([8] + [9]*y + [10]*y*y + [11]*y*y*y)*x*x)"
G = TF2("ralpha", P, 200.,800.,-6.-2.)
G.SetParameter(0, 0.2)
G.SetParameter(1, 0.)
G.SetParameter(2, 0.)
G.SetParameter(3, 0.)
G.SetParameter(4, 0.)
G.SetParameter(5, 0.)
G.SetParameter(6, 0.)
G.SetParameter(7, 0.)
G.SetParameter(8, 0.)
G.SetParameter(9, 0.)
G.SetParameter(10, 0.)
G.SetParameter(11, 0.)

G.SetLineColor(1)

ddt.Fit("ralpha")

C = TCanvas("C", "", 600, 600)
C.cd()
ddt.Draw("SURF1")
C.Print("DDTPolyFit.root")

FillEst(G, Bkgs, "MC", 1)
FillEst(G, Data, "Data", 15)



