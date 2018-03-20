# Fitter

import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys
from multiprocessing import Pool
import itertools

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

	#P = "([0] + [1]*y + [2]*y*y + [3]*y*y*y + ([4] + [5]*y + [6]*y*y + [7]*y*y*y)*x + ([8] + [9]*y + [10]*y*y + [11]*y*y*y)*x*x)"

	G = TF2("ralpha"+str(X)+str(Y), P, 200.,1000.,-6.5-2.)
	for n in range(numbering):
		G.SetParameter(n, 0.)
	return G # Makes a TF2 in 2d with X,Y order
	
def ComputeDDT(name, point, nPtBins, nRhoBins, H):
	DDT = TH2F(name, "", nRhoBins, -6.5, -2.0, nPtBins, 200, 1000)
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
	return DDT # Actual quantile evaluation

def FindAndSetMax(someset):
        maximum = 0.0
        for i in someset:
                i.SetStats(0)
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in someset:
                j.GetYaxis().SetRangeUser(0,maximum*1.35) # takes a set of histograms, sets their Range to show all of them

def MakeDDTMassComparison(DDT, Bkgs, name, skip):
	H = TH2F("H"+name, ";Jet Mass (GeV);N_{2}^{DDT}", 40, 0, 200, 40, -0.3,0.4)
	H.SetStats(0)
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
			if (T.AK8Puppijet0_N2sdb1 > 0. and M*M/PT/PT) > 0.0 and PT > 200 and PT < 1000 and M>0.:
				RHO = 2.*math.log(M/PT)
				if  RHO < -2.0 and RHO > -6.5:
					if isinstance(DDT, TH2F):
						cur_rho_index = DDT.GetXaxis().FindBin(RHO);
						cur_pt_index  = DDT.GetYaxis().FindBin(PT);
						N2DDT = T.AK8Puppijet0_N2sdb1 - DDT.GetBinContent(cur_rho_index,cur_pt_index);

					elif isinstance(DDT, TF2): N2DDT = T.AK8Puppijet0_N2sdb1 - DDT.Eval(RHO, PT)
					else: print "THIS ISN'T AN ACCEPTABLE DDT OBJECT"
					H.Fill(M,N2DDT,weight)
	return H
def MakeDDTMassComparison_Unpack(args):
	return MakeDDTMassComparison(*args)


def FillEst(DDT, Bkgs, name, skip):

	NonResPass = TH1F("NRP"+name, name+";AK8 Soft Drop Mass [GeV];events", 40, 0, 200)
	NonResFail = TH1F("NRF"+name, name+";AK8 Soft Drop Mass [GeV];events", 40, 0, 200)
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
			if (T.AK8Puppijet0_N2sdb1 > 0. and M*M/PT/PT) > 0.0 and PT > 200 and PT < 1000 and M>0.:
				RHO = 2.*math.log(M/PT)
				if  RHO < -2.0 and RHO > -6.5:
					if isinstance(DDT, TH2F):
						cur_rho_index = DDT.GetXaxis().FindBin(RHO);
						cur_pt_index  = DDT.GetYaxis().FindBin(PT);
						N2DDT = T.AK8Puppijet0_N2sdb1 - DDT.GetBinContent(cur_rho_index,cur_pt_index);

					elif isinstance(DDT, TF2): N2DDT = T.AK8Puppijet0_N2sdb1 - DDT.Eval(RHO, PT)
					else: print "THIS ISN'T AN ACCEPTABLE DDT OBJECT"

					
					if N2DDT < 0: #pass
						NonResPass.Fill(M,weight)
					else: #fail
						NonResFail.Fill(M,weight)

	NonResFail.Scale(5./95.)


	NonResPass.Sumw2()
	NonResPass.SetLineColor(1)
	NonResPass.SetFillColor(0)
	NonResPass.SetMarkerStyle(4)

	NonResFail.SetStats(0)
	NonResFail.SetLineWidth(1)
	NonResFail.GetXaxis().SetTitleSize(0.045)
	NonResFail.GetYaxis().SetTitleSize(0.045)

	return [NonResPass, NonResFail] # Estimates (and plots) pass from fail (which it also plots)

def TryWithOrder(O):
	i = O[0]
	j = O[1]
	Label = "O("+str(i)+","+str(j)+")"
	G = MakeSuperPoly(i,j)
	G.SetParameter(0, 0.2)
	G.SetLineColor(1)

	ddt.Fit("ralpha"+str(i)+str(j))

	SavePoly = TFile("SmoothDDT_O"+str(i)+str(j)+".root", "recreate")
	SavePoly.cd()
	G.Write()
	SavePoly.Write()
	SavePoly.Close()

	C = TCanvas("C"+str(i)+str(j), "", 600, 600)
	C.cd()
	ddt.Draw("SURF1")
	C.Print("DDTPolyFit_O"+str(i)+str(j)+".root")
	C.Print("DDTPolyFit_O"+str(i)+str(j)+".png")

	MC = FillEst(G, Bkgs, "MC_O"+str(i)+str(j), 1)
	DATA = FillEst(G, Data, "Data_O"+str(i)+str(j), 15)
	CompHmc = MakeDDTMassComparison(G, Bkgs, "MCMvDDT_O"+str(i)+str(j), 1)
	CompHdata = MakeDDTMassComparison(G, Data, "DataMvDDT_O"+str(i)+str(j), 15)

	CmcMd = TCanvas("CmcMd"+Label, "", 600, 600)
	CmcMd.cd()
	CompHmc.SetTitle("MC"+Label)
	CompHmc.Draw("COL")
	CmcMd.Print("MC_MassVDDT_O"+str(i)+str(j)+".png")
	CdMd = TCanvas("CdMd"+Label, "", 600, 600)
	CdMd.cd()
	CompHdata.SetTitle("Data"+Label)
	CompHdata.Draw("COL")
	CdMd.Print("Data_MassVDDT_O"+str(i)+str(j)+".png")

	return [MC[0], MC[1], DATA[0], DATA[1], Label]

if __name__ == '__main__':

	print "--- Start ---"

	RB = 19
	PB = 19

	H3 = TH3F("H3", "", RB, -6.5, -2.0, PB, 200, 1000, 750, 0., 0.75)
	H3.SetStats(0)
	H3T = TH3F("H3T", "", RB, -6.5, -2.0, PB, 200, 1000, 500, 0., 1.)
	H3T.SetStats(0)
	#Bkgs =["GJetsmvaEVv3.root", "QCDmvaEVwithExtv3.root"]
	Bkgs =["GJetsmvaEVv3.root"]
	Data = ["SinglePhotonmvaEV12av3.root"]


	o = [1,2]
	Orders = list(itertools.product(o, o))

	print "--- Env Vars Set ---"
	print "--- Will consider the following pT/rho polynomial ordres:"
	print Orders

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
					if RHO < -2.0 and RHO > -6.5 and T.AK8Puppijet0_msd > 0.:
						if T.AK8Puppijet0_N2sdb1 > 0.: 
							H3.Fill(RHO, PT, T.AK8Puppijet0_N2sdb1, weight)
						if T.AK8Puppijet0_tau21 > 0. :
							H3T.Fill(RHO, PT, T.AK8Puppijet0_tau21, weight) # #Creates the DDT mapping TH3 which we will run on.

	ddt = ComputeDDT("N2DDT", 0.05, PB, RB, H3)



	SaveDDT = TFile("MapDDT.root", "recreate")
	SaveDDT.cd()
	ddt.Write()
	SaveDDT.Write()
	SaveDDT.Close()


	del H3

	p = Pool(4)
	MvDDTresults = p.map(MakeDDTMassComparison_Unpack, [(ddt, Bkgs, "MCMvDDT", 1),(ddt, Data, "DataMvDDT", 15)])

	CmcMd = TCanvas("CmcMd", "", 600, 600)
	CmcMd.cd()
	MvDDTresults[0].SetTitle("MC")
	MvDDTresults[0].Draw("COL")
	CmcMd.Print("MC_MassVDDT_NOFIT.png")
	CdMd = TCanvas("CdMd", "", 600, 600)
	CdMd.cd()
	MvDDTresults[1].SetTitle("Data")
	MvDDTresults[1].Draw("COL")
	CdMd.Print("Data_MassVDDT_NOFIT.png")

	print "--- DDT Map created ---"

	EstsMCp = []
	EstsDatap = []
	EstsMCf = []
	EstsDataf = []

	print "--- Starting Fits ---"


	MCH = FillEst(ddt, Bkgs, "MC_nofit", 1)
	DATAH = FillEst(ddt, Data, "Data_nofit", 15)
	Results = p.map(TryWithOrder, Orders)

	Results.append([MCH[0], MCH[1], DATAH[0], DATAH[1], "NoFit"])

	print "--- Fits Done ---"

	leg = TLegend(0.475,0.25,0.89,0.89)
	leg.SetLineColor(0)
	leg.SetFillColor(0)

	colors = [1,2,4,6,8,9,12,28,32,42,46,50,7,38,44,39,3,5]

	for i in range(len(Results)):
		leg.AddEntry(Results[i][0], "pass "+Results[i][4] + " fit", "PL")
		leg.AddEntry(Results[i][1], "5/95 fail "+Results[i][4] + " fit", "L")

		Results[i][0].SetMarkerColor(colors[i])
		Results[i][0].SetLineColor(colors[i])
		Results[i][1].SetMarkerColor(colors[i])
		Results[i][1].SetLineColor(colors[i])
		Results[i][2].SetMarkerColor(colors[i])
		Results[i][2].SetLineColor(colors[i])
		Results[i][3].SetMarkerColor(colors[i])
		Results[i][3].SetLineColor(colors[i])

		EstsMCp.append(Results[i][0])
		EstsMCf.append(Results[i][1])
		EstsDatap.append(Results[i][2])
		EstsDataf.append(Results[i][3])


	print "--- Plotting ---"

	FindAndSetMax(EstsMCp+EstsMCf)
	FindAndSetMax(EstsDatap+EstsDataf)

	EstsMCf[0].SetTitle("MC")
	EstsDataf[0].SetTitle("Data")

	Cmc = TCanvas("Cmc", "", 800, 800)
	Cmc.cd()
	EstsMCf[0].Draw("hist")
	Cdata = TCanvas("Cdata", "", 800, 800)
	Cdata.cd()
	EstsDataf[0].Draw("hist")

	for i in range(len(EstsMCp)):
		Cmc.cd()
		EstsMCp[i].Draw("sameE0")
		EstsMCf[i].Draw("samehist")
		Cdata.cd()
		EstsDatap[i].Draw("sameE0")
		EstsDataf[i].Draw("samehist")

		tempL = TLegend(0.5,0.5,0.89,0.89)
		tempL.SetLineColor(0)
		tempL.SetFillColor(0)
		tempL.SetHeader("Poly "+Results[i][4])

		CTempMC = TCanvas("CTempMC"+str(i), "", 500, 500)
		CTempMC.cd()
		EstsMCp[i].Draw("E0")
		EstsMCf[i].Draw("samehist")
		CTempMC.Print("MC_FitComparison"+Results[i][4]+".png")
		CTempData = TCanvas("CTempData"+str(i), "", 500, 500)
		CTempData.cd()
		EstsDatap[i].Draw("E0")
		EstsDataf[i].Draw("samehist")
		CTempData.Print("Data_FitComparison"+Results[i][4]+".png")



	Cmc.cd()
	leg.Draw("same")
	Cdata.cd()
	leg.Draw("same")
	Cmc.Print("MC_FitComparison.png")
	Cdata.Print("Data_FitComparison.png")
