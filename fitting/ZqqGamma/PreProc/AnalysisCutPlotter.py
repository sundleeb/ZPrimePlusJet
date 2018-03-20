import ROOT
from ROOT import *
import sys
import math
import scipy
import functools
from array import array
sys.path.append('/home/marc/code/python/')
import PlottingFunctions
import RootHelperFunctions


def FindAndSetMax(someset):
        maximum = 0.0
        for i in someset:
                i.SetStats(0)
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in someset:
                j.GetYaxis().SetRangeUser(0,maximum*1.35)
        return maximum*1.35


def DoIsoPlot(F, title):
	IsoVals = TH2F("IsoPt"+title, ";Jet p_{T} (GeV); #Delta R_{j, #gamma}", 40, 0, 800, 15, 0., 3.0)
	IsoVals.SetStats(0)
	IsoValsClosest = TH2F("IsoPtC"+title, title+ ";Jet p_{T} (GeV); #Delta R_{j, #gamma}", 40, 0, 800, 15, 0., 3.0)
	IsoValsClosest.SetStats(0)

	File = TFile(F)
	T = File.Get("Events")
	n = T.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
		T.GetEntry(j)
		if T.vpho0_pt > 200 and T.AK8Puppijet0_isTightVJet ==1 and T.nmuLoose == 0 and T.ntau==0 and T.vpho0_mva > 0.2 and T.AK8Puppijet0_pt > 200. and T.AK8Puppijet0_msd > 0.:
			PT = T.AK8Puppijet0_pt
			RHO = 2.*math.log(T.AK8Puppijet0_msd/T.AK8Puppijet0_pt)
			if not (RHO < -2.0 and RHO > -7.0): continue
			Jet = TLorentzVector()
			Jet.SetPtEtaPhiM(PT, T.AK8Puppijet0_eta, T.AK8Puppijet0_phi, T.AK8Puppijet0_mass)
			Photon = TLorentzVector()
			Photon.SetPtEtaPhiM(T.vpho0_pt, T.vpho0_eta, T.vpho0_phi, 0.0)
			if not (Jet.DeltaR(Photon) > 2.2): continue
			considercsvs = []
			j1 = TLorentzVector()
			j2 = TLorentzVector()
			j3 = TLorentzVector()
			j4 = TLorentzVector()
			j1.SetPtEtaPhiM(T.AK4Puppijet0_pt, T.AK4Puppijet0_eta, T.AK4Puppijet0_phi, T.AK4Puppijet0_mass)
			j2.SetPtEtaPhiM(T.AK4Puppijet1_pt, T.AK4Puppijet1_eta, T.AK4Puppijet1_phi, T.AK4Puppijet1_mass)
			j3.SetPtEtaPhiM(T.AK4Puppijet2_pt, T.AK4Puppijet2_eta, T.AK4Puppijet2_phi, T.AK4Puppijet2_mass)
			j4.SetPtEtaPhiM(T.AK4Puppijet3_pt, T.AK4Puppijet3_eta, T.AK4Puppijet3_phi, T.AK4Puppijet3_mass)

			notisolated = False
			closest = 999.
			closept = -1.

			for J in [j1,j2,j3,j4]:
				if J.Pt() > 15:
					dR = J.DeltaR(Photon)
					IsoVals.Fill(J.Pt(), dR)
					if dR < closest:
						closest = dR
						closept = J.Pt()
			IsoValsClosest.Fill(closept, closest)

	C1 = TCanvas("C1", "", 500, 400)
	C1.cd()
	IsoVals.Draw("col")
	C1.Print("IsoPtAll_"+title+".png")

	C2 = TCanvas("C2", "", 500, 400)
	C2.cd()
	IsoValsClosest.Draw("col")
	C2.Print("IsoPtClose_"+title+".png")
					
def GetMaxCSVHist(T):
	H = TH1F("H"+T.GetName(), ";maximum ak4 jet csv;%", 50, 0, 1)
	H.SetStats(0)
	n = T.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
		T.GetEntry(j)
		if T.vpho0_pt > 200 and T.AK8Puppijet0_isTightVJet ==1 and T.nmuLoose == 0 and T.ntau==0 and T.vpho0_mva > 0.2 and T.AK8Puppijet0_pt > 200. and T.AK8Puppijet0_msd > 0.:
			PT = T.AK8Puppijet0_pt
			RHO = 2.*math.log(T.AK8Puppijet0_msd/T.AK8Puppijet0_pt)
			if not (RHO < -2.0 and RHO > -7.0): continue
			Jet = TLorentzVector()
			Jet.SetPtEtaPhiM(PT, T.AK8Puppijet0_eta, T.AK8Puppijet0_phi, T.AK8Puppijet0_mass)
			Photon = TLorentzVector()
			Photon.SetPtEtaPhiM(T.vpho0_pt, T.vpho0_eta, T.vpho0_phi, 0.0)
			if not (Jet.DeltaR(Photon) > 2.2): continue
			considercsvs = []
			considercsvs.append(T.AK4Puppijet0_csv)
			considercsvs.append(T.AK4Puppijet1_csv)
			considercsvs.append(T.AK4Puppijet2_csv)
			considercsvs.append(T.AK4Puppijet3_csv)

			if len(considercsvs) > 0:
				maxbinevent = max(considercsvs)
			else: maxbinevent = 0.
			H.Fill(maxbinevent)
	H.Scale(100/H.Integral())
	H.SetLineWidth(2)
	return H

def GetMETHist(T, C):
	H = TH1F("H"+T.GetName(), ";pfmet (GeV);%", 25, 0, 250)
	H.SetStats(0)
	n = T.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
		T.GetEntry(j)
		if T.vpho0_pt > 200 and T.nmuLoose == 0 and T.ntau==0 and T.vpho0_mva > 0.2 and T.AK8Puppijet0_pt > 200. and T.AK8Puppijet0_msd > 0.:
			PT = T.AK8Puppijet0_pt
			RHO = 2.*math.log(T.AK8Puppijet0_msd/T.AK8Puppijet0_pt)
			if not (RHO < -2.0 and RHO > -7.0): continue
			Jet = TLorentzVector()
			Jet.SetPtEtaPhiM(PT, T.AK8Puppijet0_eta, T.AK8Puppijet0_phi, T.AK8Puppijet0_mass)
			Photon = TLorentzVector()
			Photon.SetPtEtaPhiM(T.vpho0_pt, T.vpho0_eta, T.vpho0_phi, 0.0)
			H.Fill(T.pfmet)
	H.Scale(100/H.Integral())
	H.SetLineWidth(2)
	H.SetLineStyle(2)
	H.SetLineColor(C)
	return H
	
def MaxCSVComparator():
	FG = TFile("GJetsmvaEVv3withExt.root")
	TG = FG.Get("Events")
	HG  = GetMaxCSVHist(TG)
	FT = TFile("TTmvaEVv3_q.root")
	TT = FT.Get("Events")
	HT  = GetMaxCSVHist(TT)
	FS = TFile("VectorDiJet1Gammamva_M50_q.root")
	TS = FS.Get("Events")
	HS  = GetMaxCSVHist(TS)
	HG.SetLineColor(kBlue)
	HT.SetLineColor(kRed)
	HS.SetLineColor(kBlack)
	C = TCanvas("sfads", "", 600, 600)
	C.cd()
	HT.Draw("hist")
	HG.Draw("histsame")
	HS.Draw("histsame")
	C.Print("CSVComparator.png")

def METComparator():
	FS10 = TFile("VectorDiJet1Gammamva_M10_q.root")
	FS25 = TFile("VectorDiJet1Gammamva_M25_q.root")
	FS50 = TFile("VectorDiJet1Gammamva_M50_q.root")
	FS75 = TFile("VectorDiJet1Gammamva_M75_q.root")
	FS100 = TFile("VectorDiJet1Gammamva_M100_q.root")
	TS10 = FS10.Get("Events")
	TS25 = FS25.Get("Events")
	TS50 = FS50.Get("Events")
	TS75 = FS75.Get("Events")
	TS100 = FS100.Get("Events")
	HS10  = GetMETHist(TS10, kBlack)
	HS25  = GetMETHist(TS25, kBlue)
	HS50  = GetMETHist(TS50, kGreen)
	HS75  = GetMETHist(TS75, kRed)
	HS100  = GetMETHist(TS100, kViolet)
	FTT = TFile("TTmvaEVv3_q.root")
	TTT = FTT.Get("Events")
	HTT  = GetMETHist(TTT, kMagenta)
	HTT.SetLineStyle(1)
	FindAndSetMax([HS10, HS25, HS50, HS75, HS100, HTT])
	L = TLegend(0.55, 0.6, 0.89,0.89)
	L.SetFillColor(0)
	L.SetLineColor(0)
	L.AddEntry(HTT, "t#bar{t}", "L")
	L.AddEntry(HS10, "Z'_{10 GeV}", "L")
	L.AddEntry(HS25, "Z'_{25 GeV}", "L")
	L.AddEntry(HS50, "Z'_{50 GeV}", "L")
	L.AddEntry(HS75, "Z'_{75 GeV}", "L")
	L.AddEntry(HS100, "Z'_{100 GeV}", "L")

	C = TCanvas("sfads", "", 600, 600)
	C.cd()
	HS10.Draw("hist")
	HS25.Draw("histsame")
	HS50.Draw("histsame")
	HS75.Draw("histsame")
	HS100.Draw("histsame")
	HTT.Draw("histsame")
	L.Draw("same")
	C.Print("METComparator.png")

#MaxCSVComparator()
METComparator()

#DoIsoPlot("VectorDiJet1Gammamva_M10_q.root", "sig10")
#DoIsoPlot("VectorDiJet1Gammamva_M50_q.root", "sig50")
#DoIsoPlot("VectorDiJet1Gammamva_M100_q.root", "sig100")
#DoIsoPlot("TTmvaEVv3_q.root", "ttbar")
