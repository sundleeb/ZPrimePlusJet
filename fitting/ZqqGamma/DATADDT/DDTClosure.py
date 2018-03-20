import ROOT
from ROOT import *
import sys
import math
from array import array
import functools
from multiprocessing import Pool
import scipy
import PlottingFunctions
from PlottingFunctions import *
def AddCMSLumi(pad, fb, extra):

	cmsText     = "CMS";
	cmsTextFont   = 61  

	extraText   = extra
	extraTextFont = 52 

	lumiTextSize     = 0.4
	lumiTextOffset   = 0.15

	cmsTextSize      = 0.5
	cmsTextOffset    = 0.15


	H = pad.GetWh()
	W = pad.GetWw()
	l = pad.GetLeftMargin()
	t = pad.GetTopMargin()
	r = pad.GetRightMargin()
	b = pad.GetBottomMargin()
	e = 0.025

	pad.cd()

	lumiText = str(fb)+" fb^{-1} (13 TeV)"

	latex = TLatex()
	latex.SetNDC()
	latex.SetTextAngle(0)
	latex.SetTextColor(kBlack)	
	
	extraTextSize = 0.76*cmsTextSize
	
	latex.SetTextFont(42)
	latex.SetTextAlign(31) 
	latex.SetTextSize(lumiTextSize*t)	

	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)

	pad.cd()

	latex.SetTextFont(cmsTextFont)
	latex.SetTextSize(cmsTextSize*t)
	latex.SetTextAlign(11)
	latex.DrawLatex(l, 1-t+cmsTextOffset*t, cmsText)
	latex.SetTextFont(extraTextFont)
	latex.SetTextAlign(11)
	latex.SetTextSize(extraTextSize*t)
	latex.DrawLatex(l+0.11, 1-t+cmsTextOffset*t, extraText)
 

	pad.Update()


if __name__ == '__main__':
	F = TFile("DDTmaps.root")
	B = F.Get("N2DDT")
	B1 = F.Get("N2DDT1")
	B2 = F.Get("N2DDT2")
	B3 = F.Get("N2DDT3")
	S = F.Get("SmooDDT")
	S1 = F.Get("SmooDDT1")
	S2 = F.Get("SmooDDT2")
	S3 = F.Get("SmooDDT3")

	C = TCanvas("C", "", 1095, 795)
#	C.Divide(2,1)
	C.cd(1)
	S.Draw("colz")
	AddCMSLumi(gPad, 36.6, "Simulation Preliminary")
	C.Print("SmooDDT.pdf")

	Mbp = TH1F("Mbp", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mbf = TH1F("Mbf", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Msp = TH1F("Msp", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Msf = TH1F("Msf", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)

	Mb1p = TH1F("Mb1p", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mb2p = TH1F("Mb2p", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mb3p = TH1F("Mb3p", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mb1f = TH1F("Mb1f", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mb2f = TH1F("Mb2f", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mb3f = TH1F("Mb3f", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mf1p = TH1F("Mf1p", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mf2p = TH1F("Mf2p", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mf3p = TH1F("Mf3p", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mf1f = TH1F("Mf1f", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mf2f = TH1F("Mf2f", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)
	Mf3f = TH1F("Mf3f", ";Jet Soft Drop Mass (GeV)", 20, 0, 200)

	F2 = TFile("gjets.root")
	T = F2.Get("tree")
	n = T.GetEntries()
	for j in range(0, n):# Here is where we loop over all events.
		T.GetEntry(j)
		B_rho_index = B.GetXaxis().FindBin(T.rho)
		B_pt_index  = B.GetYaxis().FindBin(T.pT)
		S_rho_index = S.GetXaxis().FindBin(T.rho)
		S_pt_index  = S.GetYaxis().FindBin(T.pT)
		BDDT = T.N2 - B.GetBinContent(B_rho_index,B_pt_index)
		SDDT = T.N2 - S.GetBinContent(S_rho_index,S_pt_index)
		BDDT1 = T.N2 - B1.GetBinContent(B_rho_index,B_pt_index)
		SDDT1 = T.N2 - S1.GetBinContent(S_rho_index,S_pt_index)
		BDDT2 = T.N2 - B2.GetBinContent(B_rho_index,B_pt_index)
		SDDT2 = T.N2 - S2.GetBinContent(S_rho_index,S_pt_index)
		BDDT3 = T.N2 - B3.GetBinContent(B_rho_index,B_pt_index)
		SDDT3 = T.N2 - S3.GetBinContent(S_rho_index,S_pt_index)
		if BDDT < 0:
			Mbp.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO)
		else:
			Mbf.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO*5./95.)
		if BDDT1 < 0:
			Mb1p.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO)
		else:
			Mb1f.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO*5./95.)
		if BDDT2 < 0:
			Mb2p.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO)
		else:
			Mb2f.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO*5./95.)
		if BDDT3 < 0:
			Mb3p.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO)
		else:
			Mb3f.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO*5./95.)
		if SDDT < 0:
			Msp.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO)
		else:
			Msf.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO*5./95.)
		if SDDT1 < 0:
			Mf1p.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO)
		else:
			Mf1f.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO*5./95.)
		if SDDT2 < 0:
			Mf2p.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO)
		else:
			Mf2f.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO*5./95.)
		if SDDT3 < 0:
			Mf3p.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO)
		else:
			Mf3f.Fill(T.SDM, T.weight*T.puW*T.TrigW*T.kfNLO*5./95.)

	
	PlottingFunctions.GoodPlotFormat(Mbp, "e0markers", kBlack, 20)
	PlottingFunctions.GoodPlotFormat(Mbf, "thinline", kGray, 20)
	PlottingFunctions.GoodPlotFormat(Msp, "e0markers", kBlack, 20)
	PlottingFunctions.GoodPlotFormat(Msf, "thinline", kGray, 20)

	PlottingFunctions.GoodPlotFormat(Mb1p, "e0markers", kAzure, 20)
	PlottingFunctions.GoodPlotFormat(Mb2p, "e0markers", kPink, 20)
	PlottingFunctions.GoodPlotFormat(Mb3p, "e0markers", kSpring, 20)
	PlottingFunctions.GoodPlotFormat(Mb1f, "thinline", kAzure)
	PlottingFunctions.GoodPlotFormat(Mb2f, "thinline", kPink)
	PlottingFunctions.GoodPlotFormat(Mb3f, "thinline", kSpring)
	PlottingFunctions.GoodPlotFormat(Mf1p, "e0markers", kAzure, 20)
	PlottingFunctions.GoodPlotFormat(Mf2p, "e0markers", kPink, 20)
	PlottingFunctions.GoodPlotFormat(Mf3p, "e0markers", kSpring, 20)
	PlottingFunctions.GoodPlotFormat(Mf1f, "thinline", kAzure)
	PlottingFunctions.GoodPlotFormat(Mf2f, "thinline", kPink)
	PlottingFunctions.GoodPlotFormat(Mf3f, "thinline", kSpring)



	FindAndSetMax(Mb1p, Mb1f, Mb2p, Mb2f, Mb3p, Mb3f)
	FindAndSetMax(Mf1p, Mf1f, Mf2p, Mf2f, Mf3p, Mf3f)

	FindAndSetMax(Mbp, Mbf, Msp, Msf)

	L = TLegend(0.43,0.63,0.89,0.89)
	L.SetFillColor(0)
	L.SetLineColor(0)
	L.AddEntry(Mb1p, "Closure test using 1^{st} 3^{d} of MC events", "PL")
	L.AddEntry(Mb2p, "Closure test using 2^{nd} 3^{d} of MC events", "PL")
	L.AddEntry(Mb3p, "Closure test using 3^{d} 3^{d} of MC events", "PL")

#	C1 = TCanvas("C1", "", 800, 800)
#	C1.Divide(2,1)
#	C1.cd(1)
#	Mf1f.Draw("hist")
#	Mf2f.Draw("samehist")
#	Mf3f.Draw("samehist")
#	Mf1p.Draw("samee0")
#	Mf2p.Draw("samee0")
#	Mf3p.Draw("samee0")
#	L.Draw("same")
#	AddCMSLumi(gPad, 36.6, "Simulation Preliminary")
#	C1.Print("BadCloseSmoo.pdf")
#	C1.cd(2)
#	Mf1f.Draw("hist")
#	Mf2f.Draw("samehist")
#	Mf3f.Draw("samehist")
#	Mf1p.Draw("samee0")
#	Mf2p.Draw("samee0")
#	Mf3p.Draw("samee0")

	L = TLegend(0.5,0.89,0.5,0.89)
	L.SetLineColor(0)
	L.SetFillColor(0)
	L.AddEntry(Msp, "Events Passing the N_{2}^{DDT} Cut", "PL")
	L.AddEntry(Msf, "#frac{5}{95} #times Events Failing the N_{2}^{DDT} Cut", "L")
	C2 = TCanvas("C2", "", 800, 800)
#	C2.Divide(2,1)
	C2.cd(1)
	Msf.Draw("hist")
	Msp.Draw("samee0")
	L.Draw("same")
	AddCMSLumi(gPad, 36.6, "Simulation Preliminary")
	C2.Print("TrivialClose.pdf")
	#C2.cd(2)
	#Msf.Draw("hist")
	#Msp.Draw("samee0")
