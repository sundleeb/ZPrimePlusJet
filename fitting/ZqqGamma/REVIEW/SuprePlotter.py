import ROOT
from ROOT import *
import sys
import math
from array import array
import functools
from multiprocessing import Pool
import scipy

lumi = 29.7

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

def quickplot(File, tree, plot, var, Cut, Weight):
        temp = plot.Clone("temp")
        chain = ROOT.TChain(tree)
        chain.Add(File)
        chain.Draw(var+">>"+"temp", "("+Weight+")*("+Cut+")", "goff")
        plot.Add(temp)

class VAR:
	def __init__(self, V, B, T):
		self.B = B
		self.T = T
		self.V = V

def FindAndSetMax(someset):
        maximum = 0.0
        for i in someset:
                i.SetStats(0)
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in someset:
                j.GetYaxis().SetRangeUser(0.9,maximum*2.)

def MakeNicePlotOf(Var, S):

	cuts = "AK4Puppijet0_pt>200.&AK8Puppijet0_pt>200.&vpho0_pt>200."
	DH = TH1F("DH", "", Var.B[0], Var.B[1], Var.B[2])
	DH.Sumw2()
	DH.SetLineColor(1)
	DH.SetFillColor(0)
	DH.SetMarkerColor(1)
	DH.SetMarkerStyle(20)
	quickplot("SinglePhotonmvaEV12av3_q.root", "Events", DH, Var.V, cuts, "1.0")

	GH = TH1F("GH", "", Var.B[0], Var.B[1], Var.B[2])
	GH.SetFillColor(kMagenta-10)
	quickplot("GJetsmvaEVv3withExt.root", "Events", GH, Var.V, cuts, "scale1fb*puWeight*kfactorNLO*"+str(lumi))

	QH = TH1F("QH", "", Var.B[0], Var.B[1], Var.B[2])
	QH.SetFillColor(kBlue-7)
	quickplot("QCDmvaEVwithExtv3_q.root", "Events", QH, Var.V, cuts, "scale1fb*puWeight*kfactor*"+str(lumi))
	quickplot("QCDmvaEVwithExtv3_q.root", "Events", GH, Var.V, cuts, "scale1fb*puWeight*kfactor*"+str(lumi))

	Draxis = TH1F("Draxis", ";"+Var.T+";events", Var.B[0], Var.B[1], Var.B[2])
	Draxis.SetStats(0)
	Praxis = TH1F("Praxis", ";;Data/MC", Var.B[0], Var.B[1], Var.B[2])
	Praxis.GetYaxis().SetRangeUser(0.0,2.0)
	Praxis.SetStats(0)

	Praxis.GetXaxis().SetNdivisions(0)
	Praxis.GetYaxis().SetNdivisions(6)
	Praxis.GetYaxis().CenterTitle(True)
	Praxis.GetYaxis().SetLabelSize(85/22*Praxis.GetYaxis().GetLabelSize())
	Praxis.GetYaxis().SetTitleSize(4.4*Praxis.GetYaxis().GetTitleSize())
	Praxis.GetYaxis().SetTitleOffset(0.175)
	Draxis.GetXaxis().SetTitleOffset(1.3)

	FindAndSetMax([Draxis, DH, GH])

	C = TCanvas("C", "", 800, 800)
	plot = TPad("plot", "The pad 80% of the height",0,0.15,1,1)
	pull = TPad("pull", "The pad 20% of the height",0,0,1.0,0.15)
	plot.Draw()
	pull.Draw()
	plot.cd()
	plot.SetLogy()
	Draxis.Draw()
	GH.Draw("histsame")
	QH.Draw("histsame")
	DH.Draw("e0same")
	pull.cd()
	Praxis.Draw()
	AddCMSLumi(plot, lumi, "Preliminary")
	C.Print(S+".png")

JetPt = VAR("AK4Puppijet0_pt", [50,0,800], "AK8 Jet p_{T} (GeV)")
MakeNicePlotOf(JetPt, "test")
