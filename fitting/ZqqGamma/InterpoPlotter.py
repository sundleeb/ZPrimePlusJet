#
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

if __name__ == '__main__':
	print "-----------------------------------------------"
	BlankP = TH1F("BLANKP", ";AK8 Soft Drop Mass (GeV);a.u.", 30, 0, 150)
	BlankF = TH1F("BLANKF", ";AK8 Soft Drop Mass (GeV);a.u.", 30, 0, 150)
	BlankP.SetStats(0)
	BlankF.SetStats(0)
	Phists = []
	Fhists = []
	colors = []
	FTrue = TFile("ZGammaFermiSRTemplates.root")
	PhistsReal = []
	FhistsReal = []
	for m in [10,25,50,75,100,125]:
		pT = FTrue.Get("zqq"+str(m)+"pass_msd")
		fT = FTrue.Get("zqq"+str(m)+"fail_msd")
		pT.SetLineColor(1)
		pT.SetLineWidth(2)
		fT.SetLineColor(1)
		fT.SetLineWidth(2)
		#pT.Scale(1./pT.Integral())
		#fT.Scale(1./fT.Integral())
		PhistsReal.append(pT)
		FhistsReal.append(fT)
	FInter = TFile("InterpolatedMassPoints.root")
	Ms = [15.,20.,30.,35.,40.,45.,55.,60.0,65.,70.,80.,85.,90.0,95.,105.,110.0,115.,120.]
	i = 0
	for m in Ms:
		# new color
		
		ci = 1000 + i
		C = (float(i)/float(len(Ms)))
		print "index = " + str(ci) + ", color modifier = " + str(C)
		color = TColor(ci, 1.0-C, 0.0, 0.0+C);
		colors.append(color)
		tempP = FInter.Get("zqq"+str(int(m))+"_pass_cat1")
		tempP.Add(FInter.Get("zqq"+str(int(m))+"_pass_cat2"))
		tempP.Add(FInter.Get("zqq"+str(int(m))+"_pass_cat3"))
		tempP.Add(FInter.Get("zqq"+str(int(m))+"_pass_cat4"))
		tempP.Add(FInter.Get("zqq"+str(int(m))+"_pass_cat5"))
		#tempP.Scale(1./tempP.Integral())
		tempF = FInter.Get("zqq"+str(int(m))+"_fail_cat1")
		tempF.Add(FInter.Get("zqq"+str(int(m))+"_fail_cat2"))
		tempF.Add(FInter.Get("zqq"+str(int(m))+"_fail_cat3"))
		tempF.Add(FInter.Get("zqq"+str(int(m))+"_fail_cat4"))
		tempF.Add(FInter.Get("zqq"+str(int(m))+"_fail_cat5"))
		#tempF.Scale(1./tempF.Integral())

		tempP.SetLineColor(ci)
		tempF.SetLineColor(ci)
		
		Phists.append(tempP)
		Fhists.append(tempF)
		i += 1

	FindAndSetMax(Phists)
	FindAndSetMax(PhistsReal)
	FindAndSetMax(Phists[0], BlankP, PhistsReal[0])
	FindAndSetMax(Fhists)
	FindAndSetMax(FhistsReal)
	FindAndSetMax(Fhists[0], BlankF, FhistsReal[0])

	Cp = TCanvas("Cp", "", 800, 800)
	Cp.cd()
	BlankP.Draw()
	for h in PhistsReal:
		h.Draw("histsame")
	Cp.Print("InterpoGifStorage/BlankPass.gif")
	w = 0
	for m in Phists:
		m.SetLineWidth(4)
		m.Draw("histsame")
		Cp.Print("InterpoGifStorage/BlankPass"+str(w)+".gif")
		m.SetLineWidth(1)
		Cp.Update()
		w += 1

	Cf = TCanvas("Cf", "", 800, 800)
	Cf.cd()
	BlankF.Draw()
	for h in FhistsReal:
		h.Draw("histsame")
	Cf.Print("InterpoGifStorage/BlankFail.gif")
	w = 0
	for m in Fhists:
		m.SetLineWidth(4)
		m.Draw("histsame")
		Cf.Print("InterpoGifStorage/BlankFail"+str(w)+".gif")
		m.SetLineWidth(1)
		Cf.Update()
		w += 1




