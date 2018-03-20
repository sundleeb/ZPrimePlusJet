import ROOT,sys,math,array,os
from ROOT import *
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array

def FindAndSetMax(*someset):
        maximum = 0.0
        for i in someset:
                i.SetStats(0)
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in someset:
                j.GetYaxis().SetRangeUser(0,maximum*1.35)

#gROOT.SetBatch(kTRUE)

def plotCategory(ML,i, W):
	histograms_old = ML.Get("shapes_prefit/ch%i_%s_cat%i/%s" % (i,W,i,'tqq'))
	rags = ML.Get("norm_prefit")
	rrv_old = RooRealVar(rags.find("ch%i_%s_cat%i/%s" % (i,W,i,'tqq')))
	if rrv_old.getVal() > 0.: histograms_old.Scale(rrv_old.getVal()/histograms_old.Integral())

	histograms_new = ML.Get("shapes_fit_b/ch%i_%s_cat%i/%s" % (i,W,i,'tqq'))
	rags = ML.Get("norm_fit_b")
	rrv_new = RooRealVar(rags.find("ch%i_%s_cat%i/%s" % (i,W,i,'tqq')))
	if rrv_new.getVal() > 0.: histograms_new.Scale(rrv_new.getVal()/histograms_new.Integral())

	EO = 0.
	EN = 0.
	for i in range(1,histograms_new.GetNbinsX()+1):
		EO += histograms_old.GetBinError(i)
		EN += histograms_new.GetBinError(i)
	return (histograms_new, histograms_old, EO, EN)

	# after:




def DoAllPlots(which):
	print "START"
	ML = TFile("mlfitTTSB.root");
	oldHists = []
	newHists = []
	EO = 0.
	EN = 0.
	for i in range(4): 
		(tN, tO, tEO, tEN) = plotCategory(ML, i+1, which)
		oldHists.append(tO)
		newHists.append(tN)
		EO = TMath.Sqrt((EO*EO) + (tEO*tEO))
		EN = TMath.Sqrt((EN*EN) + (tEN*tEN))
	
	for i in range(1,4):
		oldHists[0].Add(oldHists[i])
		newHists[0].Add(newHists[i])
	oldHists[0].SetLineColor(kBlue)
	oldHists[0].SetLineWidth(2)
	oldHists[0].SetTitle(";Soft Drop Mass (GeV);Events")
	newHists[0].SetTitle(";Soft Drop Mass (GeV);Events")
	oldHists[0].SetStats(0)
	newHists[0].SetLineColor(kRed)
	newHists[0].SetLineWidth(2)
	L = TLegend(0.65,0.675,0.89,0.89)
	L.SetLineColor(0)
	L.SetFillColor(0)
	L.AddEntry(oldHists[0], "prefit t#bar{t}", "L")
	L.AddEntry(newHists[0], "postfit t#bar{t}", "L")
	C = TCanvas("C", "", 700, 700)
	C.cd()
	
	FindAndSetMax(newHists[0], oldHists[0])
	if oldHists[0].Integral() > newHists[0].Integral():
		newHists[0].SetFillColor(10)
		newHists[0].SetFillStyle(1001)
		oldHists[0].SetFillColor(kBlue)
		oldHists[0].SetFillStyle(3002)
		L.AddEntry(oldHists[0], "change in t#bar{t}", "F")
		oldHists[0].Draw("hist")
		newHists[0].Draw("histsame")
	else:
		oldHists[0].SetFillColor(10)
		oldHists[0].SetFillStyle(1001)
		newHists[0].SetFillColor(kBlue)
		newHists[0].SetFillStyle(3002)
		L.AddEntry(newHists[0], "change in t#bar{t}", "F")
		newHists[0].Draw("hist")
		oldHists[0].Draw("histsame")
		


	gPad.RedrawAxis()
	L.Draw("same")
	C.Print("TTCOMP_"+which+".png")

	print "old " + which + " " + str(oldHists[0].Integral()) + " +/- " + str(EO)
	print "new " + which + " " + str(newHists[0].Integral()) + " +/- " + str(EN)


if __name__ == '__main__':
	import optparse
	from optparse import OptionParser
	parser = OptionParser()

	parser.add_option('--W', '--which', metavar='WHC', type='string', dest='which', default="prefit")
	(Options, args) = parser.parse_args()
	directory = "OUTPUTMARCSTYLE"
	if not os.path.exists(directory):
   		os.makedirs(directory)
	DoAllPlots("pass")
	DoAllPlots("fail")

