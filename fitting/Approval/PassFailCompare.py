#
import ROOT,sys,math,array,os
from ROOT import *
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array
import PlotterFunctions
from PlotterFunctions import *


def MakeBSCOMP(F):
	print " --- --- --- --- --- --- --- --- --- --- --- "
	ML = TFile(F+"/MLF.root")
	B = TFile(F+"/base.root")
	rB = TFile(F+"/ralphabase.root")
	FolderPartsA = F.split("--")
	FolderPartsB = FolderPartsA[0].split("_p")
	PathName = FolderPartsB[0].replace("OUTPUT_", "")
	ORIGINALS = TFile("PreProc/DATASETS/"+PathName+".root")
	nptbins = ORIGINALS.Get("data_obs_pass").GetYaxis().GetNbins()
	print "      # pT bins: " + str(nptbins)
	ML_B = RooFitResult(ML.Get("fit_b"))
	ML_S = RooFitResult(ML.Get("fit_s"))
	shapes = ['wqq','zqq','tqq','qcd','zqq50']
	print "\n\n\n"
	D = []
	Ss = []
	B = []
	for j in range(nptbins):
		i = j+1
		histograms_pass_s = [];
		histograms_pass_b = [];
		for num,shapename in enumerate(shapes):
			print "looking in fit_s/ch%i_fail_cat%i/%s" % (i,i,shapename)
			histograms_pass_s.append( ML.Get("shapes_fit_s/ch%i_pass_cat%i/%s" % (i,i,shapename)))
			rags = ML.Get("norm_fit_s")
			rrv_pass = RooRealVar(rags.find("ch%i_pass_cat%i/%s" % (i,i,shapename)))
			curnorm_pass = rrv_pass.getVal()
			if curnorm_pass > 0.: histograms_pass_s[num].Scale(curnorm_pass/histograms_pass_s[num].Integral())

			print "looking in fit_b/ch%i_fail_cat%i/%s" % (i,i,shapename)
			histograms_pass_b.append( ML.Get("shapes_fit_b/ch%i_pass_cat%i/%s" % (i,i,shapename)))
			rags = ML.Get("norm_fit_b")
			rrv_pass = RooRealVar(rags.find("ch%i_pass_cat%i/%s" % (i,i,shapename)))
			curnorm_pass = rrv_pass.getVal()
			if curnorm_pass > 0.: histograms_pass_b[num].Scale(curnorm_pass/histograms_pass_b[num].Integral())
		data = convertAsymGraph(ML.Get("shapes_fit_b/ch%i_pass_cat%i/data" % (i,i)),histograms_pass_b[0],"shapes_fit_b/ch%i_pass_cat%i/data")

		###
		data.SetMarkerColor(1)
		data.SetLineColor(1)
		data.SetMarkerStyle(10)
		D.append(data)
		###

		AllB = histograms_pass_b[0].Clone("AllB")
		AllB.SetLineColor(kBlue)
		AllS = histograms_pass_s[0].Clone("AllS")
		AllS.SetLineColor(kRed)
		AllS.SetFillColor(10)
		AllS.SetFillStyle(1001)
		AllS.SetStats(0)
		AllSs = histograms_pass_s[0].Clone("AllS")
		AllSs.SetLineColor(kRed)
		AllSs.SetFillColor(kRed-4)
		AllSs.SetFillStyle(3002)
		AllSs.SetStats(0)
		AllSs.SetTitle(";Jet Soft Drop Mass (GeV);Events")
		for ii in range(1,4):
			AllB.Add(histograms_pass_b[ii])
			AllS.Add(histograms_pass_s[ii])	
			AllSs.Add(histograms_pass_s[ii])
		AllSs.Add(histograms_pass_s[4])	

		Ss.append(AllSs)
		B.append(AllB)

		Chi2B = data.Chi2Test(AllB, "UU CHI2/NDF")
		Chi2S = data.Chi2Test(AllSs, "UU CHI2/NDF")
		print Chi2B
		print Chi2S

		L = TLegend(0.6,0.6,0.89,0.89)
		L.SetFillColor(0)
		L.SetLineColor(0)
		L.AddEntry(data, "1/15^{th} 2016 Dataset", "PL")
		L.AddEntry(AllB, "B Fit (#chi^{2}/ndof = "+str(Chi2B)+")", "L")
		L.AddEntry(AllS, "B+S Fit (#chi^{2}/ndof = "+str(Chi2S)+")", "L")
		L.AddEntry(AllSs, "50 GeV Signal", "F")

		M = FindAndSetMax(data, AllB, AllSs)

		C = TCanvas()
		C.cd()
		AllSs.Draw("hist")
		AllS.Draw("histsame")
		data.Draw("e0same")
		AllB.Draw("histsame")
		L.Draw("same")
		C.Print("IndivBin"+str(i)+"Checker.png")

	Data = D[0].Clone("Data")
	BKG = B[0].Clone("BKG")
	SIG = Ss[0].Clone("SIG")
	for iii in range(1,5):
			Data.Add(D[iii])
			SIG.Add(Ss[iii])	
			BKG.Add(B[iii])
	print " ===== "
	ChiB = Data.Chi2Test(BKG, "UW CHI2")
	ChiS = Data.Chi2Test(SIG, "UW CHI2")
	print ChiB
	print ChiS
	ChiB = Data.Chi2Test(BKG, "UW CHI2/NDF")
	ChiS = Data.Chi2Test(SIG, "UW CHI2/NDF")
	print ChiB
	print ChiS
		


if __name__ == '__main__':
	print "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
	import optparse
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option('--F', '--folder', metavar='folder', type='string', dest='folder')
	(OPTIONS, args) = parser.parse_args()
	print "    " + OPTIONS.folder + "   <<< Using this file!"
	MakeBSCOMP("OUTPUT_"+OPTIONS.folder)
	