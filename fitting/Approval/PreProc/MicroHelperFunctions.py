import ROOT
from ROOT import *
import sys
import math
import scipy
import functools
from array import array

def FindAndSetMax(*args):
        maximum = 0.0
	if len(args) == 1: args = args[0]
        for i in args:
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in args:
                j.GetYaxis().SetRangeUser(0,maximum*1.35)
	return maximum

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

def quickplot(File, tree, plot, var, Cut, Weight): # Fills  a plot from a file (needs to have a TTree called "tree"...)
        temp = plot.Clone("temp") # Allows to add multiple distributions to the plot
        chain = ROOT.TChain(tree)
        chain.Add(File)
        chain.Draw(var+">>"+"temp", "("+Weight+")*("+Cut+")", "goff") 
        plot.Add(temp)

def getPUPPIweight( puppipt, puppieta):
	genCorr  = 1.
	recoCorr = 1.
	totalWeight = 1.


	PuppiWeightFile = ROOT.TFile.Open("/uscms_data/d3/mkrohn/CMSSW_8_0_12/src/HH2016/Condor_signalSamples_V25/puppiCorr.root","R")
	puppisd_corrGEN	  = PuppiWeightFile.Get("puppiJECcorr_gen")
	puppisd_corrRECO_cen = PuppiWeightFile.Get("puppiJECcorr_reco_0eta1v3")
	puppisd_corrRECO_for = PuppiWeightFile.Get("puppiJECcorr_reco_1v3eta2v5")
	PuppiWeightFile.Close()

	corrGEN = ROOT.TF1("corrGEN","[0]+[1]*pow(x*[2],-[3])",200,3500)
	corrGEN.SetParameter(0,1.00626)
	corrGEN.SetParameter(1, -1.06161)
	corrGEN.SetParameter(2,0.0799900)
	corrGEN.SetParameter(3,1.20454)
 
	corrRECO_cen = ROOT.TF1("corrRECO_cen","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",200,3500)
	corrRECO_cen.SetParameter(0,1.09302)
	corrRECO_cen.SetParameter(1,-0.000150068)
	corrRECO_cen.SetParameter(2,3.44866e-07)
	corrRECO_cen.SetParameter(3,-2.68100e-10)
	corrRECO_cen.SetParameter(4,8.67440e-14)
	corrRECO_cen.SetParameter(5,-1.00114e-17)
 
	corrRECO_for = ROOT.TF1("corrRECO_for","[0]+[1]*x+[2]*pow(x,2)+[3]*pow(x,3)+[4]*pow(x,4)+[5]*pow(x,5)",200,3500)
	corrRECO_for.SetParameter(0,1.27212)
	corrRECO_for.SetParameter(1,-0.000571640)
	corrRECO_for.SetParameter(2,8.37289e-07)
	corrRECO_for.SetParameter(3,-5.20433e-10)
	corrRECO_for.SetParameter(4,1.45375e-13)
	corrRECO_for.SetParameter(5,-1.50389e-17)
	genCorr =  corrGEN.Eval( puppipt )
	if( abs(puppieta)  < 1.3 ):
		recoCorr = corrRECO_cen.Eval( puppipt )
	else:
		recoCorr = corrRECO_for.Eval( puppipt );
#	print(genCorr,recoCorr)
	totalWeight = genCorr*recoCorr
	return totalWeight