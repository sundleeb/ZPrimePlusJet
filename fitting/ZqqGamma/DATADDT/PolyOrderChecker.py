import ROOT
from ROOT import *
import sys
import math
from array import array
import functools
from multiprocessing import Pool
import scipy
import PlottingFunctions
import RootHelperFunctions

def makePolyCheck(H, name):
	print "making poly ratio: " + name
	hp = H[-2]
	hf = H[-1]
	hr = hp.Clone("hr"+name)
	hr.SetStats(0)
	hr.Divide(hf)


	print "Chi2 = " + str(F.GetChisquare())

	C = TCanvas("C"+name, "", 800, 600)
	C.cd()
	hr.Draw("colz")
	C.Print("PolyToBeFit_"+name+".png")



F = TFile("ZGammaFermiSRTemplates.root")

Sp = F.Get("zqq50_pass_rho")
Sf = F.Get("zqq50_fail_rho")

Hp = F.Get("data_obs_pass_rho")
Htp = F.Get("tqq_pass_rho")
Hp.Add(Htp, -1.)
Hf = F.Get("data_obs_fail_rho")
Htf = F.Get("tqq_fail_rho")
Hf.Add(Htf, -1.)

Hr = Hp.Clone("PolyRatio")
Hr.Divide(Hf)
Hr.SetStats(0)
Hr.GetXaxis().SetTitleOffset(2.)
Hr.GetYaxis().SetTitleOffset(2.25)
#Hr.GetXaxis().SetTitleSize(1.2)
#Hr.GetYaxis().SetTitleSize(1.2)

#P = "(1 + [0]*y + [1]*y*y + [2]*y*y*y + ([3] + [4]*y + [5]*y*y + [6]*y*y*y)*x + ([7] + [8]*y + [9]*y*y + [10]*y*y*y)*x*x + ([11] + [12]*y + [13]*y*y + [14]*y*y*y)*x*x*x)*[15]*0.05"
P = "(1 + [1]*x + [2]*y + [3]*x*y + [4]*x*x + [5]*y*y + [6]*y*y*x + [7]*x*x*y +[8]+x*x*x + [9]*y*y*y)*[0]*0.05"
#P = "(1 + [0]*x + [1]*y + [2]*x*y + [3]*x*x + [4]*y*y + [5]*x*y*y + [6]*x*x*y + [7]*x*x*x + [8]*y*y*y)*[9]*0.05"
Fit = TF2("ralphapoly", P, 200.,800.,-7.-2.)
Hr.Fit(Fit, "")

Ndf = Fit.GetNDF()
Chi2 = Fit.GetChisquare()

#print str(Chi2/Ndf)


C = TCanvas("C", "", 800, 600)
C.cd()
Hr.Draw("SURF1")
Fit.Draw("same")
Text = '#chi^{2}/n.d.f. = '+'{:.4}'.format(Chi2/Ndf)

p1 = TPaveLabel()
p1.SetFillColor(kWhite)
xMin=0.1
yMin=0.85
xMax=0.4
yMax=0.95
p1.DrawPaveLabel(xMin,yMin,xMax,yMax, Text,"brNDC")

C.Print("PolyToBeFit.png")
