import ROOT
from ROOT import *
import sys
import math
from array import array
import functools
from multiprocessing import Pool
import scipy

mass=[10.,50.,75.,100.,125.]
xs =  [2.504, 2.127, 2.019, 1.914, 1.918]

G = TGraph(5, scipy.array(mass), scipy.array(xs))

SPLINE = TSpline3("S",G)



bigmass = []
bigXS = []

for x in  [10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.,105.,110.,115.,120.,125.]:
	XS = SPLINE.Eval(x)
	print "   ( "+ str(int(x)) + ", " + str(XS)+" )"
	bigmass.append(x)
	bigXS.append(XS)

G2 = TGraph(len(bigmass), scipy.array(bigmass), scipy.array(bigXS))

G2.SetMarkerStyle(20)
G.SetMarkerStyle(21)
G.SetMarkerColor(kRed)
G2.SetMarkerColor(kViolet)
G2.Draw("AP")
G.Draw("P")

print bigXS