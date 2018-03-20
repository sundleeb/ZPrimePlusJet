from array import array
from ROOT import *

def FindAndSetMax(*args):
        maximum = 0.0
        for i in args:
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in args:
                j.GetYaxis().SetRangeUser(0,maximum*1.35) # takes a set of histograms, sets their Range to show all of them

def GoodPlotFormat(H, *args):
	try: H.SetStats(0)
	except: print " ------------ [  No stats box found!  ]"
	if args[0] == 'thickline':
		H.SetLineColor(args[1])
		H.SetLineWidth(2)
	if args[0] == 'thinline':
		H.SetLineColor(args[1])
		H.SetLineWidth(1)
	if args[0] == 'fill':
		H.SetLineColor(args[1])
		H.SetFillColor(args[1])
		H.SetFillStyle(args[2])
	if args[0] == 'e0markers':
		H.SetLineColor(args[1])
		H.SetMarkerColor(args[1])
		H.SetMarkerStyle(args[2])
	H.GetXaxis().SetTitleSize(0.04)

def GetPull(D, E, Style):
	if Style == 'ratio':
		P = D.Clone(D.GetName()+"PullClone")
		P.GetXaxis().SetNdivisions(0)
		P.GetYaxis().SetNdivisions(2)
		P.GetYaxis().SetTitle("#frac{Data}{MC}")
		P.GetYaxis().SetLabelSize(85/15*P.GetYaxis().GetLabelSize())
		P.GetYaxis().SetTitleSize(4.2*P.GetYaxis().GetTitleSize())
		P.GetYaxis().SetTitleOffset(0.175)
		P.GetYaxis().SetRangeUser(0,2)
		P.Divide(E)
		P.SetLineWidth(2)
		P.SetLineColor(kGray+3)
		return P