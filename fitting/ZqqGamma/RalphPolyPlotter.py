import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys

import csv

p0r1  =  9.3384e-01
p0r2  =  9.9498e-02
p0r3  =  6.6309e-02
p1r0  =  -5.3971e-04
p1r1  =  1.3300e-02
p1r2  =  6.7454e-03
p1r3  =  3.4860e-06
p2r0  =  2.9652e-04
p2r1  =  6.7649e-05
p2r2  =  -7.6401e-06
p2r3  =  2.1864e-07
p3r0  =  -5.3531e-07
p3r1  =  -2.0555e-07
p3r2  =  -1.7365e-08
p3r3  =  -1.4571e-09


P = "1.6935e-02*(1 + [0]*y + [1]*y*y + [2]*y*y*y + ([3] + [4]*y + [5]*y*y + [6]*y*y*y)*x + ([7] + [8]*y + [9]*y*y + [10]*y*y*y)*x*x + ([11] + [12]*y + [13]*y*y + [14]*y*y*y)*x*x*x)/0.05"

F = TF2("ralpha", P, 200.,800.,-7.-2.)
F.SetParameter(0, p0r1)
F.SetParameter(1, p0r2)
F.SetParameter(2, p0r3)
F.SetParameter(3, p1r0)
F.SetParameter(4, p1r1)
F.SetParameter(5, p1r2)
F.SetParameter(6, p1r3)
F.SetParameter(7, p2r0)
F.SetParameter(8, p2r1)
F.SetParameter(9, p2r2)
F.SetParameter(10, p2r3)
F.SetParameter(11, p3r0)
F.SetParameter(12, p3r1)
F.SetParameter(13, p3r2)
F.SetParameter(14, p3r3)


F.GetXaxis().SetTitle("AK8 jet p_{T} (GeV)")
F.GetYaxis().SetTitle("AK8 jet #rho")
F.GetXaxis().SetTitleOffset(1.9)
F.GetYaxis().SetTitleOffset(1.9)
F.GetXaxis().SetTitleSize(0.045)
F.GetYaxis().SetTitleSize(0.045)
F.SetTitle("")

F.GetXaxis().SetRangeUser(200.,800.)
F.GetYaxis().SetRangeUser(-7.0,-2.0)
F.GetZaxis().SetRangeUser(0.0,10.0)

F.SetLineColor(1)

C = TCanvas("C", "", 800, 600)
C.cd()
F.Draw("SURF1")
C.Print("PolyOut.pdf")

