import ROOT
from ROOT import *
import sys
import math
import scipy
from array import array
import functools
from multiprocessing import Pool
import scipy

x = [10.,25.,50.,75.,100.,125.]
y = [-0.2923, 0.2045, 0.2031, 0.1395, 0.1467, 0.2046]
G = TGraph(len(x), scipy.array(x), scipy.array(y))
G.SetLineColor(kBlue)
G.SetLineWidth(3)

Blank = TH2F("B", ";Z' Mass (GeV);Bias      ", 15, 0, 150, 6, -1.5, 1.5)
B.SetStats(0)
B = TBox(0,-0.5,150,0.5)
B.SetFillColor(41)
B.SetLineColor(0)
L = TLine(0,0,150,0)
L.SetLineColor(1)
L.SetLineStyle(2)
L.SetLineWidth(2)
C = TCanvas("C", "", 750, 750)
C.cd()
Blank.Draw("col")
B.Draw("same")
L.Draw("same")
G.Draw("Lsame")
