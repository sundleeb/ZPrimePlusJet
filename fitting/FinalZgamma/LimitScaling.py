import ROOT
from ROOT import *
from optparse import OptionParser
from operator import add
import math
import sys
import time
from array import array
import os
import csv
import glob
import LimitPlottingFunctions
from LimitPlottingFunctions import *

if False:
	MC15F = GetLimsFromTXT("Limits_wC_SR15_DDT5_DATA--pseudo_p3r3.txt")
	MC10F = GetLimsFromTXT("Limits_wC_SR10_DDT5_DATA--pseudo_p3r3.txt")
	MC5F = GetLimsFromTXT("Limits_wC_SR5_DDT5_DATA--pseudo_p3r3.txt")
	MC1F = GetLimsFromTXT("Limits_wC_SR1_DDT5_DATA--pseudo_p3r3.txt")

	D15F = GetLimsFromTXT("Limits_wC_SR15_DDT5_DATA_p3r3.txt")
	D10F = GetLimsFromTXT("Limits_wC_SR10_DDT5_DATA_p3r3.txt")
	D5F = GetLimsFromTXT("Limits_wC_SR5_DDT5_DATA_p3r3.txt")
	D1F = GetLimsFromTXT("Limits_wC_SR1_DDT5_DATA_p3r3.txt")
if True: 
	MC15F = GetLimsFromTXT("Limits_SR15_DDT5_DATA--pseudo_p3r3.txt")
	MC10F = GetLimsFromTXT("Limits_SR10_DDT5_DATA--pseudo_p3r3.txt")
	MC5F = GetLimsFromTXT("Limits_SR5_DDT5_DATA--pseudo_p3r3.txt")
	MC1F = GetLimsFromTXT("Limits_SR1_DDT5_DATA--pseudo_p3r3.txt")

	D15F = GetLimsFromTXT("Limits_SR15_DDT5_DATA_p3r3.txt")
	D10F = GetLimsFromTXT("Limits_SR10_DDT5_DATA_p3r3.txt")
	D5F = GetLimsFromTXT("Limits_SR5_DDT5_DATA_p3r3.txt")
	D1F = GetLimsFromTXT("Limits_SR1_DDT5_DATA_p3r3.txt")

print "MC 15"
MCE15 = Make_gQ(MC15F[2])
print "MC 10"
MCE10 = Make_gQ(MC10F[2])
print "MC 5"
MCE5 = Make_gQ(MC5F[2])
print "MC 1"
MCE1 = Make_gQ(MC1F[2])
print "Data 15"
DaE15 = Make_gQ(D15F[2])
print "Data 10"
DaE10 = Make_gQ(D10F[2])
print "Data 5"
DaE5 = Make_gQ(D5F[2])
print "Data 1"
DaE1 = Make_gQ(D1F[2])


x = [100./10., 100./5., 100./1.]
ME10 = []
ME25 = []
ME50 = []
ME75 = []
ME100 = []
ME125 = []
DE10 = []
DE25 = []
DE50 = []
DE75 = []
DE100 = []
DE125 = []
for i in [MCE10, MCE5, MCE1]:
	ME10.append(i[0]/MCE1[0])
	ME25.append(i[1]/MCE1[1])
	ME50.append(i[2]/MCE1[2])
	ME75.append(i[3]/MCE1[3])
	ME100.append(i[4]/MCE1[4])
	ME125.append(i[5]/MCE1[5])
for i in [DaE10, DaE5, DaE1]:
	DE10.append(i[0]/DaE1[0])
	DE25.append(i[1]/DaE1[1])
	DE50.append(i[2]/DaE1[2])
	DE75.append(i[3]/DaE1[3])
	DE100.append(i[4]/DaE1[4])
	DE125.append(i[5]/DaE1[5])

Gm10 = makeAGraph(x, ME10, kRed, 2, 24)
Gm25 = makeAGraph(x, ME25, kBlue, 2, 24)
Gm50 = makeAGraph(x, ME50, kBlue, 2, 24)
Gm75 = makeAGraph(x, ME75, kBlue, 2, 24)
Gm100 = makeAGraph(x, ME100, kRed, 2, 24)
Gm125 = makeAGraph(x, ME125, kGreen, 2, 24)
Gd10 = makeAGraph(x, DE10, kRed, 1, 20)
Gd25 = makeAGraph(x, DE25, kBlue, 1, 20)
Gd50 = makeAGraph(x, DE50, kBlue, 1, 20)
Gd75 = makeAGraph(x, DE75, kBlue, 1, 20)
Gd100 = makeAGraph(x, DE100, kRed, 1, 20)
Gd125 = makeAGraph(x, DE125, kGreen, 1, 20)

JLumi = TF1("x4ths", "100**(1/4)*(1/x**(1/4))", 0.01, 100.)
JLumi.SetLineColor(kBlack)
JLumi.SetLineStyle(2)

L = TLegend(0.6,0.6,0.89,0.89)
L.SetFillColor(0)
L.SetLineColor(0)
L.AddEntry(Gm10, "10 GeV Signal, MC Exp.", "PL")
L.AddEntry(Gm50, "50 GeV Signal, MC Exp.", "PL")
#L.AddEntry(Gm100, "125 GeV Signal, MC Exp.", "PL")
L.AddEntry(Gd10, "10 GeV Signal, Data Exp.", "PL")
L.AddEntry(Gd50, "50 GeV Signal, Data Exp.", "PL")
#L.AddEntry(Gd100, "125 GeV Signal, Data Exp.", "PL")
L.AddEntry(JLumi, "expected L^{1/4} scaling", "L")

ax = TH2F("ax", ";% of dataset used;limit_{partial}/limit_{full}", 91, 9, 100, 9, 1., 4.)
ax.SetStats(0)
C = TCanvas()
C.cd()
#C.SetLogy()
C.SetLogx()
ax.Draw()
Gm10.Draw("LPsame")
Gm50.Draw("LPsame")
#Gm125.Draw("LPsame")
Gd10.Draw("LPsame")
Gd50.Draw("LPsame")
#Gd125.Draw("LPsame")
L.Draw()
JLumi.Draw("same")
C.Print("LimitScalingDemo.png")


