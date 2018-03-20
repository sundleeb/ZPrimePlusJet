import ROOT
from ROOT import *
import sys
import math
from array import array
import functools
from multiprocessing import Pool
import scipy
sys.path.append('/home/marc/code/python/')
import PlottingFunctions
import RootHelperFunctions

def ComputeDDT(name, point, nPtBins, nRhoBins, H):
	DDT = TH2F(name, "", nRhoBins, -7.0, -2.0, nPtBins, 200, 800)
	DDT.SetStats(0)
	nXb = H.GetXaxis().GetNbins()
	nYb = H.GetYaxis().GetNbins()
	for x in range(nXb):
		for y in range(nYb):
			proj = H.ProjectionZ("H3"+str(x)+str(y),x+1,x+1,y+1,y+1)
			p = array('d', [point])
			q = array('d', [0.0]*len(p))
			proj.GetQuantiles(len(p), q, p)
			DDT.SetBinContent( x+1, y+1, q[0] );
	return DDT


grho = 0
gpT = 0
nRB = 200
nPT = 200
gDRVAL = 1.0 #should 0.4
Vcut = 0.05

class DDTPoint:
	def __init__(self,r,p,n,w):
		self.rho = r
		self.pT = p
		self.N2 = n
		self.weight = w
	def GetSortNearestCrit(self, o):
		Drho = math.fabs(self.rho - o.rho)/5.
		DpT = math.fabs(self.pT - o.pT)/600.
		return (Drho*Drho) + (DpT*DpT)
	def GetSortNearestGlobal(self, pt, rho):
		Drho = math.fabs(self.rho - rho)/5.
		DpT = math.fabs(self.pT - pt)/600.
		return (Drho*Drho) + (DpT*DpT)
def SortByNearest(p1): # Compare to global (recently set) parameter
	return p1.GetSortNearestGlobal(gpT, grho)
def SortByN2(p1):
	return p1.N2
def DoOneRhoBin(rRr):
	global grho
	grho = -7. + 5.*rRr/float(nRB) + (1./float(nRB))
	Gs = []
	print "Filling rho = " + str(grho) + ""
	F = TFile("gjets.root")
	T = F.Get("tree")
	n = T.GetEntries()
	points = []
	point1 = []
	point2 = []
	point3 = []
	for j in range(0, n):# Here is where we loop over all events.
		T.GetEntry(j)
		dRHO = math.fabs(T.rho - grho)
		if dRHO < gDRVAL:
			points.append(DDTPoint(T.rho, T.pT, T.N2, T.weight*T.puW*T.TrigW*T.kfNLO))
			if j % 3 == 0:
				point1.append(DDTPoint(T.rho, T.pT, T.N2, T.weight*T.puW*T.TrigW*T.kfNLO))
			if j % 3 == 1:
				point2.append(DDTPoint(T.rho, T.pT, T.N2, T.weight*T.puW*T.TrigW*T.kfNLO))
			if j % 3 == 2:
				point3.append(DDTPoint(T.rho, T.pT, T.N2, T.weight*T.puW*T.TrigW*T.kfNLO))
	POINTSARRAY = [points, point1, point2, point3]
	for POINTS in POINTSARRAY:
		ptbins = nPT
		G = []
		for ip in range(0,ptbins):
			global gpT
			gpT = 200. + ip*600./float(ptbins) + (300./float(ptbins))
			tW = 0.
			newPoints = []
			for p in POINTS:
				nR = math.fabs(p.rho - grho)/5.
				nP = math.fabs(p.pT - gpT)/600.
				nW = 1./(0.001+math.sqrt(nR*nR + nP*nP))
				newPoints.append(DDTPoint(p.rho, p.pT, p.N2, p.weight*nW))
				tW += p.weight*nW
			sortedNP = sorted(newPoints, key=SortByN2)
			pN2 = -1.
			pW = 0.
			ID = 0.
			for p in range(len(sortedNP)):
				W = sortedNP[p].weight
				pW += W
				if pW/tW < Vcut: continue
				pN2 = sortedNP[p].N2
				ID = p
				break
			pDelta = (pW - Vcut*tW)/sortedNP[ID].weight
			BIN = [grho, gpT, sortedNP[ID].N2*pDelta + sortedNP[max(ID-1,0)].N2*(1.-pDelta)]
			G.append(BIN)
		Gs.append(G)
	return Gs


if __name__ == '__main__':
	RB = 17
	PB = 17

	H3 = TH3F("H3", ";Jet #rho;Jet p_{T} (GeV)", RB, -7.0, -2.0, PB, 200., 800., 750, 0., 0.75)
	H31 = TH3F("H31", ";Jet #rho;Jet p_{T} (GeV)", RB, -7.0, -2.0, PB, 200., 800., 750, 0., 0.75)
	H32 = TH3F("H32", ";Jet #rho;Jet p_{T} (GeV)", RB, -7.0, -2.0, PB, 200., 800., 750, 0., 0.75)
	H33 = TH3F("H33", ";Jet #rho;Jet p_{T} (GeV)", RB, -7.0, -2.0, PB, 200., 800., 750, 0., 0.75)
	H3.SetStats(0)
	H31.SetStats(0)
	H32.SetStats(0)
	H33.SetStats(0)
	F = TFile("gjets.root")
	T = F.Get("tree")
	n = T.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
		T.GetEntry(j)
		H3.Fill(T.rho, T.pT, T.N2, T.weight*T.puW*T.TrigW*T.kfNLO)
		if j % 3 == 0:
			H31.Fill(T.rho, T.pT, T.N2, T.weight*T.puW*T.TrigW*T.kfNLO)
		if j % 3 == 1:
			H32.Fill(T.rho, T.pT, T.N2, T.weight*T.puW*T.TrigW*T.kfNLO)
		if j % 3 == 2:
			H33.Fill(T.rho, T.pT, T.N2, T.weight*T.puW*T.TrigW*T.kfNLO)
	ddt = ComputeDDT("N2DDT", Vcut, PB, RB, H3)
	ddt1 = ComputeDDT("N2DDT1", Vcut, PB, RB, H31)
	ddt2 = ComputeDDT("N2DDT2", Vcut, PB, RB, H32)
	ddt3 = ComputeDDT("N2DDT3", Vcut, PB, RB, H33)


	p = Pool(22)
	SMOOTH = TH2F("SmooDDT", ";Jet #rho;Jet p_{T} (GeV)", nRB, -7.0, -2.0, nPT, 200., 800.)
	SMOOTH.SetStats(0)
	SMOOTH1 = TH2F("SmooDDT1", ";Jet #rho;Jet p_{T} (GeV)", nRB, -7.0, -2.0, nPT, 200., 800.)
	SMOOTH2 = TH2F("SmooDDT2", ";Jet #rho;Jet p_{T} (GeV)", nRB, -7.0, -2.0, nPT, 200., 800.)
	SMOOTH3 = TH2F("SmooDDT3", ";Jet #rho;Jet p_{T} (GeV)", nRB, -7.0, -2.0, nPT, 200., 800.)
	Smoothed = p.map(DoOneRhoBin, range(0,nRB))
	for S in Smoothed: # this is a rhobin
		for s in S[0]: # this is each of the 4 smooths
			print "---"
			print str(s[0]) + " (rho)"
			print str(s[1]) + " (pT)"
			BINX = SMOOTH.GetXaxis().FindBin(s[0])
			BINY = SMOOTH.GetYaxis().FindBin(s[1])
			print "bins " + str(BINX) + ',' + str(BINY)
			SMOOTH.SetBinContent(BINX, BINY, s[2])
		for s in S[1]: # this is each of the 4 smooths
			BINX = SMOOTH1.GetXaxis().FindBin(s[0])
			BINY = SMOOTH1.GetYaxis().FindBin(s[1])
			SMOOTH1.SetBinContent(BINX, BINY, s[2])
		for s in S[2]: # this is each of the 4 smooths
			BINX = SMOOTH2.GetXaxis().FindBin(s[0])
			BINY = SMOOTH2.GetYaxis().FindBin(s[1])
			SMOOTH2.SetBinContent(BINX, BINY, s[2])
		for s in S[3]: # this is each of the 4 smooths
			BINX = SMOOTH3.GetXaxis().FindBin(s[0])
			BINY = SMOOTH3.GetYaxis().FindBin(s[1])
			SMOOTH3.SetBinContent(BINX, BINY, s[2])

	Fout = TFile("DDTmaps.root", "recreate")
	Fout.cd()

	SMOOTH.Write()
	SMOOTH1.Write()
	SMOOTH2.Write()
	SMOOTH3.Write()
	ddt.Write()
	ddt1.Write()
	ddt2.Write()
	ddt3.Write()

	Fout.Write()
	Fout.Save()
	Fout.Close()

