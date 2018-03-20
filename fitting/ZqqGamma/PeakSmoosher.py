import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys

def Format(H,c ):
	H.SetStats(0)
	H.Sumw2()
	H.SetLineColor(c)
	H.SetFillColor(0)
	H.SetMarkerColor(c)
	H.SetMarkerStyle(4)


def FindAndSetMax(someset):
        maximum = 0.0
        for i in someset:
                i.SetStats(0)
                t = i.GetMaximum()
                if t > maximum:
                        maximum = t
        for j in someset:
                j.GetYaxis().SetRangeUser(0,maximum*1.35)

data = True
if not data: F = ["GJetsmvaEVv3.root", "QCDmvaEVwithExtv3.root"]
else: F = ["SinglePhotonmvaEV12av3.root"]


f_h2ddt76 = TFile("PhotonDDTs76.root")
trans_h2ddt76 = f_h2ddt76.Get("N2DDT")
trans_t2ddt76 = f_h2ddt76.Get("Tau21DDT")
f_h2ddt1412 = TFile("PhotonDDTs1412.root")
trans_h2ddt1412 = f_h2ddt1412.Get("N2DDT")
trans_t2ddt1412 = f_h2ddt1412.Get("Tau21DDT")
f_h2ddt2824 = TFile("PhotonDDTs2824.root")
trans_h2ddt2824 = f_h2ddt2824.Get("N2DDT")
trans_t2ddt2824 = f_h2ddt2824.Get("Tau21DDT")

H_TAU1 = TH1F("hTau1", ";jet #tau 21^{DDT}", 50, -0.2, 0.5)
H_TAU2 = TH1F("hTau2", ";jet #tau 21^{DDT}", 50, -0.2, 0.5)
H_TAU3 = TH1F("hTau3", ";jet #tau 21^{DDT}", 50, -0.2, 0.5)

H_DDT1 = TH1F("hDDT1", ";jet N_{2}^{DDT}", 50, -0.2, 0.5)
H_DDT2 = TH1F("hDDT2", ";jet N_{2}^{DDT}", 50, -0.2, 0.5)
H_DDT3 = TH1F("hDDT3", ";jet N_{2}^{DDT}", 50, -0.2, 0.5)

H_Mn01 = TH1F("Mn01", "two pronged mass (N_{2}^{DDT} cut);jet soft drop mass (GeV)", 40, 0, 200)
H_Mn02 = TH1F("Mn02", "two pronged mass (N_{2}^{DDT} cut);jet soft drop mass (GeV)", 40, 0, 200)
H_Mn03 = TH1F("Mn03", "two pronged mass (N_{2}^{DDT} cut);jet soft drop mass (GeV)", 40, 0, 200)

H_Mn21 = TH1F("Mn21", "two pronged mass (#tau 21^{DDT} cut);jet soft drop mass (GeV)", 40, 0, 200)
H_Mn22 = TH1F("Mn22", "two pronged mass (#tau 21^{DDT} cut);jet soft drop mass (GeV)", 40, 0, 200)
H_Mn23 = TH1F("Mn23", "two pronged mass (#tau 21^{DDT} cut);jet soft drop mass (GeV)", 40, 0, 200)

Format(H_TAU1, 1)
Format(H_DDT1, 1)
Format(H_Mn01, 1)
Format(H_Mn21, 1)
Format(H_TAU2, 2)
Format(H_DDT2, 2)
Format(H_Mn02, 2)
Format(H_Mn22, 2)
Format(H_TAU3, 3)
Format(H_DDT3, 3)
Format(H_Mn03, 3)
Format(H_Mn23, 3)

for b in F:
	infile=ROOT.TFile(b)	
	print " "
	tree= infile.Get("Events")
	nent = tree.GetEntries();

	for i in range(tree.GetEntries()):
		if (i) % 10 != 0: continue
		tree.GetEntry(i)

		if tree.neleLoose == 0 and tree.nmuLoose == 0 and tree.ntau==0 and tree.AK8Puppijet0_isTightVJet ==1:
			puweight = tree.puWeight
			fbweight = tree.scale1fb
			weight = puweight*fbweight

			if data: weight = 1.

			PT = tree.AK8Puppijet0_pt
			if not PT > 175.  and PT < 775.:
				continue
			jmsd_8 = tree.AK8Puppijet0_msd

			if not jmsd_8 > 0.:
				continue

			RHO = math.log(jmsd_8*jmsd_8/PT/PT)

			if RHO < -6.0 or RHO > -2.0: 
				continue

			cur_rho_index1 = trans_h2ddt76.GetXaxis().FindBin(RHO);
			cur_pt_index1  = trans_h2ddt76.GetYaxis().FindBin(PT);

			DDT1 = tree.AK8Puppijet0_N2sdb1 - trans_h2ddt76.GetBinContent(cur_rho_index1,cur_pt_index1);
			TauDDT1 = tree.AK8Puppijet0_tau21 - trans_t2ddt76.GetBinContent(cur_rho_index1,cur_pt_index1);


			cur_rho_index2 = trans_h2ddt1412.GetXaxis().FindBin(RHO);
			cur_pt_index2  = trans_h2ddt1412.GetYaxis().FindBin(PT);

			DDT2 = tree.AK8Puppijet0_N2sdb1 - trans_h2ddt1412.GetBinContent(cur_rho_index2,cur_pt_index2);
			TauDDT2 = tree.AK8Puppijet0_tau21 - trans_t2ddt1412.GetBinContent(cur_rho_index2,cur_pt_index2);


			cur_rho_index3 = trans_h2ddt2824.GetXaxis().FindBin(RHO);
			cur_pt_index3  = trans_h2ddt2824.GetYaxis().FindBin(PT);

			DDT3 = tree.AK8Puppijet0_N2sdb1 - trans_h2ddt2824.GetBinContent(cur_rho_index3,cur_pt_index3);
			TauDDT3 = tree.AK8Puppijet0_tau21 - trans_t2ddt2824.GetBinContent(cur_rho_index3,cur_pt_index3);

			
			H_TAU1.Fill(TauDDT1, weight)
			H_DDT1.Fill(DDT1, weight)

			if DDT1 < 0:
				H_Mn01.Fill(jmsd_8, weight)
			if TauDDT1 < 0:
				H_Mn21.Fill(jmsd_8, weight)
			
			H_TAU2.Fill(TauDDT2, weight)
			H_DDT2.Fill(DDT2, weight)

			if DDT2 < 0:
				H_Mn02.Fill(jmsd_8, weight)
			if TauDDT2 < 0:
				H_Mn22.Fill(jmsd_8, weight)
			
			H_TAU3.Fill(TauDDT3, weight)
			H_DDT3.Fill(DDT3, weight)

			if DDT3 < 0:
				H_Mn03.Fill(jmsd_8, weight)
			if TauDDT3 < 0:
				H_Mn23.Fill(jmsd_8, weight)


	print " "

C = TCanvas("C", "", 700, 700)
C.Divide(2,2)
C.cd(2)
H_TAU1.Draw("hist")
H_TAU2.Draw("histsame")
H_TAU3.Draw("histsame")
C.cd(1)
H_DDT1.Draw("histhist")
H_DDT2.Draw("histsame")
H_DDT3.Draw("histsame")
C.cd(3)
H_Mn01.Draw("")
H_Mn02.Draw("same")
H_Mn03.Draw("same")


Div = H_Mn01.Clone("afsd")
r = []
t = 2
for h in [H_Mn01, H_Mn02, H_Mn03]:
	T = h.Clone("blah"+str(t))
	T.Add(Div, -1)
	t+=1
	r.append(T)

FindAndSetMax(r)

C.cd(4)
r[0].Draw("hist")
r[1].Draw("histsame")
r[2].Draw("histsame")
C.Print("Pesmoo_PT.png")
