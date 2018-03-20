import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys

def createHist(trans_h2ddt,tag,filename,sf,lumi,mass,isdata=False):

	massbins = 40;
	masslo   = 0;
	masshi   = 200;

#	ptbins = 4;
#	ptlo = 175;
#	pthi = 775;

	binBoundaries = [200, 230, 254, 290, 362, 800]

	h_pass_ak8 = TH2F(tag+"_pass","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_fail_ak8 = TH2F(tag+"_fail","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_pass_matched_ak8 = TH2F(tag+"_pass_matched","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_pass_unmatched_ak8 = TH2F(tag+"_pass_unmatched","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_fail_matched_ak8 = TH2F(tag+"_fail_matched","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_fail_unmatched_ak8 = TH2F(tag+"_fail_unmatched","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	
	# validation
	h_pass_msd_ak8 = TH1F(tag+"pass_msd", "; AK8 m_{SD}^{PUPPI}; N", 40, 0, 200)
	h_fail_msd_ak8 = TH1F(tag+"fail_msd", "; AK8 m_{SD}^{PUPPI}; N", 40, 0, 200)
	h_pass_msd_matched_ak8 = TH1F(tag+"pass_msd_matched", "; AK8 m_{SD}^{PUPPI}; N", 40, 0, 200)
	h_pass_msd_unmatched_ak8 = TH1F(tag+"pass_msd_unmatched", "; AK8 m_{SD}^{PUPPI}; N", 40, 0, 200)
	h_fail_msd_matched_ak8 = TH1F(tag+"fail_msd_matched", "; AK8 m_{SD}^{PUPPI}; N", 40, 0, 200)
	h_fail_msd_unmatched_ak8 = TH1F(tag+"fail_msd_unmatched", "; AK8 m_{SD}^{PUPPI}; N", 40, 0, 200)

	print " "
	print filename+".root"
	infile=ROOT.TFile(filename+".root")	

	tree= infile.Get("Events")
	nent = tree.GetEntries();

	for i in range(tree.GetEntries()):

		if (i+1) % sf != 0: 
			continue
		
		tree.GetEntry(i)
		if(i % (1 * nent/100) == 0):
			sys.stdout.write("\r[" + "="*int(20*i/nent) + " " + str(round(100.*i/nent,0)) + "% done")
			sys.stdout.flush()

		puweight = tree.puWeight
		fbweight = tree.scale1fb * lumi
		weight = puweight*fbweight*sf
		if isdata: weight = puweight*fbweight

		jmsd_8 = tree.AK8Puppijet0_msd
		
		PT = tree.AK8Puppijet0_pt
		if not PT > 0.: 
			continue
		if jmsd_8 <= 3: continue

		RHO = math.log(jmsd_8*jmsd_8/PT/PT)

		if RHO < -7.0 or RHO > -2.0: 
			continue

		jtN2b1sd_8 = tree.AK8Puppijet0_N2sdb1
		cur_rho_index = trans_h2ddt.GetXaxis().FindBin(RHO);
		cur_pt_index  = trans_h2ddt.GetYaxis().FindBin(PT);
		if RHO > trans_h2ddt.GetXaxis().GetBinUpEdge( trans_h2ddt.GetXaxis().GetNbins() ): cur_rho_index = trans_h2ddt.GetXaxis().GetNbins();
		if RHO < trans_h2ddt.GetXaxis().GetBinLowEdge( 1 ): cur_rho_index = 1;
		if PT > trans_h2ddt.GetYaxis().GetBinUpEdge( trans_h2ddt.GetYaxis().GetNbins() ): cur_pt_index = trans_h2ddt.GetYaxis().GetNbins();
		if PT < trans_h2ddt.GetYaxis().GetBinLowEdge( 1 ): cur_pt_index = 1;

		DDT = jtN2b1sd_8 - trans_h2ddt.GetBinContent(cur_rho_index,cur_pt_index);

		# non resonant case
		jphi  = 9999;
		dphi  = 9999;
		dpt   = 9999;
		dmass = 9999;
		if mass > 0:
			jphi = getattr(tree,"AK8Puppijet0_phi");
			dphi = math.fabs(tree.genVPhi - jphi)
			dpt = math.fabs(tree.genVPt - PT)/tree.genVPt
			dmass = math.fabs(mass - jmsd_8)/mass
		
		# Lepton, photon veto and tight jets
		if tree.neleLoose == 0 and tree.nmuLoose == 0 and tree.ntau==0 and tree.pfmet < 100. and tree.vpho0_mva > 0.2 and tree.AK8Puppijet0_msd > 3. and tree.AK8Puppijet0_pt > 200. and tree.AK8Puppijet0_pt < 800 and tree.AK8Puppijet0_N2sdb1 > 0. and tree.AK8Puppijet0_isTightVJet ==1:

			# pass category		
			if tree.AK8Puppijet0_N2sdb1 > 0. and PT > 200. and PT < 800.:

				Jet = TLorentzVector()
				Jet.SetPtEtaPhiM(PT, tree.AK8Puppijet0_eta, tree.AK8Puppijet0_phi, tree.AK8Puppijet0_mass)
				Photon = TLorentzVector()
				Photon.SetPtEtaPhiM(tree.vpho0_pt, tree.vpho0_eta, tree.vpho0_phi, 0.0)
				if Jet.DeltaR(Photon) > 2.2:
					if DDT < 0:
						h_pass_ak8.Fill( jmsd_8, PT, weight )
						## for signal morphing
						h_pass_msd_ak8.Fill( jmsd_8, weight );
						if dphi < 0.8 and dpt < 0.5 and dmass < 0.3:
							h_pass_msd_matched_ak8.Fill( jmsd_8, weight );
							h_pass_matched_ak8.Fill( jmsd_8, PT, weight );
						else:
							h_pass_msd_unmatched_ak8.Fill( jmsd_8, weight );
							h_pass_unmatched_ak8.Fill( jmsd_8, PT, weight );
					# fail category
					if DDT > 0:
						h_fail_ak8.Fill( jmsd_8, PT, weight )
						## for signal morphing
						h_fail_msd_ak8.Fill( jmsd_8, weight );
						if dphi < 0.8 and dpt < 0.5 and dmass < 0.3:
							h_fail_msd_matched_ak8.Fill( jmsd_8, weight );
							h_fail_matched_ak8.Fill( jmsd_8, PT, weight );
						else:
							h_fail_msd_unmatched_ak8.Fill( jmsd_8, weight );	
							h_fail_unmatched_ak8.Fill( jmsd_8, PT, weight );

	hists_out = [];
	#2d histograms
	hists_out.append( h_pass_ak8 );
	hists_out.append( h_fail_ak8 );
	hists_out.append( h_pass_matched_ak8 );
	hists_out.append( h_pass_unmatched_ak8 );
	hists_out.append( h_fail_matched_ak8 );
	hists_out.append( h_fail_unmatched_ak8 );
	#1d validation histograms
	hists_out.append( h_pass_msd_ak8 );
	hists_out.append( h_pass_msd_matched_ak8 );
	hists_out.append( h_pass_msd_unmatched_ak8 );
	hists_out.append( h_fail_msd_ak8 );
	hists_out.append( h_fail_msd_matched_ak8 );
	hists_out.append( h_fail_msd_unmatched_ak8 );

	return hists_out

print "Start"

mass=[10,25,50,75,100,125, 200]
xs = [0.6927, 0.5254, 0.3312, 0.2868, 0.2473, 0.2345, 0.1859]
N = [101413, 80878, 103297,99824,96104,96628,193077]

outfile=TFile("ZGammaTemplates.root", "recreate");

print "DDT File loaded"

#lumi =36.600
lumi = 2.44
SF_tau21 =0.88

f_h2ddt = TFile("SmoothingOutputNewSUPERSmooth.root");
print "hmm"
trans_g2ddt = f_h2ddt.Get("Graph2D");
print "ah!"
trans_h2ddt = trans_g2ddt.GetHistogram()
print "oof"
trans_h2ddt.SetDirectory(0)
print "ouch"
#f_h2ddt.Close()

trans_h2ddt.Draw("colz")

data_hists = createHist(trans_h2ddt,'data_obs','data',15,1,0,True)
qcd_hists = createHist(trans_h2ddt,'qcd','nonres',1,lumi,0)
tqq_hists = createHist(trans_h2ddt,'tqq','ttbar',1,lumi,0)
wqq_hists = createHist(trans_h2ddt,'wqq','ZprimeM75',1,lumi,75.)
zqq_hists = createHist(trans_h2ddt,'zqq','ZprimeM100',1,lumi,100.)


for m in range(len(mass)):
	hs_hists = createHist(trans_h2ddt,'zqq%s'%(mass[m]),'ZprimeM%s'%(mass[m]),1,lumi*1000.,m)
	for h in hs_hists:
		h.Scale(lumi*xs[m]/N[m])
	outfile.cd()
	for h in hs_hists: h.Write();

print("Building pass/fail")	
outfile.cd()
for h in data_hists: h.Write();
for h in qcd_hists: h.Write();
for h in tqq_hists: h.Write();
for h in wqq_hists: 
	h.Scale(0.53/99824*1000.)
	h.Write();
for h in zqq_hists: 
	h.Scale(0.18/96104*1000.)
	h.Write();
outfile.Write()
outfile.Close()

'''

if T.neleLoose == 0 and T.nmuLoose == 0 and T.ntau==0 and T.pfmet < 100. and T.vpho0_mva > 0.2 and T.AK8Puppijet0_msd > 3. and T.AK8Puppijet0_pt > 200. and T.AK8Puppijet0_pt < 800 and T.AK8Puppijet0_N2sdb1 > 0.:
				if T.AK8Puppijet0_pt < (2*T.AK8Puppijet0_msd/0.8): continue
				weight = T.puWeight*T.scale1fb*T.kfactorNLO
				PT = T.AK8Puppijet0_pt
				RHO = 2.*math.log(T.AK8Puppijet0_msd/T.AK8Puppijet0_pt)
				if RHO < -2.0 and RHO > -7.0:
					Jet = TLorentzVector()
					Jet.SetPtEtaPhiM(PT, T.AK8Puppijet0_eta, T.AK8Puppijet0_phi, T.AK8Puppijet0_mass)
					Photon = TLorentzVector()
					Photon.SetPtEtaPhiM(T.vpho0_pt, T.vpho0_eta, T.vpho0_phi, 0.0)
					if Jet.DeltaR(Photon) > 2.2:
						if T.AK8Puppijet0_N2sdb1 > 0.: 
							H3.Fill(RHO, PT, T.AK8Puppijet0_N2sdb1, weight)

'''

