import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import sys



def MakeSuperPoly(X,Y):

        P = ""
        numbering = 0
        for x in range(X+1):
                ypart = "["+str(numbering)+"]"
                numbering += 1
                for y in range(1,Y+1):
                        ypart += " + ["+str(numbering)+"]"
                        for i in range(y):
                                ypart += "*y"
                        numbering += 1
                if x == 0:
                        string = ypart
                else:
                        string = " + ("+ypart+")"
                for i in range(x):
                        string += "*x"
                P += string


        print P

        G = TF2("ralpha"+str(X)+str(Y), P, 200.,1000.,-6.5-2.)
        for n in range(numbering):
                G.SetParameter(n, 0.)
        return G # Makes a TF2 in 2d with X,Y order



def createHist(DDT,tag,filename,sf,lumi,mass,isdata=False):

	massbins = 40;
	masslo   = 0;
	masshi   = 200;

#	ptbins = 4;
#	ptlo = 200;
#	pthi = 800;

	binBoundaries = [200, 230, 254, 290, 362, 1000]

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

		if (i) % sf != 0: 
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
		if jmsd_8 <= 0: jmsd_8 = 0.01

		RHO = math.log(jmsd_8*jmsd_8/PT/PT)

		if RHO < -6.0 or RHO > -2.0: 
			continue

		if not (tree.AK8Puppijet0_N2sdb1 > 0. and PT > 200. and PT < 800.):
			continue
                if isinstance(DDT, TH2F):
                        cur_rho_index = DDT.GetXaxis().FindBin(RHO);
                        cur_pt_index  = DDT.GetYaxis().FindBin(PT);
                        N2DDT = tree.AK8Puppijet0_N2sdb1 - DDT.GetBinContent(cur_rho_index,cur_pt_index);
                elif isinstance(DDT, TF2): 
			N2DDT = tree.AK8Puppijet0_N2sdb1 - DDT.Eval(RHO, PT)
                else: print "THIS ISN'T AN ACCEPTABLE DDT OBJECT"

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
		if (tree.neleLoose == 0 and tree.nmuLoose == 0 and tree.ntau==0 and tree.AK8Puppijet0_isTightVJet ==1):

			# pass category		
			if N2DDT < 0:
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
			if N2DDT > 0:
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

mass=[10,25,50,75,100,125, 200]
xs = [0.6927, 0.5254, 0.3312, 0.2868, 0.2473, 0.2345, 0.1859]
N = [101413, 80878, 103297,99824,96104,96628,193077]

outfile=TFile("ZGammaTemplatesSmooth.root", "recreate");

#lumi =36.600
lumi = 2.44
SF_tau21 =0.88

f_h2ddt = TFile("SmoothDDT.root");
trans_h2ddt = f_h2ddt.Get("N2DDT");

G = MakeSuperPoly(2,1)
G.SetParameter(0, 0.2)
G.SetLineColor(1)
trans_h2ddt.Fit("ralpha21")

data_hists = createHist(G,'data_obs','SinglePhotonmvaEV12av3',15,1,0,True)
qcd_hists = createHist(G,'qcd','NonResBkg',1,lumi,0)
tqq_hists = createHist(G,'tqq','TTmvaEVv3_1000pb_weighted',1,lumi,0)
wqq_hists = createHist(G,'wqq','VectorDiJet1Gammamva_M75',1,lumi,75.)
zqq_hists = createHist(G,'zqq','VectorDiJet1Gammamva_M100',1,lumi,100.)


for m in range(len(mass)):
	hs_hists = createHist(G,'zqq%s'%(mass[m]),'VectorDiJet1Gammamva_M%s'%(mass[m]),1,lumi*1000.,m)
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

