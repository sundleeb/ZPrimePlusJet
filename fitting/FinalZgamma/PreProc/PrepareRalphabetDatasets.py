import ROOT
from ROOT import *
import sys
import math
from array import array
import functools
from multiprocessing import Pool
import scipy

def SignalKeep(S, M):
	if S == 10:
		if M < 15. and M > 5.: return True
	if S == 25:
		if M < 37.5 and M > 15.: return True
	if S == 50:
		if M < 60. and M > 40.: return True
	if S == 75:
		if M < 90. and M > 60.: return True
	if S == 100:
		if M < 120. and M > 75.: return True
	if S == 125:
		if M < 150. and M > 70.: return True
	return False
		

def createHist(which, ddt,tag,filename, treename, sf,lumi,mass=0,isdata=False, issignal=False):
	print "This is region: " + which

	if which == "tt":
		massbins = 40;
		masslo   = 0;
		masshi   = 200;
		binBoundaries = [200, 250, 300, 350, 1000]
		rbinBoundaries = [-7.5,-7.,-6.5,-6.,-5.5,-5.,-4.5,-4.,-3.5,-3.,-2.5,-2.,]
	else:
		massbins = 40;
		masslo   = 0;
		masshi   = 200;
		binBoundaries = [200, 250, 280, 325, 380, 1000]
		rbinBoundaries = [-7.5,-7.,-6.5,-6.,-5.5,-5.,-4.5,-4.,-3.5,-3.,-2.5,-2.,]

	# HISTOGRAMS:
	h_pass_ak8 = TH2F(tag+"_pass","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_fail_ak8 = TH2F(tag+"_fail","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_pass_matched_ak8 = TH2F(tag+"_pass_matched","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_pass_unmatched_ak8 = TH2F(tag+"_pass_unmatched","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_fail_matched_ak8 = TH2F(tag+"_fail_matched","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	h_fail_unmatched_ak8 = TH2F(tag+"_fail_unmatched","; AK8 m_{SD}^{PUPPI} (GeV); AK8 p_{T} (GeV)",massbins,masslo,masshi,len(binBoundaries)-1, array('d',binBoundaries))
	
	h_pass_rho = TH2F(tag+"_pass_rho","; AK8 #rho; AK8 p_{T} (GeV)",len(rbinBoundaries)-1, array('d',rbinBoundaries),len(binBoundaries)-1, array('d',binBoundaries))
	h_fail_rho = TH2F(tag+"_fail_rho","; AK8 #rho; AK8 p_{T} (GeV)",len(rbinBoundaries)-1, array('d',rbinBoundaries),len(binBoundaries)-1, array('d',binBoundaries))

	# validation
	h_pass_msd_ak8 = TH1F(tag+"_pass_msd", "; AK8 m_{SD}^{PUPPI}; N", massbins,masslo,masshi)
	h_fail_msd_ak8 = TH1F(tag+"_fail_msd", "; AK8 m_{SD}^{PUPPI}; N", massbins,masslo,masshi)
	h_pass_msd_matched_ak8 = TH1F(tag+"_pass_msd_matched", "; AK8 m_{SD}^{PUPPI}; N", massbins,masslo,masshi)
	h_pass_msd_unmatched_ak8 = TH1F(tag+"_pass_msd_unmatched", "; AK8 m_{SD}^{PUPPI}; N", massbins,masslo,masshi)
	h_fail_msd_matched_ak8 = TH1F(tag+"_fail_msd_matched", "; AK8 m_{SD}^{PUPPI}; N", massbins,masslo,masshi)
	h_fail_msd_unmatched_ak8 = TH1F(tag+"_fail_msd_unmatched", "; AK8 m_{SD}^{PUPPI}; N", massbins,masslo,masshi)

	infile = TFile(filename+".root")
	T = infile.Get(treename)
	ne = T.GetEntries()
	for i in range(ne):
		if (i+1) % sf != 0: 
			continue
		T.GetEntry(i)
		if not issignal:
			if which == "tt":
				if math.fabs(T.PhoEta) > 2.1 or T.rho > -2.0 or T.rho < -7.5 or T.SR > 0.: continue
			elif which == "sr":
				if math.fabs(T.PhoEta) > 2.1 or T.rho > -2.0 or T.rho < -7.5 or T.SR < 0.: continue
			else:
				if math.fabs(T.PhoEta) > 2.1 or T.rho > -2.0 or T.rho < -7.5: continue
		else:
			#if not SignalKeep(mass, T.SDM*T.TheaW): continue
			if math.fabs(T.PhoEta) > 2.1 or T.rho > -2.0 or T.rho < -7.5: continue
		rho_index = ddt.GetXaxis().FindBin(T.rho)
		pt_index  = ddt.GetYaxis().FindBin(T.pT)
		DDT = T.N2 - ddt.GetBinContent(rho_index,pt_index)
		dmass = 9
		if mass > 0.:
			dmass = math.fabs(mass - T.SDM*T.TheaW)/mass
		matched = (T.V_dphi < 0.8 and T.V_dpt < 0.5 and dmass < 0.3)
		if isdata:
			W = 1.0
		else:
			W = T.weight*T.puW*T.TrigW*T.kfNLO*sf*lumi
		if DDT < 0.:
			h_pass_ak8.Fill( T.SDM*T.TheaW, T.pT, W)
			h_pass_rho.Fill( T.rho, T.pT, W)
			## for signal morphing
			h_pass_msd_ak8.Fill( T.SDM*T.TheaW, W);
			if matched:
				h_pass_msd_matched_ak8.Fill( T.SDM*T.TheaW, W);
				h_pass_matched_ak8.Fill( T.SDM*T.TheaW, T.pT, W);
			else:
				h_pass_msd_unmatched_ak8.Fill( T.SDM*T.TheaW, W);
				h_pass_unmatched_ak8.Fill( T.SDM*T.TheaW, T.pT, W);
		# fail category
		if DDT > 0.:
			h_fail_ak8.Fill( T.SDM*T.TheaW, T.pT, W)
			h_fail_rho.Fill( T.rho, T.pT, W)
			## for signal morphing
			h_fail_msd_ak8.Fill( T.SDM*T.TheaW, W);
			if matched:
				h_fail_msd_matched_ak8.Fill( T.SDM*T.TheaW, W);
				h_fail_matched_ak8.Fill( T.SDM*T.TheaW, T.pT, W);
			else:
				h_fail_msd_unmatched_ak8.Fill( T.SDM*T.TheaW, W);	
				h_fail_unmatched_ak8.Fill( T.SDM*T.TheaW, T.pT, W);

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

	hists_out.append( h_pass_rho);
	hists_out.append( h_fail_rho );

	return hists_out

def FillHists(name, FDDT, FRAC, W, OPT):
	
	outfile=TFile("DATASETS/"+name+".root", "recreate");	

	lumi = 35.9 / FRAC	
	SF_tau21 =0.85

	print "."
	f_ddt = TFile(FDDT[0]);
	f_ddt.Print()
	print "."
	ddt = f_ddt.Get(FDDT[1]);
	ddt.Print()
	ddt.SetDirectory(0)

	data_hists = createHist(W, ddt,'data_obs','data', "tree", FRAC,1,0,True)
	wqq_hists = createHist(W, ddt,'wqq','WQQ',"tree", 1,lumi,80)
	zqq_hists = createHist(W, ddt,'zqq','ZQQ',"tree", 1,lumi,91)
	qcd_hists = createHist(W, ddt,'qcd',OPT,"tree", 1,lumi)
	tqq_hists = createHist(W, ddt,'tqq','ttbar',"tree", 1,lumi)

	mass=[10,25,50,75,100,125]
	xs = [0.6927, 0.5254, 0.3312, 0.2868, 0.2473, 0.2345]
	N = [101413, 80878, 103297,99824,96104,96628]
	MVAScale = [0.9252, 0.90715, 0.9029, 0.8832, 0.87195, 0.8651]
	for m in range(len(mass)):
		hs_hists = createHist(W, ddt,'zqq%s'%(mass[m]),'ZQQpG_%s'%(mass[m]),"tree", 1,lumi,mass[m],False,True)
		for h in hs_hists:
		#	h.Scale(1000.*MVAScale[m]*xs[m]/N[m])
			h.Scale(1000.*xs[m]/N[m])
		outfile.cd()
		for h in hs_hists: h.Write();
	outfile.cd()
	for h in data_hists: h.Write();
	for h in qcd_hists: h.Write();
	for h in tqq_hists: h.Write();
	for h in wqq_hists: 
		h.Scale(1000.*1.22/1507795.125)
		h.Write();
	for h in zqq_hists: 
		h.Scale(1000.*0.569/398015.4276)
		h.Write();
	outfile.Write()
	outfile.Close()

for M in ["DATA"]:#, "MCETA"]:
	for b in ["100"]:
		for w in ["1"]:
			DDT5 = ["DDTs5_"+b+"_"+w+".root", "Smooth_N2DDT_"+M]
			DDT10 = ["DDTs10_"+b+"_"+w+".root", "Smooth_N2DDT_"+M]
			DDT15 = ["DDTs15_"+b+"_"+w+".root", "Smooth_N2DDT_"+M]
			DDT25 = ["DDTs25_"+b+"_"+w+".root", "Smooth_N2DDT_"+M]
		#	FillHists("ZGp5", DDT5, 15,"sr",'gjets')
		#	FillHists("ZGf5", DDT5, 1,"sr",'gjets')
			FillHists("ZtGp10", DDT10, 15,"",'gjets')
			FillHists("ZtGf10", DDT10, 1,"",'gjets')
		#	FillHists("ZGp15", DDT15, 15,"sr",'gjets')
		#	FillHists("ZGf15", DDT15, 1,"sr",'gjets')
			FillHists("ZGp25", DDT25, 15,"sr",'gjets')
			FillHists("ZGf25", DDT25, 1,"sr",'gjets')

#for N in [15,5,1]:
#	for M in ["MCETA", "DATA"]:
#		for V in ["5", "10", "20"]:
#			FillHists("SR"+str(N)+"_DDT"+V+"_"+M, ["DDTs"+V+".root","Smooth_N2DDT_"+M], N, "sr", 'nonres')
#			FillHists("SRG"+str(N)+"_DDT"+V+"_"+M, ["DDTs"+V+".root","Smooth_N2DDT_"+M], N, "sr", 'gjets')
#			FillHists("TCR"+str(N)+"_DDT"+V+"_"+M, ["DDTs"+V+".root","Smooth_N2DDT_"+M], N, "tt", 'nonres')
			#FillHists("BHR"+str(N)+"_DDT"+V+"_"+M, ["DDTs"+V+".root","Smooth_N2DDT_"+M], N, "full")



















