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

global TT
TT = False

def createHist(ddt,tag,filename, treename, sf,lumi,mass=0,isdata=False):

	if TT:
		massbins = 28;
		masslo   = 40;
		masshi   = 180;
	else:
		massbins = 30;
		masslo   = 0;
		masshi   = 150;

#	ptbins = 4;
#	ptlo = 175;
#	pthi = 775;

	if TT: binBoundaries = [200, 275, 350, 450, 800]
	else: binBoundaries = [200, 230, 254, 290, 362, 800]

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

	print filename+".root"
	infile=ROOT.TFile(filename+".root")	

	T = infile.Get(treename)
	ne = T.GetEntries()
	for i in range(ne):
		if (i+1) % sf != 0: 
			continue
		T.GetEntry(i)
		rho_index = ddt.GetXaxis().FindBin(T.rho)
		pt_index  = ddt.GetYaxis().FindBin(T.pT)
		DDT = T.N2 - ddt.GetBinContent(rho_index,pt_index)
		dmass = 9
		if mass > 0.:
			dmass = math.fabs(mass - T.SDM)/mass
		matched = (T.V_dphi < 0.8 and T.V_dpt < 0.5 and dmass < 0.3)
		if isdata:
			W = 1.0
		else:
			W = T.weight*T.puW*T.TrigW*T.kfNLO*sf*lumi


		if DDT < 0.:
			h_pass_ak8.Fill( T.SDM, T.pT, W)
			## for signal morphing
			h_pass_msd_ak8.Fill( T.SDM, W);
			if matched:
				h_pass_msd_matched_ak8.Fill( T.SDM, W);
				h_pass_matched_ak8.Fill( T.SDM, T.pT, W);
			else:
				h_pass_msd_unmatched_ak8.Fill( T.SDM, W);
				h_pass_unmatched_ak8.Fill( T.SDM, T.pT, W);
		# fail category
		if DDT > 0.:
			h_fail_ak8.Fill( T.SDM, T.pT, W)
			## for signal morphing
			h_fail_msd_ak8.Fill( T.SDM, W);
			if matched:
				h_fail_msd_matched_ak8.Fill( T.SDM, W);
				h_fail_matched_ak8.Fill( T.SDM, T.pT, W);
			else:
				h_fail_msd_unmatched_ak8.Fill( T.SDM, W);	
				h_fail_unmatched_ak8.Fill( T.SDM, T.pT, W);

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
	print "     -done"
	return hists_out


import optparse
from optparse import OptionParser
parser = OptionParser()

parser.add_option('--T', '--ttbar', metavar='WHC', type='string', dest='ttbar', default="no")
(Options, args) = parser.parse_args()

if Options.ttbar == "yes": TT = True
if TT: print "in ttbar mode"

mass=[10,25,50,75,100,125]
if TT: mass = [50]

xs = [0.6927, 0.5254, 0.3312, 0.2868, 0.2473, 0.2345]
N = [101413, 80878, 103297,99824,96104,96628]
if TT: outfile=TFile("ZGammaFermiTTSBTemplates.root", "recreate");
else: outfile=TFile("ZGammaFermiSRTemplates.root", "recreate");	

print "DDT File loaded"

#lumi =36.600
lumi = 2.44
SF_tau21 =0.88

f_ddt = TFile("DDTmaps.root");
ddt = f_ddt.Get("SmooDDT");
ddt.SetDirectory(0)

ddt.Draw("colz")
if TT: whichtree = 'tttree'
else: whichtree = 'tree'
data_hists = createHist(ddt,'data_obs','data', whichtree, 15,1,0,True)
wqq_hists = createHist(ddt,'wqq','ZQQpG_75',whichtree, 1,lumi,75.)
zqq_hists = createHist(ddt,'zqq','ZQQpG_100',whichtree, 1,lumi,100.)
qcd_hists = createHist(ddt,'qcd','nonres',whichtree, 1,lumi)
tqq_hists = createHist(ddt,'tqq','ttbar',whichtree, 1,lumi)


for m in range(len(mass)):
	hs_hists = createHist(ddt,'zqq%s'%(mass[m]),'ZQQpG_%s'%(mass[m]),whichtree, 1,lumi*1000.,m)
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
