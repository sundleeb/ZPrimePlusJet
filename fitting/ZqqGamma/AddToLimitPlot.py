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

# including other directories
sys.path.insert(0, '../.')

#import tdrstyle
#tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.08);
ROOT.gStyle.SetPadLeftMargin(0.11);
ROOT.gStyle.SetPadTopMargin(0.10);

ROOT.gStyle.SetPalette(1);

def Make_gQ(A):
	out = []
	xs = [0.69, 0.53, 0.33, 0.29, 0.25, 0.23, 0.19]
	#xs = [0.6927, 0.6350711229946524, 0.5789659625668448, 0.5254, 0.4753887165775401, 0.4299475935828877, 0.3900921122994652, 0.3568377540106952, 0.33119999999999994, 0.3136663144385027, 0.30261209304812836, 0.29588471443850267, 0.29133155721925136, 0.2868, 0.2806028256684492, 0.27291443422459893, 0.26437462994652405, 0.25562321711229946, 0.2473, 0.24004478288770054, 0.23449737005347593, 0.23129756577540106, 0.2310851743315508, 0.2345]
	for a in range(len(A)):
		mu = A[a]/xs[a]
		gQ = math.sqrt(mu)/6.
		out.append(gQ)
	print out
	return out

def Upsilonconstraints(x):
	
	mUps = 9.46;
	WidUps = 54e-6;
	alpha = 1/130.;
	Argus = 5.3e-2;
	DimuUps = 2.48e-2;
	RUps = Argus/DimuUps;
	mZp = x[0];

	y = ( 2*4*math.pi*alpha ) * ( math.sqrt( 3.*RUps + 1 ) + 1 ) * (1-(mZp*mZp/(mUps*mUps)))
	if y < 0.: y *= -1;

	return math.sqrt(y);

def Psiconstraints(x):
	
	mUps = 9.46;
	WidUps = 54e-6;
	alpha = 1/130.;
	Argus = 5.3e-2;
	DimuUps = 2.48e-2;
	RUps = Argus/DimuUps;
	mZp = x[0];

	Zwidth = 2.4952;
	ZwidthError = 0.0023*3; # times 3 to give 3 sigma
	relZwidthUnc = ZwidthError/Zwidth;
	sin2thetaW = 0.2312;
	sinthetaW = math.sqrt(sin2thetaW);
	costhetaW = math.sqrt(1 - sin2thetaW);
	mW = 80.385;
	vev = 246.;
	g = mW*2./vev;
	Vu = 0.25 - (4. * sin2thetaW / 6.);
	Vd = -0.25 - (2. * sin2thetaW / 6.);
	mZ = 91.18;
	mZp = x[0];d

	sinthetaW2 = 0.2312;
	sinthetaW = math.sqrt(sinthetaW2);
	costhetaW = math.sqrt(1 - sinthetaW2);
	mW = 80.;
	vev = 246.;
	g = mW*2./vev;
	e = g*sinthetaW;

	Vu = 0.25 - (4. * sinthetaW2 / 6);
	Vd = -0.25 - (2. * sinthetaW2 / 6);
	mZ = 91.18;
	mZp = x[0];

	# y = gZ
	ynum = relZwidthUnc * 3 * g * (1-(mZp*mZp/(mZ*mZ))) * (2*Vu*Vu + 3*Vd*Vd + 5/16.)
	yden = 2*0.01*costhetaW*sinthetaW*(2*Vu+3*Vd);
	# print ynum,yden,x[0],math.sqrt(ynum/yden)
	if ynum > 0: ynum *= -1.;
	y = math.sqrt(ynum/yden);
	
	###
	y = math.sqrt( math.fabs(1-(mZp*mZp/(mZ*mZ))) )

	return y;

def Zconstraints(x):

	mZ = 91.18;
	mZp = x[0];
	ynum = 4. * math.sqrt( 4. * math.pi ) * 1.96 * 1.1e-3 * ( 1-(mZp*mZp/(mZ*mZ)) );
	yden = 1.193 * 0.02;
	if ynum < 0: ynum *= -1.;
	y = math.sqrt(ynum/yden);

	return y;

def makeAGraph(listx,listy,linecolor = 1, linestyle = 1):

	a_m = array('d', []);
	a_g = array('d', []);

	for i in range(len(listx)):
		a_m.append(listx[i]);
		a_g.append(listy[i]);

	gr = ROOT.TGraph(len(listx),a_m,a_g);

	gr.SetLineColor(linecolor)
	gr.SetLineStyle(linestyle)
	gr.SetLineWidth(3)

	return gr

def makeAFillGraph(listx,listy1,listy2,linecolor = 1, fillcolor = 0, fillstyle = 0):
	print listx
	print listy1
	print listy2
	a_m = array('d', []);
	a_g = array('d', []);

	for i in range(len(listx)):
		a_m.append(listx[i]);
		a_g.append(listy1[i]);
	
	for i in range(len(listx)-1,-1,-1):
		a_m.append(listx[i]);
		a_g.append(listy2[i]);

	gr = ROOT.TGraph(2*len(listx),a_m,a_g);

	gr.SetLineColor(linecolor)
	gr.SetFillColor(fillcolor)
	gr.SetFillStyle(fillstyle)

	return gr

def MakeGB(a):
	a[0] = sqrt()

Upsfunc = ROOT.TF1('myfunc',Upsilonconstraints,0,1000,0)
Upsfunc.SetNpx(1000) 
Upsfunc.SetLineColor(30)
Upsfunc.SetLineStyle(9)
Upsfunc.SetLineWidth(2)

Psifunc = ROOT.TF1('myfunc',Psiconstraints,0,1000,0)
Psifunc.SetNpx(1000) 
Psifunc.SetLineColor(kViolet)
Psifunc.SetLineStyle(9)
Psifunc.SetLineWidth(2)

Zfunc = ROOT.TF1('myfunc',Zconstraints,0,1000,0)
Zfunc.SetNpx(1000) 
Zfunc.SetLineColor(15)
Zfunc.SetLineStyle(9)
Zfunc.SetLineWidth(2)

PLOT = TH2F("Blank", "", 50, 8, 150, 20, 0.05, 6.)
PLOT.SetStats(0)
PLOT.GetXaxis().SetTitle("Z' mass (GeV)")
PLOT.GetXaxis().SetTitleSize(0.045)
PLOT.GetYaxis().SetTitle("g_{q'}")
PLOT.GetYaxis().SetTitleSize(0.045)
PLOT.GetYaxis().SetTitleOffset(1.145)
PLOT.GetXaxis().SetTitleOffset(1.185)

#X = [10., 25., 50., 75., 100., 125.]
#MC5 = [2.53,1.69,2.53,3.16,3.25,5.30]
#MC16 = [3.56,2.35,3.48,4.48,4.61,7.65]
#MC50 = [5.29,3.51,5.13,6.78,6.97,11.76]
#MC84 = [8.10,5.47,7.98,10.6,11.00,18.65]
#MC95 = [12.19,8.51,12.34,16.3,17.0,28.9]

#X = [10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125]
X = [10,25,50,75,100,125]
results = []
with open("Limits_SR10_DDT5_DATA.txt") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        results.append(row)
MC5 = results[0]
MC16 = results[1]
MC50 = results[2]
MC84 = results[3]
MC95 = results[4]
OBS = results[5]

g5 = Make_gQ(MC5)
g16 = Make_gQ(MC16)
g50 = Make_gQ(MC50)
g84 = Make_gQ(MC84)
g95 = Make_gQ(MC95)
obs = Make_gQ(OBS)

Big = makeAFillGraph(X,g5,g95,0, kOrange, 1001)
Small = makeAFillGraph(X,g16,g84,0, 3, 1001)
Center = makeAGraph(X,g50, 1, 2)
Obs = makeAGraph(X,obs, 1, 1)

ISRJm = [50,60,70,80,90,100,115]
ISRJgQ = [0.1,0.0575,0.08,0.125,0.0925,0.09,0.14]
ISRJ = makeAGraph(ISRJm, ISRJgQ, kRed, 2)

leg = TLegend(0.3,0.11,0.6,0.3)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(Center, "Expected Limit", "L")
leg.AddEntry(Small, "Expected 1#sigma", "F")
leg.AddEntry(Big, "Expected 2#sigma", "F")
leg.AddEntry(Obs, "Observed", "L")

C = TCanvas("C", "", 800, 600)
C.cd()
C.SetLogy()
C.SetLogx()
PLOT.Draw()
Upsfunc.Draw("same")
Zfunc.Draw("same")
Big.Draw("fsame")
Small.Draw("fsame")
Center.Draw("same")
ISRJ.Draw("sameC")
Obs.Draw("sameC")

leg.Draw()
CMSLABL = TLatex()
CMSLABL.SetNDC()
CMSLABL.SetTextSize(0.045)
PRELABL = TLatex()
PRELABL.SetNDC()
PRELABL.SetTextSize(0.04)
THILABL = TLatex()
THILABL.SetNDC()
THILABL.SetTextSize(0.045)
CMSLABL.DrawLatex(0.1465,0.85,"CMS")
THILABL.DrawLatex(0.81,0.91,"#bf{13 TeV}")
PRELABL.DrawLatex(0.1465,0.812,"#bf{#it{Simulation Preliminary}}")
#C.Print("LimExp.pdf")
C.Print("LimExp.png")



