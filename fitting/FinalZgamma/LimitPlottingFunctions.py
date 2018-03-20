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

def makeAGraph(listx,listy,linecolor = 1, linestyle = 1, mkrstyle = 1):

	a_m = array('d', []);
	a_g = array('d', []);

	for i in range(len(listx)):
		a_m.append(listx[i]);
		a_g.append(listy[i]);

	gr = ROOT.TGraph(len(listx),a_m,a_g);

	gr.SetLineColor(linecolor)
	gr.SetLineStyle(linestyle)
	gr.SetMarkerStyle(mkrstyle)
	gr.SetMarkerColor(linecolor)
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

def GetLimsFromTXT(F):
	results = []
	with open(F) as csvfile:
	    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
	    for row in reader: # each row is a list
		results.append(row)
	return results



