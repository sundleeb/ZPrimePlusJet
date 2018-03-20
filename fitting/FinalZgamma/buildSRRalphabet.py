#!/usr/bin/env python

import ROOT as r,sys,math,array,os
from ROOT import *
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array
r.gROOT.Macro(os.path.expanduser('~/rootlogon.C'))
# including other directories
sys.path.insert(0, '../.')
from tools import *
from hist import *

gErrorIgnoreLevel = kFatal

global ParValDic
ParValDic = {}
def AssignRooRealVarVals(name, Start, Min, Max):
	print "checking dictionary for " + name
	print " otherwise will use: " + str(Start) + "+/- ( "+str(Min)+", "+str(Max)+")"
	if name in ParValDic:
		print "We will use: "
		S = ParValDic[name]
		E = ParValDic[name+"_err"]
		print name + " =  " + str(S) + " +/-" + str(E)
		R = r.RooRealVar(name,name, S, S-(2*E), S+(2*E))
	else:
		print "using defaults"
		R = r.RooRealVar(name,name,Start,Min,Max)
	return R



##############################################################################
##############################################################################
#### B E G I N N I N G   O F   C L A S S
##############################################################################
##############################################################################

class dazsleRhalphabetBuilder: 

	def __init__( self, hpass, hfail, inputfile, in_np, in_nr ): 

		self.InterpoMasses = []

		self.DDTCUT = 0.1

		self._hpass = hpass;
		self._hfail = hfail;
		self._inputfile = inputfile;

		self._outputName = "base.root";
		self._outfile_validation = r.TFile("validation.root","RECREATE");

		self._mass_nbins = 40;
		self._mass_lo    = 0;
		self._mass_hi    = 200;

		print "number of mass bins and lo/hi: ", self._mass_nbins, self._mass_lo, self._mass_hi;

		#polynomial order for fit
		self._poly_lNP = in_np; # 1 means linear, 2 means quadratic
		self._poly_lNR = in_nr;
		# self._poly_lNRP =2;

		self._nptbins = hpass[0].GetYaxis().GetNbins();
		self._pt_lo = hpass[0].GetYaxis().GetBinLowEdge( 1 );
		self._pt_hi = hpass[0].GetYaxis().GetBinUpEdge( self._nptbins );
	
		# define RooRealVars
		print "define RooRealVars"
		self._lMSD    = r.RooRealVar("x","x",self._mass_lo,self._mass_hi)
		self._lMSD.setBins(self._mass_nbins)		
		self._lPt     = r.RooRealVar("pt","pt",self._pt_lo,self._pt_hi);
		self._lPt.setBins(self._nptbins)
		self._lRho    = r.RooFormulaVar("rho","log(x*x/pt/pt)",r.RooArgList(self._lMSD,self._lPt))

		self._lEff    = r.RooRealVar("veff"      ,"veff"      ,0.5 ,0.,1.0)
		#self._lEffQCD = r.RooRealVar("qcdeff"    ,"qcdeff"    ,0.1 ,0.,1.0)
		self._lEffQCD = AssignRooRealVarVals("qcdeff", self.DDTCUT, 0.0, 1.0)
		self._lDM     = r.RooRealVar("dm","dm", 0.,-10,10)
		self._lShift  = r.RooFormulaVar("shift",self._lMSD.GetName()+"-dm",r.RooArgList(self._lMSD,self._lDM)) 

		self._allVars = [];
		self._allShapes = [];
		self._allData = [];
		self._allPars = [];

		self.LoopOverPtBins();
		#[175, 228, 252, 286, 352, 775][200, 228, 252, 288, 360, 800][150, 230, 254, 290, 362, 1000]
	
	def LoopOverPtBins(self):
		print "number of pt bins = ", self._nptbins;
		for ipt in range(1,self._nptbins+1):
		# for ipt in range(1,2):
			print "------- pT bin number ",ipt		
			
			# 1d histograms in each pT bin (in the order... data, w, z, qcd, top, signals)
			hpass_inPtBin = [];
			hfail_inPtBin = [];
			for ih,h in enumerate(self._hpass):
				tmppass_inPtBin = proj("cat",str(ipt),h,self._mass_nbins,self._mass_lo,self._mass_hi)
				hpass_inPtBin.append( tmppass_inPtBin )
			for ih,h in enumerate(self._hfail):
				tmpfail_inPtBin = proj("cat",str(ipt),h,self._mass_nbins,self._mass_lo,self._mass_hi); 
				hfail_inPtBin.append( tmpfail_inPtBin ) 
			
			# make RooDataset, RooPdfs, and histograms
			curptbincenter = self._hpass[0].GetYaxis().GetBinCenter(ipt);
			(pDatas,pPdfs,pHists) = self.workspaceInputs(hpass_inPtBin,hfail_inPtBin,"cat"+str(ipt),curptbincenter)
			#Get approximate pt bin value
			pPt = self._hpass[0].GetYaxis().GetBinLowEdge(ipt)+self._hpass[0].GetYaxis().GetBinWidth(ipt)*0.3;
			
			#Make the ralphabet fit for a specific pt bin
			lParHists = self.makeRhalph([hfail_inPtBin[0],hfail_inPtBin[1],hfail_inPtBin[2],hfail_inPtBin[4]],[hpass_inPtBin[0],hpass_inPtBin[1],hpass_inPtBin[2],hpass_inPtBin[4]],pPt,"cat"+str(ipt))			
			
			# #Get signals and SM backgrounds
			lPHists=[pHists[0],pHists[1],pHists[2],pHists[3]]
			lFHists=[pHists[4],pHists[5],pHists[6],pHists[7]]
			lPHists.extend(self.getSignals(hpass_inPtBin,hfail_inPtBin,"cat"+str(ipt))[0])
			lFHists.extend(self.getSignals(hpass_inPtBin,hfail_inPtBin,"cat"+str(ipt))[1])
			# #Write to file
			self.makeWorkspace(self._outputName,[pDatas[0]],lPHists,self._allVars,"pass_cat"+str(ipt),True)
			self.makeWorkspace(self._outputName,[pDatas[1]],lFHists,self._allVars,"fail_cat"+str(ipt),True)

		self._outfile_validation.Write();
		self._outfile_validation.Close();
		# for ipt in range(1,self._nptbins+1):
		# 	for imass in range(1,self._mass_nbins):
		# 		print "qcd_fail_cat%i_Bin%i flatParam" % (ipt,imass);
			
	def makeRhalph(self,iHs,iHPs,iPt,iCat):
		
		print "---- [makeRhalph]"	

		lName ="qcd";
		lUnity = r.RooConstVar("unity","unity",1.);
		lZero  = r.RooConstVar("lZero","lZero",0.);

		#Fix the pt (top) and the qcd eff
		self._lPt.setVal(iPt)
		self._lEffQCD.setVal(self.DDTCUT)# @ MARC LOOK HERE!! HERE!!!!
		self._lEffQCD.setConstant(False)

		polyArray = []
		self.buildPolynomialArray(polyArray,self._poly_lNR,self._poly_lNP,"r","p",-1.0,1.0)
		print polyArray;

		#Now build the function
		lPassBins = r.RooArgList()
		lFailBins = r.RooArgList()
		
		for i0 in range(1,self._mass_nbins+1):
			self._lMSD.setVal(iHs[0].GetXaxis().GetBinCenter(i0)) 
			lPass = self.buildRooPolyArray(self._lPt.getVal(),self._lRho.getVal(),lUnity,lZero,polyArray)
			pSum = 0
			pRes = 0

			for i1 in range(0,len(iHs)):
				pSum = pSum + iHs[i1].GetBinContent(i0) if i1 == 0 else pSum - iHs[i1].GetBinContent(i0); # subtract W/Z from data
				if i1 > 0 : pRes += iHs[i1].GetBinContent(i0)
			if pSum < 0: pSum = 0

			#5 sigma range + 10 events
			pUnc = math.sqrt(pSum)*10+10
			#pUnc = math.sqrt(pSum)*3+10
			pUnc += pRes

			#Define the failing category
			pFail = r.RooRealVar(lName+"_fail_"+iCat+"_Bin"+str(i0),lName+"_fail_"+iCat+"_Bin"+str(i0),pSum,max(pSum-pUnc,0),max(pSum+pUnc,0))
			#Now define the passing cateogry based on the failing (make sure it can't go negative)
			lArg = r.RooArgList(pFail,lPass,self._lEffQCD)
			pPass = r.RooFormulaVar(lName+"_pass_"+iCat+"_Bin"+str(i0),lName+"_pass_"+iCat+"_Bin"+str(i0),"@0*max(@1,0)*@2",lArg)
			
			print pPass.Print();
			# print pPass.GetName();
			pSumP = 0
			for i1 in range(0,len(iHPs)):
					pSumP = pSumP + iHPs[i1].GetBinContent(i0) if i1 == 0 else pSumP - iHPs[i1].GetBinContent(i0); # subtract W/Z from data 
			if pSumP < 0: pSumP = 0

			#If the number of events in the failing is small remove the bin from being free in the fit
			if pSum < 5 and pSumP < 5:
				pFail = r.RooRealVar(lName+"_pass_"+iCat+"_Bin"+str(i0),lName+"_pass_"+iCat+"_Bin"+str(i0),pSum ,0.,max(pSum,0.1))
				pFail.setConstant(True)
				pPass = r.RooRealVar(lName+"_pass_"+iCat+"_Bin"+str(i0),lName+"_pass_"+iCat+"_Bin"+str(i0),pSumP,0.,max(pSumP,0.1))
				pPass.setConstant(True)

			#Add bins to the array
			lPassBins.add(pPass)
			lFailBins.add(pFail)
			self._allVars.extend([pPass,pFail])
			self._allPars.extend([pPass,pFail])
			# print  pFail.GetName(),"flatParam",lPass#,lPass+"/("+lFail+")*@0"

		lPass  = r.RooParametricHist(lName+"_pass_"+iCat,lName+"_pass_"+iCat,self._lMSD,lPassBins,iHs[0])
		lFail  = r.RooParametricHist(lName+"_fail_"+iCat,lName+"_fail_"+iCat,self._lMSD,lFailBins,iHs[0])
		lNPass = r.RooAddition(lName+"_pass_"+iCat+"_norm",lName+"_pass_"+iCat+"_norm",lPassBins)
		lNFail = r.RooAddition(lName+"_fail_"+iCat+"_norm",lName+"_fail_"+iCat+"_norm",lFailBins)
		self._allShapes.extend([lPass,lFail,lNPass,lNFail])
		
		#Now write the wrokspace with the rooparamhist
		lWPass = r.RooWorkspace("w_pass_"+str(iCat))
		lWFail = r.RooWorkspace("w_fail_"+str(iCat))
		getattr(lWPass,'import')(lPass,r.RooFit.RecycleConflictNodes())
		getattr(lWPass,'import')(lNPass,r.RooFit.RecycleConflictNodes())
		getattr(lWFail,'import')(lFail,r.RooFit.RecycleConflictNodes())
		getattr(lWFail,'import')(lNFail,r.RooFit.RecycleConflictNodes())
		if iCat.find("1") > -1:
			lWPass.writeToFile("ralpha"+self._outputName)
		else:
			lWPass.writeToFile("ralpha"+self._outputName,False)
		lWFail.writeToFile("ralpha"+self._outputName,False)
		return [lPass,lFail]

	def buildRooPolyArray(self,iPt,iRho,iQCD,iZero,iVars):
		
		# print "---- [buildRooPolyArray]"	

		lPt  = r.RooConstVar("Var_Pt_" +str(iPt)+"_"+str(iRho), "Var_Pt_" +str(iPt)+"_"+str(iRho),(iPt))
		lRho = r.RooConstVar("Var_Rho_"+str(iPt)+"_"+str(iRho), "Var_Rho_"+str(iPt)+"_"+str(iRho),(iRho))
		lRhoArray = r.RooArgList()
		lNCount=0
		for pRVar in range(0,self._poly_lNR+1):
			lTmpArray = r.RooArgList()
			for pVar in range(0,self._poly_lNP+1):
				if lNCount == 0: lTmpArray.add(iQCD); # for the very first constant (e.g. p0r0), just set that to 1
				else: lTmpArray.add(iVars[lNCount])
				lNCount=lNCount+1
			pLabel="Var_Pol_Bin_"+str(round(iPt,2))+"_"+str(round(iRho,3))+"_"+str(pRVar)
			pPol = r.RooPolyVar(pLabel,pLabel,lPt,lTmpArray)
			print pPol.Print()
			lRhoArray.add(pPol);
			self._allVars.append(pPol)

		lLabel="Var_RhoPol_Bin_"+str(round(iPt,2))+"_"+str(round(iRho,3))
		lRhoPol = r.RooPolyVar(lLabel,lLabel,lRho,lRhoArray)
		self._allVars.extend([lPt,lRho,lRhoPol])
		return lRhoPol

	def buildPolynomialArray(self, iVars,iNVar0,iNVar1,iLabel0,iLabel1,iXMin0,iXMax0):
		print "---- [buildPolynomialArray] -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
		for i0 in range(iNVar0+1):
			for i1 in range(iNVar1+1):
				pVar = iLabel1+str(i1)+iLabel0+str(i0);
				print " ::::: ", pVar
				if pVar in ParValDic:
					print ParValDic[pVar]
				pXMin = iXMin0
				pXMax = iXMax0
				pVal  = math.pow(10,-min(i1,2))
				#pVal  = math.pow(10,-i1-i0)
				if i1 == 0: pVal  = math.pow(10,-i1-min(int(i0*0.5),1))
				pCent = 0
				pRooVar = AssignRooRealVarVals(pVar, pCent, pXMin*pVal,pXMax*pVal)
				#pRooVar = r.RooRealVar(pVar,pVar,pCent,pXMin*pVal,pXMax*pVal)
				iVars.append(pRooVar)

	def workspaceInputs(self, iHP,iHF,iBin,iPt):
		print "create the MEAT of the RooFit pieces"
		lCats = r.RooCategory("sample","sample") 
		lCats.defineType("pass",1) 
		lCats.defineType("fail",0) 
		lPData = r.RooDataHist("data_obs_pass_"+iBin,"data_obs_pass_"+iBin,r.RooArgList(self._lMSD),iHP[0])
		lFData = r.RooDataHist("data_obs_fail_"+iBin,"data_obs_fail_"+iBin,r.RooArgList(self._lMSD),iHF[0])
		lData  = r.RooDataHist("comb_data_obs","comb_data_obs",r.RooArgList(self._lMSD),r.RooFit.Index(lCats),r.RooFit.Import("pass",lPData),r.RooFit.Import("fail",lFData)) 

		lW    = self.rooTheHistFunc([iHP[1],iHF[1]],"wqq",iBin)
		lZ    = self.rooTheHistFunc([iHP[2],iHF[2]],"zqq",iBin)
		ltop  = self.rooTheHistFunc([iHP[4],iHF[4]],"tqq",iBin)		
		lQCD  = self.rooTheHistFunc([iHP[3],iHF[3]],"qcd",iBin)

		lTotP = r.RooAddPdf("tot_pass"+iBin,"tot_pass"+iBin,r.RooArgList(lQCD[0]))
		lTotF = r.RooAddPdf("tot_fail"+iBin,"tot_fail"+iBin,r.RooArgList(lQCD[1]))
		lEWKP = r.RooAddPdf("ewk_pass"+iBin,"ewk_pass"+iBin,r.RooArgList(lW[2],lZ[2],ltop[2]))
		lEWKF = r.RooAddPdf("ewk_fail"+iBin,"ewk_fail"+iBin,r.RooArgList(lW[3],lZ[3],ltop[3]))
		
		lTot  = r.RooSimultaneous("tot","tot",lCats) 
		lTot.addPdf(lTotP,"pass") 
		lTot.addPdf(lTotF,"fail")     
		self._allData.extend([lPData,lFData])
		self._allShapes.extend([lTotP,lTotF,lEWKP,lEWKF])

		## find out which to make global
		## RooDataHist (data), then RooAbsPdf (qcd,ewk), then RooHistPdf of each electroweak
		return ([lPData,lFData],[lTotP,lTotF,lEWKP,lEWKF],[lW[4],lZ[4],ltop[4],lQCD[4],lW[5],lZ[5],ltop[5],lQCD[5]])

	def rooTheHistFunc(self,iH,iLabel="w",iBin="_cat0"):

		# normalization
		lNTot   = r.RooRealVar (iLabel+"norm"+iBin,iLabel+"norm"+iBin,(iH[0].Integral()+iH[1].Integral()),0.,5.*(iH[0].Integral()+iH[1].Integral()))
		lNPass  = r.RooFormulaVar(iLabel+"fpass"+iBin ,iLabel+"norm"+iBin+"*(veff)"  ,r.RooArgList(lNTot,self._lEff))
		lNFail  = r.RooFormulaVar(iLabel+"fqail"+iBin ,iLabel+"norm"+iBin+"*(1-veff)",r.RooArgList(lNTot,self._lEff))
		# shapes
		lPData  = r.RooDataHist(iLabel+"_pass_"+iBin,iLabel+"_pass_"+iBin, r.RooArgList(self._lMSD),iH[0])
		lMData  = r.RooDataHist(iLabel+"_fail_"+iBin,iLabel+"_fail_"+iBin, r.RooArgList(self._lMSD),iH[1]) 
		lP      = r.RooHistPdf (iLabel+"passh"+iBin,iLabel+"passh"+iBin, r.RooArgList(self._lShift),r.RooArgList(self._lMSD),lPData,0)
		lF      = r.RooHistPdf (iLabel+"failh"+iBin,iLabel+"failh"+iBin, r.RooArgList(self._lShift),r.RooArgList(self._lMSD),lMData,0)
		# extended likelihood from normalization and shape above
		lEP     = r.RooExtendPdf(iLabel+"_passe_" +iBin,iLabel+"pe" +iBin,lP,lNPass)
		lEF     = r.RooExtendPdf(iLabel+"_faile_" +iBin,iLabel+"fe" +iBin,lF,lNFail)
		
		lHist   = [lP,lF,lEP,lEF,lPData,lMData]
		self._allVars.extend([lNTot,lNPass,lNFail])
		self._allShapes.extend(lHist)
		return lHist	

	def getSignals(self,iHP,iHF,iBin):
			
		lPSigs  = []
		lFSigs  = []
		lPHists = [] 
		lFHists = [] 
		lVars=[10,25,50,75,100,125]
		for i0 in range(0,len(lVars)):
			lSig = self.rooTheHistFunc([iHP[i0+5],iHF[i0+5]],"zqq"+str(lVars[i0]),iBin)
			lPSigs.append(lSig[4])
			lFSigs.append(lSig[5])

		return (lPSigs,lFSigs)		

	def makeWorkspace(self,iOutput,iDatas,iFuncs,iVars,iCat="cat0",iShift=True):
		
		lW = r.RooWorkspace("w_"+str(iCat))

		# get the pT bin
		ipt = iCat[-1:];

		sigMassesForInterpolation_10to25 = [];
		shapeForInterpolation_central_10to25 = [];
		shapeForInterpolation_scaleUp_10to25 = [];
		shapeForInterpolation_scaleDn_10to25 = [];
		shapeForInterpolation_smearUp_10to25 = [];
		shapeForInterpolation_smearDn_10to25 = [];

		sigMassesForInterpolation_25to50 = [];
		shapeForInterpolation_central_25to50 = [];
		shapeForInterpolation_scaleUp_25to50 = [];
		shapeForInterpolation_scaleDn_25to50 = [];
		shapeForInterpolation_smearUp_25to50 = [];
		shapeForInterpolation_smearDn_25to50 = [];

		sigMassesForInterpolation_50to75 = [];
		shapeForInterpolation_central_50to75 = [];
		shapeForInterpolation_scaleUp_50to75 = [];
		shapeForInterpolation_scaleDn_50to75 = [];
		shapeForInterpolation_smearUp_50to75 = [];
		shapeForInterpolation_smearDn_50to75 = [];

		sigMassesForInterpolation_75to100 = [];
		shapeForInterpolation_central_75to100 = [];
		shapeForInterpolation_scaleUp_75to100 = [];
		shapeForInterpolation_scaleDn_75to100 = [];
		shapeForInterpolation_smearUp_75to100 = [];
		shapeForInterpolation_smearDn_75to100 = [];

		sigMassesForInterpolation_100to125 = [];
		shapeForInterpolation_central_100to125 = [];
		shapeForInterpolation_scaleUp_100to125 = [];
		shapeForInterpolation_scaleDn_100to125 = [];
		shapeForInterpolation_smearUp_100to125 = [];
		shapeForInterpolation_smearDn_100to125 = [];


		self._outfile_validation.cd();			

		for pFunc in iFuncs:
			
			ptbin = ipt;
			process = pFunc.GetName().split("_")[0];
			cat     = pFunc.GetName().split("_")[1];
			mass    = 0.;

			#### bbb
			hout = [];
			histDict = {}
			if 'tqq' in process or 'wqq' in process or 'zqq' in process: 
				tmph = self._inputfile.Get(process+'_'+cat).Clone(process+'_'+cat)
				tmph_up = self._inputfile.Get(process+'_'+cat).Clone(process+'_'+cat+'_'+'mcstatUp')
				tmph_down = self._inputfile.Get(process+'_'+cat).Clone(process+'_'+cat+'_'+'mcstatDown')
				# tmph.Scale(getSF(process,cat,self._inputfile))
				# tmph_up.Scale(getSF(process,cat,self._inputfile))
				# tmph_down.Scale(getSF(process,cat,self._inputfile))
				tmph_mass = proj('cat',str(ipt),tmph,self._mass_nbins,self._mass_lo,self._mass_hi)      
				tmph_mass_up = proj('cat',str(ipt),tmph_up,self._mass_nbins,self._mass_lo,self._mass_hi)
				tmph_mass_down = proj('cat',str(ipt),tmph_down,self._mass_nbins,self._mass_lo,self._mass_hi)
				for i in range(1,tmph_mass_up.GetNbinsX()+1):
					mcstatup = tmph_mass_up.GetBinContent(i) + tmph_mass_up.GetBinError(i)
					mcstatdown = max(0.,tmph_mass_down.GetBinContent(i) - tmph_mass_down.GetBinError(i))
					tmph_mass_up.SetBinContent(i,mcstatup)
					tmph_mass_down.SetBinContent(i,mcstatdown)                     
				tmph_mass.SetName(pFunc.GetName())                      
				tmph_mass_up.SetName(pFunc.GetName()+'_'+pFunc.GetName().replace('_','')+'mcstatUp')                
				tmph_mass_down.SetName(pFunc.GetName()+'_'+pFunc.GetName().replace('_','')+'mcstatDown')
				histDict[pFunc.GetName()] = tmph_mass
				histDict[pFunc.GetName()+'_'+pFunc.GetName().replace('_','')+'mcstatUp'] = tmph_mass_up
				histDict[pFunc.GetName()+'_'+pFunc.GetName().replace('_','')+'mcstatDown'] = tmph_mass_down
				uncorrelate(histDict,'mcstat')
				for key, myhist in histDict.iteritems():
					if 'mcstat' in key:
						print key
						hout.append(myhist)
				for h in hout:
					tmprdh = r.RooDataHist(h.GetName(),h.GetName(),r.RooArgList(self._lMSD),h)
					getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
					# validation
					self._outfile_validation.cd()
					h.Write()

			if iShift and ("wqq" in process or "zqq" in process):
				print process
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				print "+++++++++++++"
				if process == "wqq": mass = 80.;
				elif process == "zqq": mass = 91.;
				else: mass = float(process[3:])

				print process, mass;			

				####
				# get the matched and unmatched hist
				print "getting the matched and unmatched hists"
				tmph_matched = self._inputfile.Get(process+"_"+cat+"_matched");
				tmph_unmatched = self._inputfile.Get(process+"_"+cat+"_unmatched");
				tmph_mass_matched = proj("cat",str(ipt),tmph_matched,self._mass_nbins,self._mass_lo,self._mass_hi);
				tmph_mass_unmatched = proj("cat",str(ipt),tmph_unmatched,self._mass_nbins,self._mass_lo,self._mass_hi);

				print "start loop over bins"
					
				#####
				# smear/shift the matched
				hist_container = hist( [mass],[tmph_mass_matched] );	
				mass_shift = 1.0;
				mass_shift_unc = 0.15; # This is 5 sigma shift!  Change the card accordingly
				res_shift = 1.10;
				res_shift_unc = 0.1;
				# get new central value
				shift_val = mass - mass*mass_shift;
				tmp_shifted_h = hist_container.shift( tmph_mass_matched, shift_val);
				# get new central value and new smeared value
				smear_val = res_shift - 1.;
				tmp_smeared_h =  hist_container.smear( tmp_shifted_h[0] , smear_val)
				hmatched_new_central = tmp_smeared_h[0];
				if smear_val <= 0.: hmatched_new_central = tmp_smeared_h[1];
				# get shift up/down
				shift_unc = mass*mass_shift*mass_shift_unc;
				hmatchedsys_shift = hist_container.shift( hmatched_new_central, shift_unc);
				# get res up/down
				hmatchedsys_smear = hist_container.smear( hmatched_new_central, res_shift_unc);	
				#####
				# add back the unmatched 
				print "adding back the unmatched events"
				hmatched_new_central.Add(tmph_mass_unmatched);
				hmatchedsys_shift[0].Add(tmph_mass_unmatched);
				hmatchedsys_shift[1].Add(tmph_mass_unmatched);
				hmatchedsys_smear[0].Add(tmph_mass_unmatched);
				hmatchedsys_smear[1].Add(tmph_mass_unmatched);
				hmatched_new_central.SetName(pFunc.GetName());
				hmatchedsys_shift[0].SetName(pFunc.GetName()+"_scaleUp");
				hmatchedsys_shift[1].SetName(pFunc.GetName()+"_scaleDown");
				hmatchedsys_smear[0].SetName(pFunc.GetName()+"_smearUp");
				hmatchedsys_smear[1].SetName(pFunc.GetName()+"_smearDown");
				hout = [hmatched_new_central,hmatchedsys_shift[0],hmatchedsys_shift[1],hmatchedsys_smear[0],hmatchedsys_smear[1]];
	
				print "match the jets to NOT the W/Z"

				if mass > 0 and process != "wqq" and process != "zqq":
					print mass
					if mass < 30 and mass > 5:
						sigMassesForInterpolation_10to25.append(mass);     
						shapeForInterpolation_central_10to25.append(hmatched_new_central) 
						shapeForInterpolation_scaleUp_10to25.append(hmatchedsys_shift[0]) 
						shapeForInterpolation_scaleDn_10to25.append(hmatchedsys_shift[1])  
						shapeForInterpolation_smearUp_10to25.append(hmatchedsys_smear[0])  
						shapeForInterpolation_smearDn_10to25.append(hmatchedsys_smear[1]) 
					if mass < 55 and mass > 20:
						sigMassesForInterpolation_25to50.append(mass);     
						shapeForInterpolation_central_25to50.append(hmatched_new_central) 
						shapeForInterpolation_scaleUp_25to50.append(hmatchedsys_shift[0]) 
						shapeForInterpolation_scaleDn_25to50.append(hmatchedsys_shift[1])  
						shapeForInterpolation_smearUp_25to50.append(hmatchedsys_smear[0])  
						shapeForInterpolation_smearDn_25to50.append(hmatchedsys_smear[1]) 
					if mass < 80 and mass > 45:
						sigMassesForInterpolation_50to75.append(mass);     
						shapeForInterpolation_central_50to75.append(hmatched_new_central) 
						shapeForInterpolation_scaleUp_50to75.append(hmatchedsys_shift[0]) 
						shapeForInterpolation_scaleDn_50to75.append(hmatchedsys_shift[1])  
						shapeForInterpolation_smearUp_50to75.append(hmatchedsys_smear[0])  
						shapeForInterpolation_smearDn_50to75.append(hmatchedsys_smear[1]) 
					if mass < 105 and mass > 70:
						sigMassesForInterpolation_75to100.append(mass);     
						shapeForInterpolation_central_75to100.append(hmatched_new_central) 
						shapeForInterpolation_scaleUp_75to100.append(hmatchedsys_shift[0]) 
						shapeForInterpolation_scaleDn_75to100.append(hmatchedsys_shift[1])  
						shapeForInterpolation_smearUp_75to100.append(hmatchedsys_smear[0])  
						shapeForInterpolation_smearDn_75to100.append(hmatchedsys_smear[1]) 
					if mass < 130 and mass > 95:
						sigMassesForInterpolation_100to125.append(mass);     
						shapeForInterpolation_central_100to125.append(hmatched_new_central) 
						shapeForInterpolation_scaleUp_100to125.append(hmatchedsys_shift[0]) 
						shapeForInterpolation_scaleDn_100to125.append(hmatchedsys_shift[1])  
						shapeForInterpolation_smearUp_100to125.append(hmatchedsys_smear[0])  
						shapeForInterpolation_smearDn_100to125.append(hmatchedsys_smear[1])

				for h in hout:
					print process
					print h.GetName()					
					h.Write();
					tmprdh = RooDataHist(h.GetName(),h.GetName(),r.RooArgList(self._lMSD),h)
					getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
					if h.GetName().find("scale") > -1:
						pName=h.GetName().replace("scale","scalept")
						tmprdh = RooDataHist(pName,pName,r.RooArgList(self._lMSD),h)
						getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())

			else: 
				getattr(lW,'import')(pFunc,r.RooFit.RecycleConflictNodes())
				
		# do the signal interpolation
		print "---------------------------------------------------------------"
		#print len(sigMassesForInterpolation), sigMassesForInterpolation
		print iCat
		morphedHistContainer_central_10to25 = hist(sigMassesForInterpolation_10to25,shapeForInterpolation_central_10to25);
		morphedHistContainer_scaleUp_10to25 = hist(sigMassesForInterpolation_10to25,shapeForInterpolation_scaleUp_10to25);
		morphedHistContainer_scaleDn_10to25 = hist(sigMassesForInterpolation_10to25,shapeForInterpolation_scaleDn_10to25);
		morphedHistContainer_smearUp_10to25 = hist(sigMassesForInterpolation_10to25,shapeForInterpolation_smearUp_10to25);
		morphedHistContainer_smearDn_10to25 = hist(sigMassesForInterpolation_10to25,shapeForInterpolation_smearDn_10to25);

		morphedHistContainer_central_25to50 = hist(sigMassesForInterpolation_25to50,shapeForInterpolation_central_25to50);
		morphedHistContainer_scaleUp_25to50 = hist(sigMassesForInterpolation_25to50,shapeForInterpolation_scaleUp_25to50);
		morphedHistContainer_scaleDn_25to50 = hist(sigMassesForInterpolation_25to50,shapeForInterpolation_scaleDn_25to50);
		morphedHistContainer_smearUp_25to50 = hist(sigMassesForInterpolation_25to50,shapeForInterpolation_smearUp_25to50);
		morphedHistContainer_smearDn_25to50 = hist(sigMassesForInterpolation_25to50,shapeForInterpolation_smearDn_25to50);

		morphedHistContainer_central_50to75 = hist(sigMassesForInterpolation_50to75,shapeForInterpolation_central_50to75);
		morphedHistContainer_scaleUp_50to75 = hist(sigMassesForInterpolation_50to75,shapeForInterpolation_scaleUp_50to75);
		morphedHistContainer_scaleDn_50to75 = hist(sigMassesForInterpolation_50to75,shapeForInterpolation_scaleDn_50to75);
		morphedHistContainer_smearUp_50to75 = hist(sigMassesForInterpolation_50to75,shapeForInterpolation_smearUp_50to75);
		morphedHistContainer_smearDn_50to75 = hist(sigMassesForInterpolation_50to75,shapeForInterpolation_smearDn_50to75);

		morphedHistContainer_central_75to100 = hist(sigMassesForInterpolation_75to100,shapeForInterpolation_central_75to100);
		morphedHistContainer_scaleUp_75to100 = hist(sigMassesForInterpolation_75to100,shapeForInterpolation_scaleUp_75to100);
		morphedHistContainer_scaleDn_75to100 = hist(sigMassesForInterpolation_75to100,shapeForInterpolation_scaleDn_75to100);
		morphedHistContainer_smearUp_75to100 = hist(sigMassesForInterpolation_75to100,shapeForInterpolation_smearUp_75to100);
		morphedHistContainer_smearDn_75to100 = hist(sigMassesForInterpolation_75to100,shapeForInterpolation_smearDn_75to100);


		morphedHistContainer_central_100to125 = hist(sigMassesForInterpolation_100to125,shapeForInterpolation_central_100to125);
		morphedHistContainer_scaleUp_100to125 = hist(sigMassesForInterpolation_100to125,shapeForInterpolation_scaleUp_100to125);
		morphedHistContainer_scaleDn_100to125 = hist(sigMassesForInterpolation_100to125,shapeForInterpolation_scaleDn_100to125);
		morphedHistContainer_smearUp_100to125 = hist(sigMassesForInterpolation_100to125,shapeForInterpolation_smearUp_100to125);
		morphedHistContainer_smearDn_100to125 = hist(sigMassesForInterpolation_100to125,shapeForInterpolation_smearDn_100to125);

		interpolatedMasses = [15.,20.,30.,35.,40.,45.,55.,60.,65.,70.,80.,85.,90.,95.,105.,110.,115.,120.]
		#interpolatedMasses = []

		for m in interpolatedMasses:
			if m < 25.:
				print "    "
				print "-------------------------------------------------------"
				print "       interpolating mass " + str(m)
				htmp_central = morphedHistContainer_central_10to25.morph(m);
				htmp_scaleUp = morphedHistContainer_scaleUp_10to25.morph(m);
				htmp_scaleDn = morphedHistContainer_scaleDn_10to25.morph(m);
				htmp_smearUp = morphedHistContainer_smearUp_10to25.morph(m);
				htmp_smearDn = morphedHistContainer_smearDn_10to25.morph(m);
				htmp_central.SetName("zqq%i_%s" % (int(m),iCat));
				self.InterpoMasses.append([htmp_central, iCat, int(m)])
				htmp_scaleUp.SetName("zqq%i_%s_scaleUp" % (int(m),iCat)); 
				htmp_scaleDn.SetName("zqq%i_%s_scaleDown" % (int(m),iCat));
				htmp_smearUp.SetName("zqq%i_%s_smearUp" % (int(m),iCat));
				htmp_smearDn.SetName("zqq%i_%s_smearDown" % (int(m),iCat));
				hout = [htmp_central,htmp_scaleUp,htmp_scaleDn,htmp_smearUp,htmp_smearDn];
				for h in hout:
					print h.GetName()
					#self.RhoSafetyCheck(h, ipt)
					h.Write();
					tmprdh = RooDataHist(h.GetName(),h.GetName(),r.RooArgList(self._lMSD),h)
					getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
					if h.GetName().find("scale") > -1:
						pName=h.GetName().replace("scale","scalept")
						tmprdh = RooDataHist(pName,pName,r.RooArgList(self._lMSD),h)
						getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
			elif m < 50.:
				print "    "
				print "-------------------------------------------------------"
				print "       interpolating mass " + str(m)
				htmp_central = morphedHistContainer_central_25to50.morph(m);
				htmp_scaleUp = morphedHistContainer_scaleUp_25to50.morph(m);
				htmp_scaleDn = morphedHistContainer_scaleDn_25to50.morph(m);
				htmp_smearUp = morphedHistContainer_smearUp_25to50.morph(m);
				htmp_smearDn = morphedHistContainer_smearDn_25to50.morph(m);
				htmp_central.SetName("zqq%i_%s" % (int(m),iCat));
				self.InterpoMasses.append([htmp_central, iCat, int(m)])
				htmp_scaleUp.SetName("zqq%i_%s_scaleUp" % (int(m),iCat)); 
				htmp_scaleDn.SetName("zqq%i_%s_scaleDown" % (int(m),iCat));
				htmp_smearUp.SetName("zqq%i_%s_smearUp" % (int(m),iCat));
				htmp_smearDn.SetName("zqq%i_%s_smearDown" % (int(m),iCat));
				hout = [htmp_central,htmp_scaleUp,htmp_scaleDn,htmp_smearUp,htmp_smearDn];
				for h in hout:
					print h.GetName()
					#self.RhoSafetyCheck(h, ipt)
					h.Write();
					tmprdh = RooDataHist(h.GetName(),h.GetName(),r.RooArgList(self._lMSD),h)
					getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
					if h.GetName().find("scale") > -1:
						pName=h.GetName().replace("scale","scalept")
						tmprdh = RooDataHist(pName,pName,r.RooArgList(self._lMSD),h)
						getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
			elif m < 75.:
				print "    "
				print "-------------------------------------------------------"
				print "       interpolating mass " + str(m)
				htmp_central = morphedHistContainer_central_50to75.morph(m);
				htmp_scaleUp = morphedHistContainer_scaleUp_50to75.morph(m);
				htmp_scaleDn = morphedHistContainer_scaleDn_50to75.morph(m);
				htmp_smearUp = morphedHistContainer_smearUp_50to75.morph(m);
				htmp_smearDn = morphedHistContainer_smearDn_50to75.morph(m);
				htmp_central.SetName("zqq%i_%s" % (int(m),iCat));
				self.InterpoMasses.append([htmp_central, iCat, int(m)])
				htmp_scaleUp.SetName("zqq%i_%s_scaleUp" % (int(m),iCat)); 
				htmp_scaleDn.SetName("zqq%i_%s_scaleDown" % (int(m),iCat));
				htmp_smearUp.SetName("zqq%i_%s_smearUp" % (int(m),iCat));
				htmp_smearDn.SetName("zqq%i_%s_smearDown" % (int(m),iCat));
				hout = [htmp_central,htmp_scaleUp,htmp_scaleDn,htmp_smearUp,htmp_smearDn];
				for h in hout:
					print h.GetName()
					#self.RhoSafetyCheck(h, ipt)
					h.Write();
					tmprdh = RooDataHist(h.GetName(),h.GetName(),r.RooArgList(self._lMSD),h)
					getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
					if h.GetName().find("scale") > -1:
						pName=h.GetName().replace("scale","scalept")
						tmprdh = RooDataHist(pName,pName,r.RooArgList(self._lMSD),h)
						getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
			elif m < 100.:
				print "    "
				print "-------------------------------------------------------"
				print "       interpolating mass " + str(m)
				htmp_central = morphedHistContainer_central_75to100.morph(m);
				htmp_scaleUp = morphedHistContainer_scaleUp_75to100.morph(m);
				htmp_scaleDn = morphedHistContainer_scaleDn_75to100.morph(m);
				htmp_smearUp = morphedHistContainer_smearUp_75to100.morph(m);
				htmp_smearDn = morphedHistContainer_smearDn_75to100.morph(m);
				htmp_central.SetName("zqq%i_%s" % (int(m),iCat));
				self.InterpoMasses.append([htmp_central, iCat, int(m)])
				htmp_scaleUp.SetName("zqq%i_%s_scaleUp" % (int(m),iCat)); 
				htmp_scaleDn.SetName("zqq%i_%s_scaleDown" % (int(m),iCat));
				htmp_smearUp.SetName("zqq%i_%s_smearUp" % (int(m),iCat));
				htmp_smearDn.SetName("zqq%i_%s_smearDown" % (int(m),iCat));
				hout = [htmp_central,htmp_scaleUp,htmp_scaleDn,htmp_smearUp,htmp_smearDn];
				for h in hout:
					print h.GetName()
					#self.RhoSafetyCheck(h, ipt)
					h.Write();
					tmprdh = RooDataHist(h.GetName(),h.GetName(),r.RooArgList(self._lMSD),h)
					getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
					if h.GetName().find("scale") > -1:
						pName=h.GetName().replace("scale","scalept")
						tmprdh = RooDataHist(pName,pName,r.RooArgList(self._lMSD),h)
						getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())

			else:
				print "    "
				print "-------------------------------------------------------"
				print "       interpolating mass " + str(m)
				htmp_central = morphedHistContainer_central_100to125.morph(m);
				htmp_scaleUp = morphedHistContainer_scaleUp_100to125.morph(m);
				htmp_scaleDn = morphedHistContainer_scaleDn_100to125.morph(m);
				htmp_smearUp = morphedHistContainer_smearUp_100to125.morph(m);
				htmp_smearDn = morphedHistContainer_smearDn_100to125.morph(m);
				htmp_central.SetName("zqq%i_%s" % (int(m),iCat));
				self.InterpoMasses.append([htmp_central, iCat, int(m)])
				htmp_scaleUp.SetName("zqq%i_%s_scaleUp" % (int(m),iCat)); 
				htmp_scaleDn.SetName("zqq%i_%s_scaleDown" % (int(m),iCat));
				htmp_smearUp.SetName("zqq%i_%s_smearUp" % (int(m),iCat));
				htmp_smearDn.SetName("zqq%i_%s_smearDown" % (int(m),iCat));
				hout = [htmp_central,htmp_scaleUp,htmp_scaleDn,htmp_smearUp,htmp_smearDn];
				for h in hout:
					print h.GetName()
					#self.RhoSafetyCheck(h, ipt)
					h.Write();
					tmprdh = RooDataHist(h.GetName(),h.GetName(),r.RooArgList(self._lMSD),h)
					getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())
					if h.GetName().find("scale") > -1:
						pName=h.GetName().replace("scale","scalept")
						tmprdh = RooDataHist(pName,pName,r.RooArgList(self._lMSD),h)
						getattr(lW,'import')(tmprdh, r.RooFit.RecycleConflictNodes())

		for pData in iDatas:
			getattr(lW,'import')(pData,r.RooFit.RecycleConflictNodes())

		if iCat.find("pass_cat1") == -1:
			lW.writeToFile(iOutput,False)
		else:
			lW.writeToFile(iOutput)	
		# lW.writeToFile(iOutput)	

	def GetInterpoMasses(self):
		print "WE MADE SO MANY MASS POINTS:"
		print str(len(self.InterpoMasses))
		Fout = TFile("InterpolatedMassPoints.root", "recreate")
		Fout.cd()
		for i in self.InterpoMasses:
			i[0].Write()
		Fout.Write()
		Fout.Save()
		Fout.Close()

##############################################################################
##############################################################################
#### E N D   O F   C L A S S
##############################################################################
##############################################################################

def main(options,args):
	# Load the input histograms
	# 	- 2D histograms of pass and fail mass,pT distributions
	# 	- for each MC sample and the data
	f  = r.TFile(options.input);
	(hpass,hfail) = loadHistograms(f,options.pseudo);

	# Build the workspacees
	G = dazsleRhalphabetBuilder(hpass,hfail,f, options.np, options.nr);
	G.GetInterpoMasses()

##-------------------------------------------------------------------------------------
def loadHistograms(f1,pseudo):

	hpass = [];
	hfail = [];

	lHP1 = f1.Get("wqq_pass")
	print 'wqq_pass ', lHP1.Integral() 
	lHF1 = f1.Get("wqq_fail")
	print 'wqq_fail ', lHF1.Integral()
	lHP2 = f1.Get("zqq_pass")
	print 'zqq_pass ', lHP2.Integral()
	lHF2 = f1.Get("zqq_fail")
	print 'zqq_fail ', lHF2.Integral()
	lHP3 = f1.Get("qcd_pass")
	print 'qcd_pass ', lHP3.Integral()
	lHF3 = f1.Get("qcd_fail")
	print 'qcd_fail ', lHF3.Integral()
	lHP4 = f1.Get("tqq_pass")
	print 'tqq_pass ', lHP4.Integral()
	lHF4 = f1.Get("tqq_fail")
	scale=[1.0,0.8,0.75,0.7,0.6,0.5,0.5]
	for i0 in range(1,lHF4.GetNbinsX()+1):
		for i1 in (1,lHF4.GetNbinsY()+1):
			lHP4.SetBinContent(i0,i1,lHP4.GetBinContent(i0,i1)*scale[i1])
			lHF4.SetBinContent(i0,i1,lHF4.GetBinContent(i0,i1)*scale[i1])

	
	print 'tqq_fail ', lHF4.Integral()
	print 'total mc pass ', lHP1.Integral()+lHP2.Integral()+lHP3.Integral()+lHP4.Integral()
	print 'total mc fail ', lHF1.Integral()+lHF2.Integral()+lHF3.Integral()+lHF4.Integral()
  
	if pseudo:
		lHP0 = lHP3.Clone("data_obs_pass")
		lHF0 = lHF3.Clone("data_obs_fail")
		lHF0.Add(lHF1)
		lHF0.Add(lHF2)
		lHF0.Add(lHF4)
		lHP0.Add(lHP1)
		lHP0.Add(lHP2)
		lHP0.Add(lHP4)
		print 'pass ', lHP0.Integral()
		print 'fail ', lHF0.Integral()

	else:
		lHP0 = f1.Get("data_obs_pass")
		lHF0 = f1.Get("data_obs_fail")

	hpass.extend([lHP0,lHP1,lHP2])
	hfail.extend([lHF0,lHF1,lHF2])
	hpass.extend([lHP3,lHP4])
	hfail.extend([lHF3,lHF4])

	#signals
	masses=[10,25,50,75,100,125]
	for mass in masses:
		print mass
		hpass.append(f1.Get("zqq"+str(mass)+"_pass"))
		hfail.append(f1.Get("zqq"+str(mass)+"_fail"))

	for lH in (hpass+hfail):
		lH.SetDirectory(0)	
		print lH.GetName(), lH.Integral()

	print "lengths = ", len(hpass), len(hfail)
	print hpass;
	print hfail;
	return (hpass,hfail);
	# return (hfail,hpass);

##-------------------------------------------------------------------------------------
if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
	parser.add_option('--pseudo', action='store_true', dest='pseudo', default =False,help='data = MC', metavar='isData')
	parser.add_option('--freezefitpar', action='store_true', dest='freezefirpar', default =False,help='do not let the polynopmial parameters change from their initialized values', metavar='FFP')
	parser.add_option('--input', dest='input', default = 'SR15_DDT5_DATA',help='directory with data', metavar='idir')
        parser.add_option('--np', dest="np", type=int,default=3, help='degree poly pt')
        parser.add_option('--nr', dest="nr", type=int,default=3, help='degree poly rho')
	parser.add_option('--dic', dest='dic', default = 'NODIC',help='initialize with these values', metavar='idic')

	(options, args) = parser.parse_args()

	options.input = "PreProc/DATASETS/"+options.input+".root"

	if options.dic != 'NODIC':
		os.system('cp ' + options.dic + ' ParValDictionary.py')
		import ParValDictionary
		ParValDic = ParValDictionary.POSTFITVALS
	print ParValDic
	import tdrstyle
	tdrstyle.setTDRStyle()
	r.gStyle.SetPadTopMargin(0.10)
	r.gStyle.SetPadLeftMargin(0.16)
	r.gStyle.SetPadRightMargin(0.10)
	r.gStyle.SetPalette(1)
	r.gStyle.SetPaintTextFormat("1.1f")
	r.gStyle.SetOptFit(0000)
	r.gROOT.SetBatch()
	
	main(options,args)
	if options.dic != 'NODIC':
		os.system('rm ParValDictionary.py')
##-------------------------------------------------------------------------------------
