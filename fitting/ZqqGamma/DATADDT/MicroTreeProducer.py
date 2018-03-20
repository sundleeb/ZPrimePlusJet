import ROOT
from ROOT import *
import sys
import math
import scipy
import functools
from array import array
sys.path.append('/home/marc/code/python/')
import PlottingFunctions
import RootHelperFunctions

class MicroTree:
	def __init__(self, name, File, isData, scaleByMVA):
		self.name = name
		self.isData = isData
		self.scaleByMVA = scaleByMVA
		self.File = File
		self.TrigF = TFile("triggerjetht.root")
		self.TrigH = self.TrigF.Get("PhoPt2_clone")
		self.__book__()
		self.Ev = 0.
		self.TTv = 0.
		self.Wp = 0.
		self.FillMicroTree()
	def GetTrigWeight(self, PT):
		b = self.TrigH.FindFixBin(PT)
		return self.TrigH.GetEfficiency(b)
	def __book__(self):
		self.f = ROOT.TFile( self.name + ".root", "recreate" )
		self.f.cd()
		self.tree = ROOT.TTree("tree", "tree")
		self.tttree = ROOT.TTree("tttree", "tttree")
		# PT
		self.pT = array('f', [-1.0])
		self.addBranch('pT', self.pT, self.tree)
		self.addBranch('pT', self.pT, self.tttree)
		# PhoEta
		self.PhoEta = array('f', [-1.0])
		self.addBranch('PhoEta', self.PhoEta, self.tree)
		self.addBranch('PhoEta', self.PhoEta, self.tttree)
		# Photon PT
		self.phopT = array('f', [-1.0])
		self.addBranch('phopT', self.phopT, self.tree)
		self.addBranch('phopT', self.phopT, self.tttree)
		# JetMass
		self.SDM = array('f', [-1.0])
		self.addBranch('SDM', self.SDM, self.tree)
		self.addBranch('SDM', self.SDM, self.tttree)
		# PrM
		self.PrM = array('f', [-1.0])
		self.addBranch('PrM', self.PrM, self.tree)
		self.addBranch('PrM', self.PrM, self.tttree)
		# rho
		self.rho = array('f', [-1.0])
		self.addBranch('rho', self.rho, self.tree)
		self.addBranch('rho', self.rho, self.tttree)
		#Jg mass
		self.JgM = array('f', [-1.0])
		self.addBranch('JgM', self.JgM, self.tree)
		self.addBranch('JgM', self.JgM, self.tttree)
		#Jet Charge
		self.JC = array('f', [-1.0])
		self.addBranch('JC', self.JC, self.tree)
		self.addBranch('JC', self.JC, self.tttree)
		#N2
		self.N2 = array('f', [-1.0])
		self.addBranch('N2', self.N2, self.tree)
		self.addBranch('N2', self.N2, self.tttree)
		#Deta
		self.Deta = array('f', [-1.0])
		self.addBranch('Deta', self.Deta, self.tree)
		self.addBranch('Deta', self.Deta, self.tttree)
		#V_dphi
		self.V_dphi = array('f', [-1.0])
		self.addBranch('V_dphi', self.V_dphi, self.tree)
		self.addBranch('V_dphi', self.V_dphi, self.tttree)
		#V_dpt
		self.V_dpt = array('f', [-1.0])
		self.addBranch('V_dpt', self.V_dpt, self.tree)
		self.addBranch('V_dpt', self.V_dpt, self.tttree)

		# Weights:
		self.weight = array('f', [0.0])
		self.addBranch('weight', self.weight, self.tree)
		self.addBranch('weight', self.weight, self.tttree)
		self.TrigW = array('f', [0.0])
		self.addBranch('TrigW', self.TrigW, self.tree)
		self.addBranch('TrigW', self.TrigW, self.tttree)
		self.puW = array('f', [0.0])
		self.addBranch('puW', self.puW, self.tree)
		self.addBranch('puW', self.puW, self.tttree)
		self.kf = array('f', [0.0])
		self.addBranch('puW', self.puW, self.tree)
		self.addBranch('puW', self.puW, self.tttree)
		self.kfNLO = array('f', [0.0])
		self.addBranch('kfNLO', self.kfNLO, self.tree)
		self.addBranch('kfNLO', self.kfNLO, self.tttree)

	def FillMicroTree(self):
		File = TFile(self.File)
		print " "
		print " "
		print "---- Starting:" + self.name
		print "Path = " + self.File
		eta4 = 0.
		eta5 = 0.
		T = File.Get("Events")
		self.n = T.GetEntries()
		for j in range(0, self.n): # Here is where we loop over all events.
			if j % 50000 == 0 or j == 1:
				percentDone = float(j) / float(self.n) * 100.0
				print 'Processing '+self.name+' {0:10.0f}/{1:10.0f} : {2:5.2f} %'.format(j, self.n, percentDone )
			T.GetEntry(j)
			self.Wp = self.Wp + (T.scale1fb*T.puWeight)
			if  T.vpho0_mva > 0.2 and T.vpho0_pt > 200 and T.AK8Puppijet0_pt > 200. and T.AK8Puppijet0_msd > 0.:
				PT = T.AK8Puppijet0_pt
				RHO = 2.*math.log(T.AK8Puppijet0_msd/T.AK8Puppijet0_pt)
				if not (RHO < -1.0 and RHO > -9.0): continue
				Jet = TLorentzVector()
				Jet.SetPtEtaPhiM(PT, T.AK8Puppijet0_eta, T.AK8Puppijet0_phi, T.AK8Puppijet0_mass)
				Photon = TLorentzVector()
				Photon.SetPtEtaPhiM(T.vpho0_pt, T.vpho0_eta, T.vpho0_phi, 0.0)
				considercsvs = []
				j1 = TLorentzVector()
				j2 = TLorentzVector()
				j3 = TLorentzVector()
				j4 = TLorentzVector()
				j1.SetPtEtaPhiM(T.AK4Puppijet0_pt, T.AK4Puppijet0_eta, T.AK4Puppijet0_phi, T.AK4Puppijet0_mass)
				j2.SetPtEtaPhiM(T.AK4Puppijet1_pt, T.AK4Puppijet1_eta, T.AK4Puppijet1_phi, T.AK4Puppijet1_mass)
				j3.SetPtEtaPhiM(T.AK4Puppijet2_pt, T.AK4Puppijet2_eta, T.AK4Puppijet2_phi, T.AK4Puppijet2_mass)
				j4.SetPtEtaPhiM(T.AK4Puppijet3_pt, T.AK4Puppijet3_eta, T.AK4Puppijet3_phi, T.AK4Puppijet3_mass)


				if j1.DeltaR(Photon) < 1.5:
					considercsvs.append(T.AK4Puppijet0_csv)
				if j2.DeltaR(Photon) < 1.5:
					considercsvs.append(T.AK4Puppijet1_csv)
				if j3.DeltaR(Photon) < 1.5:
					considercsvs.append(T.AK4Puppijet2_csv)
				if j4.DeltaR(Photon) < 1.5:
					considercsvs.append(T.AK4Puppijet3_csv)

				if len(considercsvs) > 0:
					maxbinevent = max(considercsvs)
				else: maxbinevent = 0.


				closest = 999.
				closept = -1.

				isolated = True
				for J in [j1,j2,j3,j4]:
					if J.Pt() > 5:
						dR = J.DeltaR(Photon)
						if dR < closest:
							closest = dR
							closept = J.Pt()
				if closept > 30. and closest < 0.8:
					isolated = False
				if T.pfmet < 100. and maxbinevent < 0.5426 and T.neleLoose == 0 and isolated and T.AK8Puppijet0_isTightVJet ==1 and Jet.DeltaR(Photon) > 0.0 :
					whichtree = self.tree
					self.Ev = self.Ev + (T.scale1fb*T.puWeight)
				else:
					whichtree = self.tttree
					self.TTv = self.TTv + (T.scale1fb*T.puWeight)

				if self.scaleByMVA: whichtree = self.tree

				self.pT[0] = PT
				self.PhoEta[0] = T.vpho0_eta
				self.phopT[0] = T.vpho0_pt
				self.SDM[0] = T.AK8Puppijet0_msd
				self.PrM[0] = T.AK8Puppijet0_mprun
				self.rho[0] = RHO
				self.JgM[0] = (Jet+Photon).M()
				self.JC[0] = T.AK8Puppijet0_charge
				self.N2[0] = T.AK8Puppijet0_N2sdb1
				self.Deta[0] = T.vpho0_eta - T.AK8Puppijet0_eta
				
				self.weight[0] = T.scale1fb
				self.puW[0] = T.puWeight
				self.TrigW[0] = self.GetTrigWeight(T.vpho0_pt)
				self.kf[0] = T.kfactor
				self.kfNLO[0] = T.kfactorNLO
				if not self.isData:
					self.V_dphi[0] = math.fabs(T.genVPhi - T.AK8Puppijet0_phi)
					self.V_dpt[0] = math.fabs(T.genVPt - PT)/T.genVPt
				else:
					self.V_dphi[0] = -1.
					self.V_dpt[0] = -1.

				whichtree.Fill()

		File.Close()
		print str(self.Ev) + "  ("+str(100*self.Ev/self.Wp)+"%) events in SR"
		print str(self.TTv) + " ("+str(100*self.TTv/self.Wp)+"%) events in SB"

		self.f.cd()
		self.f.Write()
		self.f.Close()

	def addBranch(self, name, var, T): 
		T.Branch(name, var, name+'/F')
	def __del__(self):
		print "done!"



if __name__ == '__main__':

	V10 = MicroTree("ZQQpG_10", "/eos/uscms/store/user/rapte/baconbits/VectorDiJet1Gammamva_M10_noeta.root", False, False)
	V25 = MicroTree("ZQQpG_25", "/eos/uscms/store/user/rapte/baconbits/VectorDiJet1Gammamva_M25_noeta.root", False, False)
	V50 = MicroTree("ZQQpG_50", "/eos/uscms/store/user/rapte/baconbits/VectorDiJet1Gammamva_M50_noeta.root", False, False)
	V75 = MicroTree("ZQQpG_75", "/eos/uscms/store/user/rapte/baconbits/VectorDiJet1Gammamva_M75_noeta.root", False, False)
	V100 = MicroTree("ZQQpG_100", "/eos/uscms/store/user/rapte/baconbits/VectorDiJet1Gammamva_M100_noeta.root", False, True)
	V125 = MicroTree("ZQQpG_125", "/eos/uscms/store/user/rapte/baconbits/VectorDiJet1Gammamva_M125_noeta.root", False, True)


	TT = MicroTree("ttbar", "/eos/uscms/store/user/rapte/baconbits/TTmvaEVv3_noeta.root", False, False)
	QCD = MicroTree("qcd", "/eos/uscms/store/user/rapte/baconbits/QCDmvaEVwithExtv3_noeta.root", False, False)
	GJET = MicroTree("gjets", "/eos/uscms/store/user/rapte/baconbits/GJetsmvaEVv3withExt_noeta.root", False, False)

	Data = MicroTree("data", "/eos/uscms/store/user/rapte/baconbits/SinglePhotonmvaEV12av3_noeta.root", True, False)
