#
import ROOT
from ROOT import *
from optparse import OptionParser
import random
import scipy
import sys
import glob, os

for file in glob.glob("PreProc/DATASETS/Z*.root"):
	print file
	rooF = file.split("/")[2]
	print rooF
	tag = rooF.split(".")[0]
	print tag
	#os.system("./BKGOnlyFit.sh " + tag + " 3 3 --pseudo")
	os.system("python ResultsPlotter.py --folder OUTPUT_"+tag+"--pseudo_p3r3/")
