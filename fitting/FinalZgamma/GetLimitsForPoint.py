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
import csv

F = TFile("higgsCombineTest.Asymptotic.mH120.root")
L = F.Get("limit")
L.Print()

rows = []
with open(sys.argv[1], 'rb') as f:
	reader = csv.reader(f)
	for row in reader:
		rows.append(row)

for i in range (0,6):
	L.GetEntry(i)
	rows[i].append(L.limit)

for r in rows:
	print r

with open(sys.argv[1], 'wb') as f:
	writer = csv.writer(f)
	writer.writerows(rows)
