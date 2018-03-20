#!/usr/bin/env python

import ROOT as r,sys,math,array,os
from optparse import OptionParser
import random

import tdrstyle
tdrstyle.setTDRStyle()
import limit
from limit import *

if __name__ == "__main__":
	options = parser()
	print options

	bias('cards_33all.txt','cards_33all.txt',options.toys,50.,"fitbase50")
