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

	#ftest('cards_F33.txt','cards_F32.txt',options.toys,options.sig,'f32')
	#ftest('cards_F33.txt','cards_F23.txt',options.toys,options.sig,'f23')
	#ftest('cards_F33.txt','cards_F34.txt',options.toys,options.sig,'f34')
	
	#ftest('cards_F33.txt','cards_F43.txt',options.toys,options.sig,'f43')
	bias('cards_SRall.txt','cards_SRall.txt',options.toys, 50.,"fitbase125")

	#goodness('cards_F33.txt',options.toys,"goodness50")
