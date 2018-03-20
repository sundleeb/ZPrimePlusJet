#!/usr/bin/env python

import ROOT as r,sys,math,array,os
from optparse import OptionParser
import random
import sys

import CombineFunctions
from CombineFunctions import *

if __name__ == "__main__":
	options = parser()
	print options

	#ftest('cards_F33.txt','cards_F32.txt',options.toys,options.sig,'f32')
	#ftest('cards_F33.txt','cards_F23.txt',options.toys,options.sig,'f23')
	#ftest('cards_F33.txt','cards_F34.txt',options.toys,options.sig,'f34')
	#ftest('cards_F33.txt','cards_F43.txt',options.toys,options.sig,'f43')
	bias('card_BIAS.txt','card_BIAS.txt', 215, 35., sys.argv[1])

	#goodness('cards_F33.txt',options.toys,"goodness50")
