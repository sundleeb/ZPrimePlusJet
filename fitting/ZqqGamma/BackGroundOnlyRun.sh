#!/bin/bash

#python buildRhalphabet_TT.py -b --np 2 --nr 3
#python makeCards.py --np 2 --nr 3 --nmass 16 --ptbins 4
#cd cards
#combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt > ../cards_TTSBall.txt
#cd ..
#combine -M MaxLikelihoodFit cards_TTSBall.txt --saveWithUncertainties --saveShapes -v 2 --rMin -50 --rMax 50
#python ResultsPlotter.py --W fit_b --O outputTTSB --P 4
#python ResultsPlotter.py --W prefit --O outputTTSB --P 4
#mv mlfit.root mlfitTTSB.root
#python PlotTTEffect.py


python buildRhalphabet.py -b --np 3 --nr 3
python makeCards.py --np 3 --nr 3
cd cards
combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt card_rhalphabet_cat5.txt > ../cards_SRall.txt
cd ..
combine -M MaxLikelihoodFit cards_SRall.txt --saveWithUncertainties --saveShapes -v 2 --rMin -50 --rMax 50
python ResultsPlotter.py --W fit_b --O outputSR
python ResultsPlotter.py --W fit_s --O outputBS
python ResultsPlotter.py --W prefit --O outputSR
mv mlfit.root mlfitSR.root
