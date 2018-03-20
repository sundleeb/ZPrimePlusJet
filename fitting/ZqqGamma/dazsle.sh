#!/bin/sh

python buildRhalphabet.py -b --np 2 --nr 2
python makeCards.py --np 2 --nr 2
cd cards
combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt card_rhalphabet_cat5.txt > ../cards_22all.txt
cd ..
combine -M MaxLikelihoodFit cards_22all.txt --saveWithUncertainties --saveShapes -v 2 --rMin -50 --rMax 50
