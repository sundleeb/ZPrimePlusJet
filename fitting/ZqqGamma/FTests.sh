#!/bin/bash

python buildRhalphabet.py -b --np 3 --nr 3
python makeCards.py --np 3 --nr 3
cd cards
combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt card_rhalphabet_cat5.txt > ../cards_F33.txt
cd ..

python buildRhalphabet.py -b --np 2 --nr 3
python makeCards.py --np 2 --nr 3
cd cards
combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt card_rhalphabet_cat5.txt > ../cards_F23.txt
cd ..

python buildRhalphabet.py -b --np 3 --nr 2
python makeCards.py --np 3 --nr 2
cd cards
combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt card_rhalphabet_cat5.txt > ../cards_F32.txt
cd ..

python buildRhalphabet.py -b --np 3 --nr 4
python makeCards.py --np 3 --nr 4
cd cards
combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt card_rhalphabet_cat5.txt > ../cards_F34.txt
cd ..

python buildRhalphabet.py -b --np 4 --nr 3
python makeCards.py --np 4 --nr 3
cd cards
combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt card_rhalphabet_cat5.txt > ../cards_F43.txt
cd ..

python FTest.py
ls

