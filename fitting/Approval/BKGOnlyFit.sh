echo "------------------------"
mkdir -p OUTPUT_$1$4_p$2r$3
python buildSRRalphabet.py -b --np $2 --nr $3 --input $1 $4 #--dic OUTPUT_$1$4_p$2r$3/
python makeSRCards.py --np $2 --nr $3
cd cards
combineCards.py card_rhalphabet_cat1.txt card_rhalphabet_cat2.txt card_rhalphabet_cat3.txt card_rhalphabet_cat4.txt card_rhalphabet_cat5.txt > ../cards_$1$4_p$2r$3.txt
cd ..
combine -M MaxLikelihoodFit cards_$1$4_p$2r$3.txt --saveWithUncertainties --saveShapes -v 2 --rMin 0.00001 --rMax 20
mv mlfit.root OUTPUT_$1$4_p$2r$3/MLF.root
mv base.root OUTPUT_$1$4_p$2r$3/base.root
mv ralphabase.root OUTPUT_$1$4_p$2r$3/ralphabase.root
mv cards_$1$4_p$2r$3.txt OUTPUT_$1$4_p$2r$3/card.txt
python ResultsPlotter.py --folder OUTPUT_$1$4_p$2r$3

